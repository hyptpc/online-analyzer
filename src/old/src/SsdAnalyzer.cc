// -*- C++ -*-

#include "SsdAnalyzer.hh"

#include <string>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include "UnpackerManager.hh"

#include "DetectorID.hh"
#include "HodoParamMan.hh"

namespace
{
  const std::string& class_name("SsdAnalyzer");
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const HodoParamMan& gHodoParamMan = HodoParamMan::GetInstance();
}

//______________________________________________________________________________
SsdAnalyzer::SsdAnalyzer( void )
  : m_ssdt(NumOfSegSSDT),
    m_ssdct(NumOfSegSSDT),
    m_adc(NumOfLayersSSD1+NumOfLayersSSD2,
	  std::vector<Double_t>(NumOfSegSSD1, -9999.) ),
    m_tdc(NumOfLayersSSD1+NumOfLayersSSD2,
	  std::vector<Double_t>(NumOfSegSSD1, -9999.) ),
    m_de(NumOfLayersSSD1+NumOfLayersSSD2,
	 std::vector<Double_t>(NumOfSegSSD1, -9999.) ),
    m_time(NumOfLayersSSD1+NumOfLayersSSD2,
	   std::vector<Double_t>(NumOfSegSSD1, -9999.)),
    m_chisqr(NumOfLayersSSD1+NumOfLayersSSD2,
	     std::vector<Double_t>(NumOfSegSSD1, -9999.))
{
}

//______________________________________________________________________________
SsdAnalyzer::~SsdAnalyzer( void )
{
}

//______________________________________________________________________________
Bool_t
SsdAnalyzer::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  // SSDT
  {
    static const Int_t k_device = gUnpacker.get_device_id("SSDT");
    static const Int_t k_tdc    = gUnpacker.get_data_id("SSDT","tdc");
    for(Int_t seg=0; seg<NumOfSegSSDT; ++seg){
      Int_t nhit = gUnpacker.get_entries(k_device, 0, seg*2, 0, k_tdc);
      if(nhit==0) continue;
      Int_t tdc  = gUnpacker.get(k_device, 0, seg*2, 0, k_tdc);
      Double_t time;
      gHodoParamMan.GetTime(DetIdSSDT, 0, seg*2, 0, tdc, time);
      m_ssdt[seg]  = tdc;
      m_ssdct[seg] = time;
    }
#if 0
    gUnpacker.dump_data_device(k_device);
#endif
  }

  // SSD1/SSD2
  {
    static const Int_t k_device[2] = {
      gUnpacker.get_device_id("SSD1"),
      gUnpacker.get_device_id("SSD2")
    };
    static const Int_t k_adc[2]    = {
      gUnpacker.get_data_id("SSD1","adc"),
      gUnpacker.get_data_id("SSD2","adc")
    };

    for(Int_t i=0; i<2; ++i){
      for(Int_t l=0; l<NumOfLayersSSD1; ++l){
	for(Int_t seg=0; seg<NumOfSegSSD1; ++seg){
	  Int_t nhit_a = gUnpacker.get_entries(k_device[i], l, seg, 0, k_adc[i]);
	  if( nhit_a==0 )
	    continue;
	  if( nhit_a != NumOfSamplesSSD ){
	    std::cerr << "#W SSD" << i+1
		      << " layer:" << l
		      << " seg:" << seg
		      << " the number of samples is wrong : "
		      << nhit_a << "/" << NumOfSamplesSSD << std::endl;
	    continue;
	  }

	  std::vector<Double_t> adc(NumOfSamplesSSD);
	  std::vector<Double_t> tdc(NumOfSamplesSSD);
	  Double_t peak_height = -1.;
	  Double_t peak_sample = -1.;
	  for(Int_t m=0; m<NumOfSamplesSSD; ++m){
	    adc[m] = gUnpacker.get(k_device[i], l, seg, 0, k_adc[i], m);
	    tdc[m] = (m+1)*25.;
	    if( peak_height < adc[m] ){
	      peak_height = adc[m];
	      peak_sample = m;
	    }
	  }

	  static TGraph graph;
	  graph.Set(0);
	  //( NumOfSamplesSSD, &(tdc[0]), &(adc[0]) );

	  m_adc[l+i*NumOfLayersSSD1][seg] = peak_height;
	  m_tdc[l+i*NumOfLayersSSD1][seg] = peak_sample;
	  Double_t pedestal = adc[0];
	  for(Int_t m=0; m<NumOfSamplesSSD; ++m){
	    adc[m] -= pedestal;
	    graph.SetPoint(m, tdc[m], adc[m]);
	  }

	  /*** SSD Waveform Fitting ***************************
	   *
	   *  f(x) = a * (x-b) * exp(-(x-b)/c)
	   *    a : scale factor
	   *    b : start timing
	   *    c : decay constant
	   *
	   *  f'(x)          = a/c * (-x+b+c) * exp(-(x-b)/c)
	   *  f'(x)|x=b+c    = 0.
	   *  f(b+c)         = a * c * exp(-1)
	   *  Sf(x)dx|b->inf = a * c^2
	   *
	   *  peak time = b + c
	   *  amplitude = a * c * exp(-1)
	   *  integral  = a * c^2
	   *
	   ****************************************************/

	  Double_t xmin =  40.;
	  Double_t xmax = 210.;
	  static TF1 func("func", "[0]*(x-[1])*exp(-(x-[1])/[2])", xmin, xmax);
	  func.SetParameter( 0, adc[3]*std::exp(1)/60. );
	  func.SetParLimits( 0, 0., 50000.*std::exp(-1.) );
	  func.SetParameter( 1, 40. );
	  func.SetParLimits( 1, 10., 100. );
	  func.FixParameter( 2, 50. );
	  graph.Fit("func", "RQ");
	  Double_t p[3];
	  func.GetParameters(p);

	  Double_t peak_time = p[1] + p[2];
	  // Double_t amplitude = p[0] * p[2] * std::exp(-1.);
	  Double_t de        = p[0]*p[2]*p[2];
	  Double_t chisqr    = func.GetChisquare() / func.GetNDF();
	  m_de[l+i*NumOfLayersSSD1][seg]   = de;
	  m_time[l+i*NumOfLayersSSD1][seg] = peak_time + m_ssdct[1];
	  m_chisqr[l+i*NumOfLayersSSD1][seg] = chisqr;
	}//for(seg)
      }//for(l)
#if 0
      gUnpacker.dump_data_device(k_device[i]);
#endif
    }
  }

  return true;
}
