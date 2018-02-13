// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <vector>

#include <TGFileBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>

#include <UnpackerManager.hh>

#include "Controller.hh"
#include "user_analyzer.hh"

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCAnalyzerOld.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "HistMaker.hh"
#include "MacroBuilder.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    const HistMaker&       gHist     = HistMaker::getInstance();
    const UnpackerManager& gUnpacker = GUnpacker::get_instance();
    const DCGeomMan&       gGeom     = DCGeomMan::GetInstance();

    const Double_t dist_FF = 1200.;

    std::vector<TH1*> hptr_array;

    enum HistName
      {
	// Inside Target
	i_Zn495, i_Zn338, i_Zn167_5, i_Z0, i_Z150, i_Z240,
	// Detector positions
	i_SFT, i_EG, i_DC1, i_SAC,
	NHist
      };
    static const Double_t TgtOrg_plus[] =
      {
	// Inside Target
	-495., -338., -167.5, 0, 150., 240.,
	// Detector
	400., 540., 733., 866.
      };

    static const std::string posName[] =
      {
	// Inside Target
	"Tgt Z=-495", "Tgt Z=-338", "Tgt Z=-167.5", "Tgt Z=0",
	"Tgt Z=150",  "Tgt Z=240",
	// Detectors
	"SFT", "End guard", "DC1", "KbAC"
      };

    static const Double_t posSSD[] =
      {
	-100, 100, 100, 100
      };

    static Double_t FF_plus = 0;
  }

//____________________________________________________________________________
Int_t
process_begin( const std::vector<std::string>& argv )
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  gConfMan.InitializeDCGeomMan();
  gConfMan.InitializeDCTdcCalibMan();
  gConfMan.InitializeDCDriftParamMan();
  gConfMan.InitializeUserParamMan();
  if( !gConfMan.IsGood() ) return -1;
  // unpacker and all the parameter managers are initialized at this stage

  FF_plus = UserParamMan::GetInstance().GetParameter("FF_PLUS", 0);
  std::cout << "#D : FF+" << FF_plus << std::endl;

  // Make tabs
  hddaq::gui::Controller& gCon = hddaq::gui::Controller::getInstance();
  TGFileBrowser *tab_hist  = gCon.makeFileBrowser("Hist");
  TGFileBrowser *tab_macro = gCon.makeFileBrowser("Macro");

  // Add macros to the Macro tab
  tab_macro->Add(macro::Get("clear_all_canvas"));
  tab_macro->Add(macro::Get("clear_canvas"));
  tab_macro->Add(macro::Get("split22"));
  tab_macro->Add(macro::Get("split32"));
  tab_macro->Add(macro::Get("split33"));
  tab_macro->Add(macro::Get("dispBeamProfile_e40"));
  tab_macro->Add(macro::Get("dispBcOutFF"));
  tab_macro->Add(macro::Get("dispSSD1Profile"));

  // Add histograms to the Hist tab
  HistMaker& gHist = HistMaker::getInstance();
  //BcOut
  {
    TList *sub_dir = new TList;
    const char* nameSubDir = "BcOut";
    sub_dir->SetName(nameSubDir);
    Int_t unique_id = gHist.getUniqueID(kMisc, 0, kHitPat);
    // Profile X
    for(Int_t i = 0; i<NHist; ++i){
      char* title = i==i_Z0?
	Form("%s_X %s (FF+%d)", nameSubDir, posName[i].c_str(), (Int_t)FF_plus) :
	Form("%s_X %s", nameSubDir, posName[i].c_str());
      sub_dir->Add(gHist.createTH1(unique_id++, title,
				    400,-200,200,
				    "x position [mm]", ""));
    }
    // Profile Y
    for(Int_t i = 0; i<NHist; ++i){
      char* title = i==i_Z0?
	Form("%s_Y %s (FF+%d)", nameSubDir, posName[i].c_str(), (Int_t)FF_plus) :
	Form("%s_Y %s", nameSubDir, posName[i].c_str());
      sub_dir->Add(gHist.createTH1(unique_id++, title,
				    200,-100,100,
				    "y position [mm]", ""));
    }
    tab_hist->Add(sub_dir);
    // Profile XY
    for(Int_t i = 0; i<NHist; ++i){
      char* title = i==i_Z0?
	Form("%s_XY %s (FF+%d)", nameSubDir, posName[i].c_str(), (Int_t)FF_plus) :
	Form("%s_YY %s", nameSubDir, posName[i].c_str());
      sub_dir->Add(gHist.createTH2(unique_id++, title,
				   400,-200,200, 200,-100,100,
				   "x position [mm]", "y position [mm]"));
    }
    tab_hist->Add(sub_dir);
  }

  //SSD1
  {
    TList *sub_dir = new TList;
    const char* nameSubDir = "SSD1";
    sub_dir->SetName(nameSubDir);
    Int_t unique_id = gHist.getUniqueID(kMisc, 0, kHitPat2D);
    const char* nameLayer[] = { "Y0", "X0", "Y1", "X1" };
    for(Int_t l=0; l<NumOfLayersSSD1; ++l){
      char* title = Form("%s_HitPos_%s FF+%d", nameSubDir, nameLayer[l], (Int_t)posSSD[l]);
      sub_dir->Add(gHist.createTH1(unique_id++, title,
				   200,-100,100,
				   "Position [mm]", ""));
    }
    tab_hist->Add(sub_dir);
  }

  // Set histogram pointers to the vector sequentially.
  // This vector contains both TH1 and TH2.
  // Then you need to do down cast when you use TH2.
  if(0 != gHist.setHistPtr(hptr_array)){return -1;}

  gStyle->SetOptStat(1110);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetStatW(.42);
  gStyle->SetStatH(.3);

  return 0;
}

//____________________________________________________________________________
Int_t
process_end( void )
{
  return 0;
}

//____________________________________________________________________________
Int_t
process_event( void )
{
  static const Int_t xpos_id  = gHist.getSequentialID(kMisc, 0, kHitPat);
  static const Int_t ypos_id  = gHist.getSequentialID(kMisc, 0, kHitPat, NHist+1);
  static const Int_t xypos_id = gHist.getSequentialID(kMisc, 0, kHitPat, NHist*2+1);

  debug::ObjectCounter::Check();

  EventAnalyzer event;
  event.DecodeRawData();
  event.DecodeDCRawHits();

  Int_t ntBcOut = event.TrackSearchBcOut();
  for( Int_t i=0; i<ntBcOut; ++i ){
    const DCLocalTrack* track = event.GetTrackBcOut(i);
    if( !track ) continue;
    Double_t x0 = track->GetX0(); Double_t y0 = track->GetY0();
    // Double_t u0 = track->GetU0(); Double_t v0 = track->GetV0();
    for( Int_t j=0; j<NHist; ++j ){
      Double_t z = dist_FF+FF_plus+TgtOrg_plus[j];
      Double_t x = track->GetX(z); Double_t y = track->GetY(z);
      hptr_array[xpos_id+j]->Fill(x);
      hptr_array[ypos_id+j]->Fill(y);
      hptr_array[xypos_id+j]->Fill(x, y);
    }
  }

  return 0;

  /////////// BcOut
  {
    DCRHC BcOutAna(DetIdBcOut);
    bool BcOutTrack = BcOutAna.TrackSearch(9);

    if(BcOutTrack){
      for(Int_t i = 0; i<NHist; ++i){
	Double_t xpos = BcOutAna.GetPosX(dist_FF+FF_plus+TgtOrg_plus[i]);
	Double_t ypos = BcOutAna.GetPosY(dist_FF+FF_plus+TgtOrg_plus[i]);
	hptr_array[xpos_id+i]->Fill(xpos);
	hptr_array[ypos_id+i]->Fill(ypos);
	hptr_array[xypos_id+i]->Fill(xpos, ypos);
      }
    }
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  ////////// SSD1
  {
    // data type
    static const Int_t k_device = gUnpacker.get_device_id("SSD1");
    static const Int_t k_adc    = gUnpacker.get_data_id("SSD1","adc");
    static const Int_t k_flag   = gUnpacker.get_data_id("SSD1","flag");
    // sequential id
    static const Int_t hit_id = gHist.getSequentialID(kMisc, 0, kHitPat2D);

    bool chit_flag[NumOfLayersSSD1][NumOfSegSSD1];

    for(Int_t l=0; l<NumOfLayersSSD1; ++l){
      for(Int_t seg=0; seg<NumOfSegSSD1; ++seg){
	chit_flag[l][seg] = false;
	// ADC
	Int_t nhit_a = gUnpacker.get_entries(k_device, l, seg, 0, k_adc);
	if( nhit_a==0 ) continue;
	if( nhit_a != NumOfSamplesSSD ){
	  std::cerr << "#W SSD1 layer:" << l << " seg:" << seg
		    << " the number of samples is wrong : "
		    << nhit_a << "/" << NumOfSamplesSSD << std::endl;
	  continue;
	}
	Int_t  adc[nhit_a];
	bool slope[nhit_a-1];
	Int_t  peak_height   = -1;
	Int_t  peak_position = -1;
	for(Int_t m=0; m<nhit_a; ++m){
	  adc[m] = gUnpacker.get(k_device, l, seg, 0, k_adc, m);
	  adc[m] -= adc[0];
	  if( m>0 )
	    slope[m] = adc[m]>adc[m-1];
	  if( adc[m]>peak_height ){
	    peak_height   = adc[m];
	    peak_position = m;
	  }
	}
	// Zero Suppression Flag
	Int_t  nhit_flag = gUnpacker.get_entries(k_device, l, seg, 0, k_flag);
	bool  hit_flag = false;
	if(nhit_flag != 0){
	  Int_t flag = gUnpacker.get(k_device, l, seg, 0, k_flag);
	  if(flag==1) hit_flag = true;
	}
	chit_flag[l][seg] = hit_flag &&
	  slope[0] &&
	  slope[1] && slope[2] && !slope[4] && !slope[5] && !slope[6];
	  //&& ( peak_height > 350 );
	if( peak_height>=0 && peak_position>=1 && chit_flag[l][seg] ){
	  Double_t wpos = gGeom.CalcWirePosition( l+7, seg+1 );
	  hptr_array[hit_id +l]->Fill( wpos );
	}
      }//for(seg)
    }//for(l)
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  return 0;
}

}
