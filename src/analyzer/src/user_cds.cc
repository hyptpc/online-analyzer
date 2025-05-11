// -*- C++ -*-
// Author: Tomonori Takahashi

#include <iomanip>
#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TGFileBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <THttpServer.h>
#include <TKey.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TText.h>
#include <TTimeStamp.h>
#include <TFile.h>
#include <TTree.h>

#include "Controller.hh"
#include "HttpServer.hh"
#include "Updater.hh"

#include "user_analyzer.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "DAQNode.hh"
#include "filesystem_util.hh"

#include "ConfMan.hh"
#include "HistMaker.hh"
#include "DetectorID.hh"
#include "MacroBuilder.hh"

#define DEBUG    0
#define FLAG_DAQ 1

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    std::vector<TH1*> hptr_array;
    HistMaker&   gHist = HistMaker::getInstance();
    HttpServer&    gHttp = HttpServer::GetInstance();
    TText text;
    TText end;
  }

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{
  ConfMan::getInstance().initialize(argv);
  // gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  // gStyle->SetStatW(.42);
  // gStyle->SetStatH(.35);
  gStyle->SetStatW(.32);
  gStyle->SetStatH(.25);
  gStyle->SetPalette(55);
  
  // unpacker and all the parameter managers are initialized at this stage
  
  gHttp.SetPort(8088);
  gHttp.Open();

  gHttp.Register(gHist.createHodo(kCDH,"CDH",36,2,2048,-0.5,4095.5,1e4,0,1e7));
  gHttp.Register(gHist.createBLDC(kCDC,"CDC",118,16,true));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer (TCanvas)
  // for CDH
  gHttp.Register(http::QDC(kCDH,"CDH",0,36,6,6,0,1000));
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_U",0,36,6,6));
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_D",1,36,6,6));
  gHttp.Register(http::MHTDCMeanTime(kCDH,"_CDH",0,36,6,6));
  gHttp.Register(http::MHTDCHitPatMulti(kCDH,"CDH"));
  // for CDC
  for(int l=0; l<8; l++){
    gHttp.Register(http::CDCHitPat(kCDC,"CDC",l*16), "CDC");
  }
  for(int l=0; l<8; l++){
    gHttp.Register(http::CDCMulti(kCDC,"CDC",l*16), "CDC");
  }
  for(int l=0; l<8; l++){
    gHttp.Register(http::CDCTDC(kCDC,"CDC",l*16), "CDC");
  }

  for(int l=0;l<118;l++){ 
    gHttp.Register(http::BLDCWIRE(kCDC,"CDC",l,16,4,4), "CDC");
  }
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  
  return 0;
}
  
  //____________________________________________________________________________
  int
  process_end( void )
  {
  hptr_array.clear();
  return 0;
}

//____________________________________________________________________________
int
process_event( void )
{
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  static HistMaker&     gHist     = HistMaker::getInstance();
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  // Counter TDC ------------------------------------------------------------
  {
    static const int nmtdc=1;
    DetectorType mtdc[nmtdc]={kCDH};
    TString mtdc_name[nmtdc]={"CDH"};
    const int nsegs[nmtdc]={36};
    for(int i=0;i<nmtdc;i++)
      {
	// data type
	DetectorType kDET=mtdc[i];
	const int k_device = gUnpacker.get_device_id(mtdc_name[i].Data());
	const int k_l    = 1;
	const int k_t    = 2;      
	const int ndata=2;
	const int data[ndata]={k_l,k_t};
	const int type[ndata]={kTDC,kTDC2D};
	//	std::cout<<mtdc_name[i]<<"  "<<k_device<<std::endl;
	// TDC & HitPat & Multi
	int hid=-1;
	int multiplicity    = 0;
	for(int seg = 0; seg<nsegs[i]; ++seg){
	  for(int ud=0; ud<2; ++ud){	 
	    for(int idata=0; idata<ndata; ++idata){ 
	      int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	      //	      std::cout<<seg<<"  "<<ud<<"  "<<idata<<"  "<<nhit<<std::endl;
	      if(data[idata]==k_l&&nhit>0){
		hid  = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
		hptr_array[hid]->Fill(seg);	    
		// mul[ud]++;
		// ntdc[ud]=nhit;	  
	      }
	      for( int m=0; m<nhit; ++m ){
		int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
		hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
		hptr_array[hid]->Fill(tdc);	    
	      }
	    }
	  }
	}
      }
  }
  //return 0;
  // QDC&PADC ------------------------------------------------------------
  {
    // data type
    static const int nQDC=1;
    DetectorType adc[nQDC]={kCDH};
    TString adc_name[nQDC]={"CDH"};
    const int nsegs[nQDC]={36};
    for(int i=0;i<nQDC;i++){
      DetectorType kDET=adc[i];
      const int k_device = gUnpacker.get_device_id(adc_name[i].Data());
      for(int l = 0; l<nsegs[i]; ++l){
	int nhit = gUnpacker.get_entries(k_device, 0, l, 0, 0);
	if( nhit ){
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, 0, l, 0, 0 , m);
	    int hid  = gHist.getSequentialID(kDET, 0, kADC, l+1);
	    //	    std::cout<<k_device<<"  "<<kDET<<"  "<<i<<"  "<<l<<"  "<<uid<<"  "<<hid<<std::endl;
	    hptr_array[hid]->Fill(tdc);	    
	     //if(bhd_tdc[l]>0){
	     //  hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, l+1);
	     //  hptr_array[hid]->Fill(tdc);	    
	     //}
	  }
	}
      }    
    }
  }
  //return 0;

  // Chamber ------------------------------------------------------------
  const int nchm=1;
  DetectorType chm[nchm]={kCDC};
  TString chm_name[nchm]={"CDC"};
  const int nwires[nchm]={16};
  const int nlayers[nchm]={118};
  for(int i=0;i<nchm;i++)
    {
      // data type
      DetectorType kDET=chm[i];
      const int k_device = gUnpacker.get_device_id(chm_name[i].Data());
      const int k_l    = 0;
      const int k_t    = 1;
      
      // TDC & HitPat & Multi
      TString hname;
      int hid=-1;
      for(int l = 0; l<nlayers[i]; ++l){
	int multiplicity    = 0;
	for(int w = 0; w<nwires[i]; ++w){
	  int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	  if( nhit==0 ) continue;
	  // This wire fired at least one times.
	  ++multiplicity;
	  hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	  hptr_array[hid]->Fill(w, nhit);
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	    hid=gHist.getSequentialID(kDET,l+1,kTDC,w+1);
	    hptr_array[hid]->Fill(tdc);
	    hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	    hptr_array[hid]->Fill(tdc);
	  }	
	  nhit = gUnpacker.get_entries(k_device, l, 0, w, k_t);
	  if( nhit==0 ) continue;
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, l, 0, w, k_t , m);
	    hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	    hptr_array[hid]->Fill(tdc);
	  } 	
	}
	hid=gHist.getSequentialID(kDET,0,kMulti,l+1);
	hptr_array[hid]->Fill(multiplicity);
      }
    }
  //Scaler
  // const int k_device = gUnpacker.get_device_id("Scaler");	
  // int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
  // for(int module=0;module<1;module
  // for( int m=0; m<nhit; ++m ){
  //   if
  //   int tdc = gUnpacker.get(k_device, module, 0, channel);
  // }	    
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return 0;
}

}

