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
#include "PsMaker.hh"
#include "MacroBuilder.hh"
#include "UserParamMan.hh"
#include "HodoParamMan.hh"
#include "MatrixParamMan.hh"

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
    HodoParamMan&   gHodo      = HodoParamMan::GetInstance();
    TText text;
    TText end;
  }

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{
  ConfMan& gConfMan = ConfMan::getInstance();
  gConfMan.initialize(argv);
  gConfMan.initializeHodoParamMan();

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
  
  gHttp.SetPort(8083);
  gHttp.Open();

  int tdcnbins=3000;
  double tdcmin=0;
  double tdcmax=2e6;
  gHttp.Register(gHist.createHodo(kBHD,"BHD",16,2,1024,0,1024,1024,0,2048));
  gHttp.Register(gHist.createHodo(kT0, "T0",5,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT0new,"T0new",1,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::QDC(kBHD,"_BHD_U",0,16,4,4,0,1000),"BHD");
  gHttp.Register(http::QDC(kBHD,"_BHD_D",1,16,4,4,0,1000),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_U",0,16,4,4),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_D",1,16,4,4),"BHD");

  gHttp.Register(http::QDC(kT0,"_T0_U",0,5,3,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,5,3,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");

  gHttp.Register(http::QDC(kT0new,"_T0new_U",0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::QDC(kT0new,"_T0new_D",1,1,1,1,0,1000),"T0new");
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
  // Hodoscope ------------------------------------------------------------
  {
    // data type
    const int nhodo=3;
    DetectorType hodo[nhodo]={kBHD,kT0,kT0new};
    TString hodo_name[nhodo]={"BHD","T0","T0new"};
    const int nsegs[nhodo]={16,5,1};
    const int nud[nhodo]={2,2,2};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
      const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
      for(int seg = 0; seg<nsegs[i]; ++seg){
	for(int ud=0; ud<nud[i]; ++ud){
	  for(int idata=0; idata<ndata; ++idata){
	    int hid;
	    int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	      hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	      //		std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	      hptr_array[hid]->Fill(tdc);	    
	      if(data[idata]==k_adc&&gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
		hid  = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
		hptr_array[hid]->Fill(tdc);	    
	      } 
	      if(data[idata]==k_leading){
		double time;
		gHodo.GetTime(1,0,seg,ud,tdc,time);
		std::cout<<tdc<<"  "<<time<<std::endl;
	      }
	    } // nhit
	  } // data
	}//ud
      }//seg
    } //hodo
  }
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  
  return 0;
}
  
}

