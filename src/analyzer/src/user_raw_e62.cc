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
  
  gHttp.SetPort(8081);
  gHttp.Open();

  gHttp.Register(gHist.createQDC(kQDC1,"QDC1",32,1,1024,0,1024));
  gHttp.Register(gHist.createQDC(kQDC2,"QDC2",32,1,1024,0,1024));
  gHttp.Register(gHist.createQDC(kQDC3,"QDC3",32,1,1024,0,1024));
  gHttp.Register(gHist.createQDC(kSDD1,"SDD1",32,1,1024,0,4096));
  gHttp.Register(gHist.createMHTDC(kHULTDC1,"HULTDC1",32,2048));
  gHttp.Register(gHist.createMHTDC(kHULTDC2,"HULTDC2",32,2048));
  gHttp.Register(gHist.createMHTDC(kHULTDC3,"HULTDC3",64,6553600));
  gHttp.Register(gHist.createBLDC(kBPC,"BPC",8,15,true));
  gHttp.Register(gHist.createBLDC(kSDC,"SDC",8,16,true));
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  gHttp.Register(gHist.createBLDC(kFDC,"FDC",6,64,true));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::QDC(kQDC1,"QDC1",0,32,8,4,0,1000));
  gHttp.Register(http::QDC(kQDC2,"QDC2",0,32,8,4,0,1000));
  gHttp.Register(http::QDC(kQDC3,"QDC3",0,32,8,4,0,1000));
  gHttp.Register(http::QDC(kSDD1,"SDD1",0,32,8,4,0,4000));
  gHttp.Register(http::MHTDCTDC(kHULTDC1,"HULTDC1",0,32,8,4));
  gHttp.Register(http::MHTDCTDC(kHULTDC2,"HULTDC2",0,32,8,4));
  gHttp.Register(http::MHTDCTDC(kHULTDC3,"HULTDC3",0,64,8,8));
  // gHttp.Register(http::BLDCHitPat(kBLC1a,"_BLC1a",8));
  // gHttp.Register(http::BLDCMulti(kBLC1a,"_BLC1a",8));
  // gHttp.Register(http::BLDCHitPat(kBLC1b,"_BLC1b",8));
  // gHttp.Register(http::BLDCMulti(kBLC1b,"_BLC1b",8));
  // gHttp.Register(http::BLDCHitPat(kBLC2a,"_BLC2a",8));
  // gHttp.Register(http::BLDCMulti(kBLC2a,"_BLC2a",8));
  // gHttp.Register(http::BLDCHitPat(kBLC2b,"_BLC2b",8));
  // gHttp.Register(http::BLDCMulti(kBLC2b,"_BLC2b",8));

  for(int l=0;l<8;l++){ 
    gHttp.Register(http::BLDCWIRE(kBPC,"BPC",l,15,5,3), "BPC");
    gHttp.Register(http::BLDCWIRE(kSDC,"SDC",l,16,4,4), "SDC");
    gHttp.Register(http::BLDCWIRE(kBLC1a,"BLC1a",l,32,8,4), "BLC1a");
    gHttp.Register(http::BLDCWIRE(kBLC1b,"BLC1b",l,32,8,4), "BLC1b");
    gHttp.Register(http::BLDCWIRE(kBLC2a,"BLC2a",l,32,8,4), "BLC2a");
    gHttp.Register(http::BLDCWIRE(kBLC2b,"BLC2b",l,32,8,4), "BLC2b");
  }
  for(int l=0;l<6;l++){ 
    gHttp.Register(http::BLDCWIRE(kFDC,"FDC",l,64,8,8), "FDC");
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
  int bhd_tdc[32];
  {
    const int nmtdc=3;
    DetectorType mtdc[nmtdc]={kHULTDC1,kHULTDC2,kHULTDC3};
    TString mtdc_name[nmtdc]={"HULTDC1","HULTDC1","HULTDC3"};
    const int nsegs[nmtdc]={32,32,64};
    for(int i=0;i<nmtdc;i++)
      {
	// data type
	DetectorType kDET=mtdc[i];
	const int k_device = gUnpacker.get_device_id(mtdc_name[i].Data());
	const int k_l    = 0;
	const int k_t    = 1;      
	// TDC & HitPat & Multi
	int hid=-1;
	int multiplicity    = 0;
	for(int l = 0; l<nsegs[i]; ++l){
	  if(i==0)	  bhd_tdc[l]=-1;
	  int nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_l);
	  if( nhit ){
	    ++multiplicity;
	    hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	    hptr_array[hid]->Fill(l, nhit);
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, l, 0, k_l, m);
	      hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	      hptr_array[hid]->Fill(tdc);
	      if(i==0&&m==0) bhd_tdc[l]=tdc;	      
	    }	
	  }
	  nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_t);
	  if( nhit ){
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, l, 0, k_t, m);
	      hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	      hptr_array[hid]->Fill(tdc);	    
	    } 	
	  }
	}      
	hid=gHist.getSequentialID(kDET,0,kMulti,1);
	hptr_array[hid]->Fill(multiplicity);
      }
  }

  // QDC&PADC ------------------------------------------------------------
  {
    // data type
    DetectorType adc[4]={kQDC1,kQDC2,kQDC3,kSDD1};
    TString adc_name[4]={"QDC1","QDC2","QDC3","SDD1"};
    const int nsegs[4]={32,32,32,32};
    for(int i=0;i<4;i++){
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
	    if(i==2&&bhd_tdc[l]>0){
	      hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, l+1);
	      hptr_array[hid]->Fill(tdc);	    
	    }
	  }
	}
      }    
    }
  }
  // Chamber ------------------------------------------------------------
  const int nchm=7;
  DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b, kBPC,kFDC,kSDC};
  TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b", "BPC","FDC","SDC"};
  const int nwires[nchm]={32,32,32,32,15,64,16};
  const int nlayers[nchm]={8,8,8,8,8,6,8};
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

