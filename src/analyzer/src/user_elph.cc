// -*- C++ -*-

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
#include <TRandom3.h>

#include "HodoAnalyzer.hh"
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

#define DEBUG    0

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    std::vector<TH1*> hptr_array;
    HistMaker&   gHist = HistMaker::getInstance();
    HttpServer&    gHttp = HttpServer::GetInstance();
    //    HodoParamMan&   gHodo      = HodoParamMan::GetInstance();
    HodoAnalyzer *hodoAna;
    TFile *root_file;
  }

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{
  root_file=new TFile("tmp.root","recreate");
  
  ConfMan& gConfMan = ConfMan::getInstance();
  gConfMan.initialize(argv);
  // gConfMan.initializeHodoParamMan();
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
  
  gHttp.SetPort(8080);
  gHttp.Open();

  gHttp.Register(gHist.createBHT(kBHT,"BHT"  ,8, 2, 3000,0,5e6));
    
  int tdcnbins=4096;
  double tdcmin=0;
  double tdcmax=4096;
  int adcnbins=1024;
  double adcmin=0;
  double adcmax=4096;
  gHttp.Register(gHist.createHodo(kCDH,"CNC"      ,2,2,adcnbins,adcmin,adcmax,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT0, "T0"        ,1,2,adcnbins,adcmin,adcmax,tdcnbins,tdcmin,tdcmax)); // trigup
  gHttp.Register(gHist.createHodo(kDEF,"Ref"       ,2,2,adcnbins,adcmin,adcmax,tdcnbins,tdcmin,tdcmax)); //reference
  gHttp.Register(gHist.createHodo(kFinger,"Trigger",1,1,adcnbins,adcmin,adcmax,tdcnbins,tdcmin,tdcmax));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_U",0,8,4,2),"BHT");
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_D",1,8,4,2),"BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_U",8,0,4,2), "BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_D",8,1,4,2), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_U",8,0,4,2), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_D",8,1,4,2), "BHT");
  gHttp.Register(http::MHTDCHitPatMulti(kBHT,"_BHT"),"BHT");

  gHttp.Register(http::QDC(kCDH,"_CNC_U",0,0,2,2,2,0,1000),"CNC");
  gHttp.Register(http::QDC(kCDH,"_CNC_D",1,0,2,2,2,0,1000),"CNC");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CNC_U",0,2,2,2),"CNC");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CNC_D",1,2,2,2),"CNC");
  gHttp.Register(http::MHTDCMeanTime(kCDH,"_CNC",0,2,2,2),"CNC");

  gHttp.Register(http::QDC(kDEF,"_Ref_U",0,0,2,2,2,0,1000),"Ref");
  gHttp.Register(http::QDC(kDEF,"_Ref_D",1,0,2,2,2,0,1000),"Ref");
  gHttp.Register(http::MHTDCTDC(kDEF,"_Ref_U",0,2,2,2),"Ref");
  gHttp.Register(http::MHTDCTDC(kDEF,"_Ref_D",1,2,2,2),"Ref");
  gHttp.Register(http::MHTDCMeanTime(kDEF,"_Ref",0,2,2,2),"Ref");

  gHttp.Register(http::QDC(kT0,"_T0_U",0,0,1,2,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,0,1,2,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,1,2,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,1,2,2),"T0");
  gHttp.Register(http::MHTDCMeanTime(kT0,"_T0",0,1,2,2),"T0");
  
  gHttp.Register(http::QDC(kFinger,"_Trigger",0,0,2,2,2,0,1000),"Finger");
  gHttp.Register(http::MHTDCTDC(kFinger,"_Trigger",0,2,2,2),"Finger");


  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }

  gHttp.Begin();  
  std::cout<<__func__<<std::endl;
  return 0;
}
  
  //____________________________________________________________________________
  int
  process_end( void )
  {
    root_file->cd();
    for(int i=0;i<hptr_array.size();i++)
      if(hptr_array[i]) hptr_array[i]->Write();
    root_file->Close();    
    hptr_array.clear();
    return 0;
  }

//____________________________________________________________________________
int
process_event( void )
{
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  static HistMaker&     gHist     = HistMaker::getInstance();
  // RawData  *rawData=new RawData;
  // rawData->DecodeHits();
  // hodoAna=new HodoAnalyzer;
  // hodoAna->DecodeRawHits( rawData );

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  // BHT ------------------------------------------------------------
  {
    // data type
    int hid;
    DetectorType kDET=kBHT;
    const int k_l=0;
    const int k_t=1;
    const int k_device = gUnpacker.get_device_id("BHT");
    int multiplicity=0;
    int mul[2]={0,0};
    int segud[2]={};
    for(int seg = 0; seg<8; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<2; ++ud){	  
	int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, k_l);
	if( nhit==0 ) continue;
	segud[ud]=seg;
	// This wire fired at least one times.
	++multiplicity;
	hid=gHist.getSequentialID(kDET,0,kHitPat,ud+1);
	hptr_array[hid]->Fill(seg, nhit);
	mul[ud]++;
	ntdc[ud]=nhit;	  
	std::vector< int > leading_array;
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, k_l, m);
	  hid=gHist.getSequentialID(kDET,ud,kTDC,seg+1);
	  //	  if(ud==0&&seg==0) std::cout<<seg<<"  "<<m<<"  "<<tdc<<"  "<<hid<<std::endl;
	  hptr_array[hid]->Fill(tdc);
	  leading_array.push_back(tdc);
	}	
	// traling
	nhit = gUnpacker.get_entries(k_device, 0, seg, ud, k_t);
	if( nhit==0 ) continue;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, k_t , m);
	  hid=gHist.getSequentialID(kDET,ud,kTDC2D,seg+1);
	  hptr_array[hid]->Fill(tdc);
	  // for TOT
	  if(m<leading_size){
	    hid=gHist.getSequentialID(kDET,ud,kTOT,seg+1);
	    hptr_array[hid]->Fill(leading_array.at(m) - tdc);
	    hid=gHist.getSequentialID(kDET,ud,kADC2D,seg+1);
	    hptr_array[hid]->Fill(leading_array.at(m), leading_array.at(m) - tdc);
	  }
	}
      }//ud
      if(ntdc[0]>0&&ntdc[1]>0){
	hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	hptr_array[hid]->Fill(seg);	    
	multiplicity++;
      }
    }//seg
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
    hptr_array[hid]->Fill(mul[0]);	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
    hptr_array[hid]->Fill(mul[1]);	    
    if(multiplicity>0){
      hid  = gHist.getSequentialID(kDET, 0, kHitPat2D, 0); 
      hptr_array[hid]->Fill(segud[0],segud[1]);	    
    }
  } //bht
  
  // Hodoscope ------------------------------------------------------------
  int hid;
  {
    // data type
    const int nhodo=4;
    int Cid[nhodo]={DetIdCDH,DetIdDEF,DetIdT0,DetIdFinger}; 
    DetectorType hodo[nhodo]={kCDH, kDEF, kT0,kFinger};
    TString hodo_name[nhodo]={"CNC","DEF","T0","Trigger"};
    const int nsegs[nhodo]={2,2,1,1};
    const int nud[nhodo]={2,2,2,1};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};
    // hodo rawdata    
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
      const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
      int multiplicity=0;
      int mul[2]={0,0};
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	for(int ud=0; ud<nud[i]; ++ud){	  
	  for(int idata=0; idata<ndata; ++idata){
	    int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	    if(data[idata]==k_leading&&nhit>0){
	      hid  = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
	      hptr_array[hid]->Fill(seg);	    
	      mul[ud]++;
	      ntdc[ud]=nhit;	  
	    }
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	      hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	      hptr_array[hid]->Fill(tdc);	    
	      if(data[idata]==k_adc&&gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
		hid  = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
		hptr_array[hid]->Fill(tdc);	  		
	      } 
	    } // nhit
	  } // data
	}//ud
	if(ntdc[0]>0&&ntdc[1]>0){
	  hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);	    
	  multiplicity++;
	}
      }//seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(multiplicity);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(mul[1]);	    
    } //hodo
  }
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  if(hodoAna) delete hodoAna;  
  //  if(rawData) delete rawData;
  return 0;
}
  
}

