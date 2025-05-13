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
#define FLAG_DAQ 1
#define WIRE 1 // for Canvas only

namespace analyzer
{
using namespace hddaq;
using namespace hddaq::unpacker;
namespace
{
std::vector<TH1*> hptr_array;
HistMaker&   gHist = HistMaker::getInstance();
HttpServer&    gHttp = HttpServer::GetInstance();
TText text;
TText end;
TString outputname="tmp.root";

// chambers
const int nchm=4;
DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b};
TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b"};
const int nwires[nchm]={32,32,32,32};
const int nlayers[nchm]={8,8,8,8};
   
// hodoscopes
const int nhodo=5;
DetectorType hodo[nhodo]={kT0,kBAC,kSAC,kKVC1,kBH2};
TString hodo_name[nhodo]={"T0","BAC","SAC","KVC1","BH2"};
const int nsegs[nhodo]={5,5,9,5,12};
const int nud[nhodo]={2,1,1,3,3};
const int k_adc=0;
const int k_leading=1;
const int k_trailing=2;
const int ndata=3;
const int data[ndata]={k_leading,k_trailing,k_adc};
const int type[ndata]={kTDC,kTDC2D,kADC};

//    int nbinsqdc=512;
int nbinsqdc=1024;
int nbinshrtdc=int(2e3);
double hrtdcmin=0.0e6;
double hrtdcmax=3.0e6;    
double hrtdcmin2=0.0e6;
double hrtdcmax2=1.0e6;    

TString flagnames[16]={"SpillStart","SpillEnd",
  "Beam","Pion","Kaon2","Kaon3",
  "KxCDH1","KxCDH2","KxCDH3","KxCDH1xG","KaonxG",
  "PixCDH","PixPbF2","ExPbF2",
  "CDH cosmic","clock(10s)"};

//For FADC Output
    
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
  int port=8081;
  if(argv.size()==4){
    outputname=argv.at(3);
    port=8089;
  }
  if(argv.size()==5){
    outputname=argv.at(3);
    port=std::strtoull(argv.at(4).c_str(), nullptr, 10);
  }
  
  gHttp.SetPort(port);
  gHttp.Open();


  // Hodoscopes
  //  gHttp.Register(gHist.createHodo(kBHD		,"BHT"  ,16, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createBHT(kBHT, "BHT"  ,63, 2, nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kT0 ,"T0"   , 5, 2, nbinsqdc/4,-0.5,511.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kBH2, "BH2",12 , 3, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kBAC, "BAC", 5, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHTOF, "HTOF", 34, 5,
				  nbinsqdc/2, -0.5, 2047.5,
				  nbinshrtdc, hrtdcmin, hrtdcmax));
  gHttp.Register(gHist.createHodo(kKVC1, "KVC1", 5, 3, 2048,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kKVC2, "KVC2", 5, 5,
				  2048, -0.5, 2047.5,
				  nbinshrtdc, hrtdcmin, hrtdcmax));
  gHttp.Register(gHist.createHodo(kSAC, "SAC", 9, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  
  gHttp.Register(gHist.createBVH());
  gHttp.Register(gHist.createT1());
  gHttp.Register(gHist.createT2());
  gHttp.Register(gHist.createTriggerFlag());
  // Chambers
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  // BHT
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_U",0,63,8,8),"BHT");
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_D",1,63,8,8),"BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_U",63,0,8,8), "BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_D",63,1,8,8), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_U",63,0,8,8), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_D",63,1,8,8), "BHT");
  gHttp.Register(http::MHTDCHitPatMulti(kBHT,"_BHT"),"BHT");

  // Hodoscopes
  gHttp.Register(http::QDC(kT0,"_T0_U",0,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");
  gHttp.Register(http::MHTDCHitPatMulti(kT0,"_T0"),"T0");

  gHttp.Register(http::BH2ADCU());
  gHttp.Register(http::BH2ADCD());
  gHttp.Register(http::BH2TDCU());
  gHttp.Register(http::BH2TDCD());
  gHttp.Register(http::BH2TDCMT());
  gHttp.Register(http::BAC());
  gHttp.Register(http::HTOFADC());
  gHttp.Register(http::HTOFTDCMT());
  gHttp.Register(http::HTOFTDC2());
  gHttp.Register(http::KVC1());
  gHttp.Register(http::KVC2());
  gHttp.Register(http::SAC());
  gHttp.Register(http::BVHTDCTOT());
  gHttp.Register(http::T1T2());
  gHttp.Register(http::HitPat());
  gHttp.Register(http::Multiplicity());

  // //gHttp.Register(http::MHTDCTDC(kBAC,"_BAC",0,1,1,1),"BAC");
  // //gHttp.Register(http::QDC(kSAC,"_SAC",0,0,9,3,2,0,2000),"BAC");
  // //gHttp.Register(http::MHTDCTDC(kSAC,"_SAC",0,1,1,1),"BAC");

  // // gHttp.Register(http::QDC(kDEF,"_DEF_U",0,0,4,2,2,0,2000),"DEF");
  // // gHttp.Register(http::QDC(kDEF,"_DEF_D",1,0,4,2,2,0,2000),"DEF");
  // // gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_U",0,4,2,2),"DEF");
  // // gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_D",1,4,2,2),"DEF");
  // // gHttp.Register(http::MHTDCHitPatMulti(kDEF,"_DEF"),"DEF");

  // // gHttp.Register(http::QDC(kBTC,"_BTC_L",0,0,2,2,2,0,1000),"BTC");
  // // gHttp.Register(http::QDC(kBTC,"_BTC_R",1,0,2,2,2,0,1000),"BTC");
  // // gHttp.Register(http::MHTDCTDC(kBTC,"_BTC_L",0,2,2,2),"BTC");
  // // gHttp.Register(http::MHTDCTDC(kBTC,"_BTC_R",1,2,2,2),"BTC");
  // // gHttp.Register(http::MHTDCHitPatMulti(kBTC,"_BTC"),"BTC");

  // // gHttp.Register(http::MHTDCTDC(kT98PMT,"_T98PMT",0,1,2,2),"T98");

  // // gHttp.Register(http::MHTDCTDC(kT98MPPC,"_T98MPPC",0,4,2,2),"T98");
  // // gHttp.Register(http::MHTDCHitPatMulti(kT98MPPC,"_T98MPPC"),"T98");
  // Chambers except for CDC
  for(int i=0;i<nchm;i++){
    const std::string tmpstr=Form("_%s",chm_name[i].Data());
    gHttp.Register(http::BLDCHitPat(chm[i],tmpstr,8),chm_name[i]); 
    gHttp.Register(http::BLDCMulti(chm[i],tmpstr,8),chm_name[i]);
    gHttp.Register(http::BLDCTDC(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTOT(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTDCvsTOT(chm[i],tmpstr,8), chm_name[i]);
    for(int l=0;l<8;l++){ 
      gHttp.Register(http::BLDCWIRE(chm[i],tmpstr,l,nwires[i],8,4), chm_name[i]);
    }
  }

  // // TriggerFlag
  // gHttp.Register(http::MHTDCTDC(kTriggerFlag,"_TriggerFlag",0,32,8,4),"TriggerFlag");
  // gHttp.Register(http::MHTDCHitPatMulti(kTriggerFlag,"_TriggerFlag"),"TriggerFlag");
  // {
  //   int hid1 = gHist.getSequentialID(kTriggerFlag, 0, kHitPat, 1);  
  //   hptr_array[hid1]->GetXaxis()->SetTitle("");
  //   for( Int_t i=0; i<32; ++i ){
  //     int hid2 = gHist.getSequentialID(kTriggerFlag, 0, kTDC, i+1);
  //     hptr_array[hid2]->SetTitle(Form("%s_%s",hptr_array[hid2]->GetTitle(),flagnames[i].Data()));
  //     hptr_array[hid1]->GetXaxis()->SetBinLabel(i+1,flagnames[i]);
  //   }
  // }		   

  //=== set directory ===//
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  //=== set directory ===//


  //  std::cout << "Start " <<std::endl;
  gHttp.Begin();  
  return 0;
}
  
//____________________________________________________________________________
int
process_end( void )
{
  TFile *root_file=new TFile(outputname,"recreate");
  for(int i=0;i<hptr_array.size();i++)
    if(hptr_array[i]) hptr_array[i]->Write();

  TString pdfname=outputname.ReplaceAll(".root",".pdf");
  bool INIT=false;
  TIter next(gROOT->GetListOfCanvases());
  TCanvas* c=0;
  while((c=(TCanvas*)next()) ){
    std::cout<<c->GetName()<<std::endl;
    c->Write();
    if(!INIT){
      c->Print(pdfname+"[");
      INIT=true;
    }
    c->Print(pdfname);
  }
  c=new TCanvas();
  c->Print(pdfname+"]");
  delete c;
  root_file->Close();    
    
  hptr_array.clear();
  return 0;
}

//____________________________________________________________________________
int
process_event( void )
{
  // gSystem->ProcessEvents();
  // return 0;
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  static HistMaker&     gHist     = HistMaker::getInstance();
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  //  std::cout << "Start-1 " <<std::endl;

  static Int_t run_number = -1;
  if(run_number != gUnpacker.get_run_number()){
    for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
      hptr_array[i]->Reset();
    }
    run_number = gUnpacker.get_run_number();
  }
  auto event_number = gUnpacker.get_event_number();
  //  std::cout << "Start-2 " <<std::endl;

  // for (auto& h : hptr_array){
  //   h->SetTitle(h->GetName() + TString(Form(" run%05d", run_number)));
  // }


  Int_t hid;
  bool COSMIC=false;
  bool CLOCK=false;
  
  std::bitset<NumOfSegTFlag> trigger_flag;
  { // TriggerFlag
    const Int_t k_device = gUnpacker.get_device_id("TriggerFlag");
    const Int_t k_tdc = gUnpacker.get_data_id("TriggerFlag", "leading");
    for(Int_t seg=0; seg<NumOfSegTFlag; ++seg){
      Bool_t has_hit = false;
      { // TDC
	Int_t n = gUnpacker.get_entries(k_device, 0, seg, 0, k_tdc);
	for(Int_t m=0; m<n; ++m){
	  UInt_t tdc = gUnpacker.get(k_device, 0, seg, 0, k_tdc, m);
	  hid = gHist.getSequentialID(kTriggerFlag, 0, kTDC, seg+1);
	  hptr_array[hid]->Fill(tdc);
	  has_hit = true;
	}
      }
      if(has_hit){
	hid = gHist.getSequentialID(kTriggerFlag, 0, kHitPat);
	hptr_array[hid]->Fill(seg);
	trigger_flag.set(seg);
      }
    }
  }

  if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
    return 0;

  //if(COSMIC&&!CLOCK) return 0;
  // BHT ------------------------------------------------------------
  {
    // data type    
    DetectorType kDET=kBHT;
    const int k_l=0;
    const int k_t=1;
    const int k_device = gUnpacker.get_device_id("BHT");
    int multiplicity=0;
    int mul[2]={0,0};
    int segud[2]={};
    for(int seg = 0; seg<63; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<2; ++ud){	  
	int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, k_l);
	if( nhit==0 ) continue;
	segud[ud]=seg;
	// This wire fired at least one times.
	++multiplicity;
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
          //          if(tdc>1.15e6&&tdc<1.2e6){
          hid=gHist.getSequentialID(kDET,0,kHitPat,ud+1);
          hptr_array[hid]->Fill(seg);
          //          }
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
  //  std::cout << "Start -3" <<std::endl;
    
  // Hodoscope ------------------------------------------------------------

  {
    // T0    
    int i=0;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    int multiplicity=0;
    int mul[2]={0,0};
    int segud[2]={-1,-1};
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<nud[i]; ++ud){	  
        for(int idata=0; idata<ndata; ++idata){
          int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
          if(data[idata]==k_leading && nhit>0){
            hid = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
            hptr_array[hid]->Fill(seg);
            segud[ud]=seg;
            mul[ud]++;
            ntdc[ud]=nhit;	  
          }
          for( int m=0; m<nhit; ++m ){
            int val = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
            hid = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
            hptr_array[hid]->Fill(val);	    
            if(data[idata]==k_adc &&
               gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
              hid = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(val);
            } 
          } // nhit
        } // data
      }//ud
      if(ntdc[0]>0 && ntdc[1]>0){
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
      hptr_array[hid]->Fill(segud[0], segud[1]);	    
    }
  } //hodo


  {
    // BAC    
    int i=1;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    int mul=0;
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc=0;	  
      for(int idata=0; idata<ndata; ++idata){
	int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	if( data[idata]==k_leading && nhit>0){
	  //hid = gHist.getSequentialID(kDET, 0, kHitPat, 1);
	  //hptr_array[hid]->Fill(seg);
	  ntdc=nhit;	  
	  mul=1;
	}
	for( int m=0; m<nhit; ++m ){
	  int val = gUnpacker.get(k_device, 0, seg, 0, data[idata] , m);
	  hid = gHist.getSequentialID(kDET, 0, type[idata], seg+1);
	  hptr_array[hid]->Fill(val);	    
	  if(data[idata]==k_adc &&
	     gUnpacker.get_entries(k_device, 0, seg, 0, k_leading)){
	    hid = gHist.getSequentialID(kDET, 0, kADCwTDC, seg+1);
	    hptr_array[hid]->Fill(val);
	  } 
	} // nhit
      } // data
      // if( ntdc>0 ){
      //   hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
      //   hptr_array[hid]->Fill(mul);	    
      // }
    }//seg	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(mul);	    	    
  } //hodo

  {
    // SAC    
    int i=2;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    int mul=0;
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc=0;	  
      for(int idata=0; idata<ndata; ++idata){
	int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	if( data[idata]==k_leading && nhit>0){
	  //hid = gHist.getSequentialID(kDET, 0, kHitPat, 1);
	  //hptr_array[hid]->Fill(seg);
	  ntdc=nhit;	  
	  mul=1;
	}
	for( int m=0; m<nhit; ++m ){
	  int val = gUnpacker.get(k_device, 0, seg, 0, data[idata] , m);
	  hid = gHist.getSequentialID(kDET, 0, type[idata], seg+1);
	  hptr_array[hid]->Fill(val);	    
	  if(data[idata]==k_adc &&
	     gUnpacker.get_entries(k_device, 0, seg, 0, k_leading)){
	    hid = gHist.getSequentialID(kDET, 0, kADCwTDC, seg+1);
	    hptr_array[hid]->Fill(val);
	  } 
	} // nhit
      } // data
      // if( ntdc>0 ){
      //   hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
      //   hptr_array[hid]->Fill(mul);	    
      // }
    }//seg	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(mul);	    	    
  } //hodo

  {
    // KVC1    
    int i=3;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    int multiplicity=0;
    int mul[3]={0,0,0};
    int ML=0;
    int segud[3]={-1,-1,-1};
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc[3]={0,0,0};
      for(int ud=0; ud<nud[i]; ++ud){	  
        for(int idata=0; idata<ndata; ++idata){
          int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
          if( data[idata]==k_leading && nhit>0){
            hid = gHist.getSequentialID(kDET, 0, kHitPat, ud);
            hptr_array[hid]->Fill(seg);
            segud[ud]=seg;
            ntdc[ud]=nhit;	  
            mul[ud]++;
            if(seg!=4){
              ML++;
            }
          }
          for( int m=0; m<nhit; ++m ){
            int val = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
            hid = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
            hptr_array[hid]->Fill(val);	    
            if(data[idata]==k_adc &&
               gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
              hid = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(val);
            } 
          } // nhit
        } // data
      }//ud
      // if( ntdc[0]>0 && ntdc[1]>0){
      //   hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
      //   hptr_array[hid]->Fill(seg);
      //   multiplicity++;
      // }
    }//seg
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(ML);	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
    hptr_array[hid]->Fill(mul[0]);	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
    hptr_array[hid]->Fill(mul[1]);	    
    if( multiplicity>0){
      hid  = gHist.getSequentialID(kDET, 0, kHitPat2D, 0); 
      hptr_array[hid]->Fill(segud[0], segud[1]);	    
    }
  } //hodo

  { // KVC2
    const int k_device = gUnpacker.get_device_id("KVC2");
    const int k_adc = gUnpacker.get_data_id("KVC2", "adc");
    const int k_tdc = gUnpacker.get_data_id("KVC2", "leading");
    const Int_t n_seg = 4;
    const Int_t n_ch = 5;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kKVC2, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kKVC2, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kKVC2, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit && ch == 4){
	  hid = gHist.getSequentialID(kKVC2, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kKVC2, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  {
    // BH2    
    int i=4;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    int mul=0;
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<nud[i]; ++ud){	  
        for(int idata=0; idata<ndata; ++idata){
          int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
          if( ud==2 && data[idata]==k_leading && nhit>0){
            //hid = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
            //hptr_array[hid]->Fill(seg);
            mul=1;
            ntdc[ud]=nhit;	  
          }
          for( int m=0; m<nhit; ++m ){
            int val = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
            hid = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
            hptr_array[hid]->Fill(val);	    
            if(data[idata]==k_adc &&
               gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
              hid = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(val);
            } 
          } // nhit
        } // data
      }//ud
    }//seg	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(mul);	    	    
  } //hodo

  { // HTOF
    const int k_device = gUnpacker.get_device_id("HTOF");
    const int k_adc = gUnpacker.get_data_id("HTOF", "adc");
    const int k_tdc = gUnpacker.get_data_id("HTOF", "leading");
    const Int_t n_seg = 34;
    const Int_t n_ch = 5;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kHTOF, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kHTOF, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kHTOF, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit && ch == 2){
	  hid = gHist.getSequentialID(kHTOF, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kHTOF, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  { ///// BVH
    static const auto k_device = gUnpacker.get_device_id("BVH");
    static const auto k_leading = gUnpacker.get_data_id("BVH", "leading");
    static const auto k_trailing = gUnpacker.get_data_id("BVH", "trailing");
    static const auto tdc_hid = gHist.getSequentialID(kBVH, 0, kTDC);
    static const auto tot_hid = gHist.getSequentialID(kBVH, 0, kADC);
    static const auto hit_hid = gHist.getSequentialID(kBVH, 0, kHitPat);
    static const auto mul_hid = gHist.getSequentialID(kBVH, 0, kMulti);
    Int_t multi = 0;
    for (Int_t seg=0; seg<NumOfSegBVH; ++seg) {
      Bool_t is_in_range = false;
      for (Int_t i=0, n=gUnpacker.get_entries(k_device, 0, seg, 0, k_leading); i<n; ++i) {
	auto tdc = gUnpacker.get(k_device, 0, seg, 0, k_leading, i);
	auto tdc_t = gUnpacker.get(k_device, 0, seg, 0, k_trailing, i);
	auto tot = tdc - tdc_t;
	hptr_array[tdc_hid + seg]->Fill(tdc);
	hptr_array[tot_hid + seg]->Fill(tot);
	is_in_range = true;
      }
      if (is_in_range) {
	++multi;
	hptr_array[hit_hid]->Fill(seg);
      }
    }
    hptr_array[mul_hid]->Fill(multi);
  }

  { ///// T1
    const int k_device = gUnpacker.get_device_id("T1");
    const int k_adc = gUnpacker.get_data_id("T1", "adc");
    const int k_tdc = gUnpacker.get_data_id("T1", "leading");
    const Int_t n_ch = 2;
    Int_t multi = 0;
    for (Int_t seg=0; seg<NumOfSegT1; ++seg) {
      std::bitset<n_ch> has_hit;
      for(int ch=0; ch<n_ch; ++ch){
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kT1, 0, kTDC, seg*n_ch+ch+1);
            hptr_array[hid]->Fill(tdc);
	    has_hit.set(ch);
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kT1, 0, kADC, seg*n_ch+ch+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit[ch]){
              hid = gHist.getSequentialID(kT1, 0, kADCwTDC, seg*n_ch+ch+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
      }
      if(has_hit.count() == n_ch){
	hid = gHist.getSequentialID(kT1, 0, kHitPat);
	hptr_array[hid]->Fill(seg);
	multi++;
      }
    }
    hid = gHist.getSequentialID(kT1, 0, kMulti);
    hptr_array[hid]->Fill(multi);
  }

  { ///// T2
    const int k_device = gUnpacker.get_device_id("T2");
    const int k_adc = gUnpacker.get_data_id("T2", "adc");
    const int k_tdc = gUnpacker.get_data_id("T2", "leading");
    const Int_t n_ch = 2;
    Int_t multi = 0;
    for (Int_t seg=0; seg<NumOfSegT2; ++seg) {
      std::bitset<n_ch> has_hit;
      for(int ch=0; ch<n_ch; ++ch){
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kT2, 0, kTDC, seg*n_ch+ch+1);
            hptr_array[hid]->Fill(tdc);
	    has_hit.set(ch);
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kT2, 0, kADC, seg*n_ch+ch+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit[ch]){
              hid = gHist.getSequentialID(kT2, 0, kADCwTDC, seg*n_ch+ch+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
      }
      if(has_hit.count() == n_ch){
	hid = gHist.getSequentialID(kT2, 0, kHitPat);
	hptr_array[hid]->Fill(seg);
	multi++;
      }
    }
    hid = gHist.getSequentialID(kT2, 0, kMulti);
    hptr_array[hid]->Fill(multi);
  }

  // Aerogel 
  {
    // double ac_sum=0;
    // DetectorType kDET=kBAC;
    // const int k_device = gUnpacker.get_device_id("BAC");
    // const int nseg = 5;
    // const int k_ac_leading=1;
    // bool ACHIT=false;
    // // tdc leading
    // int nhit=gUnpacker.get_entries(k_device, 0, 4, 0, k_ac_leading);
    // if(nhit>0){
    //   for( int m=0; m<nhit; ++m ){
    //     int tdc = gUnpacker.get(k_device, 0, 4, 0, k_ac_leading , m);
    //     hid  = gHist.getSequentialID(kDET, 0, kTDC, 1);
    //     //	if(tdc<6.7e5&&tdc>6.3e5) ACHIT=true;
    //     if(tdc<1060&&tdc>1020) ACHIT=true;
    //     hptr_array[hid]->Fill(tdc);	    
    //   }
    // }
    // // adc
    // for(int seg=0; seg<nseg; ++seg){
    //   if (gUnpacker.get_entries(k_device, 0, seg, 0, 0)) {
    //     int adc = gUnpacker.get(k_device, 0, seg, 0, 0);
    //     hid = gHist.getSequentialID(kDET, 0, kADC, seg+1);
    //     //  std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
    //     // ac_sum+=tdc;
    //     hptr_array[hid]->Fill(adc);	    
    //     if(ACHIT){	   
    //       hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, seg+1);
    //       hptr_array[hid]->Fill(adc);	    
    //     } 
    //   } // nhit
    // }//ch

    // double ac_sum=0;
    // DetectorType kDET=kSAC;
    // const int k_device = gUnpacker.get_device_id("SAC");
    // const int nseg = 9;
    // const int k_ac_leading=1;
    // bool ACHIT=false;
    // // tdc leading
    // int nhit=gUnpacker.get_entries(k_device, 0, 4, 0, k_ac_leading);
    // if(nhit>0){
    //   for( int m=0; m<nhit; ++m ){
    //     int tdc = gUnpacker.get(k_device, 0, 4, 0, k_ac_leading , m);
    //     hid  = gHist.getSequentialID(kDET, 0, kTDC, 1);
    //     //	if(tdc<6.7e5&&tdc>6.3e5) ACHIT=true;
    //     //if(tdc<1060&&tdc>1020) ACHIT=true;
    //     hptr_array[hid]->Fill(tdc);	    
    //   }
    // }
    // // adc
    // for(int seg=0; seg<nseg; ++seg){
    //   if (gUnpacker.get_entries(k_device, 0, seg, 0, 0)) {
    //     int adc = gUnpacker.get(k_device, 0, seg, 0, 0);
    //     hid = gHist.getSequentialID(kDET, 0, kADC, seg+1);
    //     //  std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
    //     // ac_sum+=tdc;
    //     hptr_array[hid]->Fill(adc);	    
    //     // if(ACHIT){	   
    //     //   hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, seg+1);
    //     //   hptr_array[hid]->Fill(adc);	    
    //     } 
    //   } // nhit
    // }//ch

    // hid = gHist.getSequentialID(kDET, 0, kADC, 5);
    // hptr_array[hid]->Fill(ac_sum);	    
    // if(ACHIT){
    //   hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, 5);
    //   hptr_array[hid]->Fill(ac_sum);	    
    // }
  }
  // Chamber ------------------------------------------------------------
  // BLDCs
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

	  // leading
	  int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	  if( nhit==0 ) continue;
	  // This wire fired at least one times.
	  ++multiplicity;
	  hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	  hptr_array[hid]->Fill(w, nhit);
	  std::vector< int > leading_array(nhit);
	  int leading_size = nhit;
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	    hid=gHist.getSequentialID(kDET,l+1,kTDC,w+1);
	    hptr_array[hid]->Fill(tdc);
	    hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	    hptr_array[hid]->Fill(tdc);
	    leading_array[m] = tdc;
	  }	

	  // trailing
	  nhit = gUnpacker.get_entries(k_device, l, 0, w, k_t);
	  if( nhit==0 ) continue;
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, l, 0, w, k_t , m);
	    hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	    hptr_array[hid]->Fill(tdc);
	    // for TOT
	    if(m<leading_size){
	      hid=gHist.getSequentialID(kDET,0,kTOT,l+1);
	      hptr_array[hid]->Fill(leading_array[m] - tdc);
	      hid=gHist.getSequentialID(kDET,0,kADC2D,l+1);
	      hptr_array[hid]->Fill(leading_array[m], leading_array[m] - tdc);
	    }
	  } 
	} // for(int w = 0; w<nwires[i]; ++w){
	hid=gHist.getSequentialID(kDET,0,kMulti,l+1);
	hptr_array[hid]->Fill(multiplicity);
      }
    }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  //update
  if(gUnpacker.get_counter()%100 == 0){
    auto prev_level = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;
    http::UpdateCounterEfficiency();
    //http::UpdateBcOutEfficiency();
    //http::UpdateSdcInOutEfficiency();
    // http::UpdateT0PeakFitting();
    //http::UpdateTOTPeakFitting();
    //http::UpdateAFTEfficiency();
    gErrorIgnoreLevel = prev_level;
  }

  gSystem->ProcessEvents();
  return 0;
}

}

