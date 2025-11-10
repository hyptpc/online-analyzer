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
#include "TpcPadHelper.hh"
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
const auto& gUser     = UserParamMan::GetInstance();
auto& gTpcPad = TpcPadHelper::GetInstance();

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
// const int nhodo=5;
// DetectorType hodo[nhodo]={kT0,kBAC,kSAC,kKVC1,kBH2};
// TString hodo_name[nhodo]={"T0","BAC","SAC","KVC1","BH2"};
// const int nsegs[nhodo]={5,5,9,5,12};
// const int nud[nhodo]={2,1,1,3,3};
const int nhodo=9;
DetectorType hodo[nhodo]={kT0,kBH2,kBAC,kHTOF,kKVC,kT1,kCVC,kSAC3,kSFV};
TString hodo_name[nhodo]={"T0","BH2","BAC","HTOF","KVC","T1","CVC","SAC3","SFV"};
const int nsegs[nhodo]={5,15,5,34,9,1,8,1,5};
const int nud[nhodo]={2,3,1,3,3,1,3,1,1};
const int k_adc=0;
const int k_leading=1;
const int k_trailing=2;
const int ndata=3;
const int data[ndata]={k_leading,k_trailing,k_adc};
const int type[ndata]={kTDC,kTDC2D,kADC};

//    int nbinsqdc=512;
int nbinsqdc=4098;
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
  ConfMan::getInstance().initializeUserParamMan();
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
  gHttp.Register(gHist.createHodo(kBH2, "BH2", 15, 3, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kBAC, "BAC",  5, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHTOF(kHTOF, "HTOF", 34, 3,
				  nbinsqdc/2, -0.5, 2047.5,
				  nbinshrtdc, hrtdcmin, hrtdcmax));
  gHttp.Register(gHist.createHodo(kKVC, "KVC", 8, 5,
				  2048, -0.5, 2047.5,
				  nbinshrtdc, hrtdcmin, hrtdcmax));
  gHttp.Register(gHist.createHodo(kT1, "T1", 1, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kCVC, "CVC", 8, 3, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kSAC3, "SAC3", 1, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kSFV, "SFV", 5, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kCOBO, "COBO", 8, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,16.75e6,16.8e6));
  
  gHttp.Register(gHist.createTriggerFlag());
  // Chambers
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  
  // TPC
  gHttp.Register(gHist.createTPC(true));
  
  
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
  gHttp.Register(http::HTOFADCU());
  gHttp.Register(http::HTOFADCD());
  gHttp.Register(http::HTOFADCSUM());
  gHttp.Register(http::HTOFTDCU());
  gHttp.Register(http::HTOFTDCD());
  gHttp.Register(http::HTOFTDCHT());
  gHttp.Register(http::KVC());
  gHttp.Register(http::T1());
  gHttp.Register(http::CVCADC());
  gHttp.Register(http::CVCTDC());
  gHttp.Register(http::SAC3());
  gHttp.Register(http::SFV());
  gHttp.Register(http::COBO());
  gHttp.Register(http::HitPat());
  gHttp.Register(http::Multiplicity());

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
  // TFile *root_file=new TFile(outputname,"recreate");
  // for(int i=0;i<hptr_array.size();i++)
  //   if(hptr_array[i]) hptr_array[i]->Write();

  // TString pdfname=outputname.ReplaceAll(".root",".pdf");
  // bool INIT=false;
  // TIter next(gROOT->GetListOfCanvases());
  // TCanvas* c=0;
  // while((c=(TCanvas*)next()) ){
  //   std::cout<<c->GetName()<<std::endl;
  //   c->Write();
  //   if(!INIT){
  //     c->Print(pdfname+"[");
  //     INIT=true;
  //   }
  //   c->Print(pdfname);
  // }
  // c=new TCanvas();
  // c->Print(pdfname+"]");
  // delete c;
  // root_file->Close();    
    
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

  for (auto& h : hptr_array){
    h->SetTitle(h->GetName() + TString(Form(" run%05d", run_number)));
  }


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
    const UInt_t tdc_min = gUser.GetParameter("BHT_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("BHT_TDC", 1);
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
	ntdc[ud] = 0;
	std::vector< int > leading_array;
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, k_l, m);
	  hid=gHist.getSequentialID(kDET,ud,kTDC,seg+1);
	  //	  if(ud==0&&seg==0) std::cout<<seg<<"  "<<m<<"  "<<tdc<<"  "<<hid<<std::endl;
	  hptr_array[hid]->Fill(tdc);
	  leading_array.push_back(tdc);
	  if (tdc_min < tdc && tdc < tdc_max) {
	    hid = gHist.getSequentialID(kDET,0,kHitPat,ud+1);
	    hptr_array[hid]->Fill(seg);
	    ntdc[ud]++;
	  }
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
    const UInt_t tdc_min = gUser.GetParameter("T0_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("T0_TDC", 1);
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
    // BH2    
    int i=1;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    const UInt_t tdc_min = gUser.GetParameter("BH2_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("BH2_TDC", 1);
    int mul=0;
    for(int seg = 0; seg<nsegs[i]; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<nud[i]; ++ud){	  
        for(int idata=0; idata<ndata; ++idata){
          int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
          if( ud==2 && data[idata]==k_leading && nhit>0){
            //hid = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
            //hptr_array[hid]->Fill(seg);
            // mul=1;
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
            if (ud==2 && data[idata]==k_leading) {
	      if (tdc_min < val && val < tdc_max) {
		hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
		hptr_array[hid]->Fill(seg);
		mul++;
	      }
	    }
          } // nhit
        } // data
      }//ud
    }//seg	    
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(mul);	    	    
  } //hodo


  {
    // BAC    
    int i=2;
    DetectorType kDET=hodo[i];
    const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
    const UInt_t tdc_min = gUser.GetParameter("BAC_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("BAC_TDC", 1);
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

  { // HTOF
    DetectorType kDET=kHTOF;
    const int k_device = gUnpacker.get_device_id("HTOF");
    const int k_adc = gUnpacker.get_data_id("HTOF", "adc");
    const int k_tdc = gUnpacker.get_data_id("HTOF", "leading");
    const UInt_t tdc_min = gUser.GetParameter("HTOF_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("HTOF_TDC", 1);
    const Int_t n_seg = 34;
    const Int_t n_ch = 3;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      Bool_t has_hit_ud[3] = {false, false, false};
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit){	
	  has_hit_ud[ch] = true;
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, ch);
	  hptr_array[hid]->Fill(seg);
	}
      }
      if (has_hit_ud[0] && has_hit_ud[1]) {
	hid = gHist.getSequentialID(kDET, 0, kHitPat, 3);
	hptr_array[hid]->Fill(seg);
	multiplicity++;
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo
  
  { // KVC
    DetectorType kDET=kKVC;
    const int k_device = gUnpacker.get_device_id("KVC");
    const int k_adc = gUnpacker.get_data_id("KVC", "adc");
    const int k_tdc = gUnpacker.get_data_id("KVC", "leading");
    const UInt_t tdc_min = gUser.GetParameter("KVC_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("KVC_TDC", 1);
    const Int_t n_seg = 8;
    const Int_t n_ch = 5;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit && ch == 4){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  { // T1
    DetectorType kDET=kT1;
    const int k_device = gUnpacker.get_device_id("T1");
    const int k_adc = gUnpacker.get_data_id("T1", "adc");
    const int k_tdc = gUnpacker.get_data_id("T1", "leading");
    const UInt_t tdc_min = gUser.GetParameter("T1_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("T1_TDC", 1);
    const Int_t n_seg = 1;
    const Int_t n_ch = 1;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  { // CVC
    DetectorType kDET=kCVC;
    const int k_device = gUnpacker.get_device_id("CVC");
    const int k_adc = gUnpacker.get_data_id("CVC", "adc");
    const int k_tdc = gUnpacker.get_data_id("CVC", "leading");
    const UInt_t tdc_min = gUser.GetParameter("CVC_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("CVC_TDC", 1);
    const Int_t n_seg = 8;
    const Int_t n_ch = 3;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit && ch == 2){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  { // SAC3
    DetectorType kDET=kSAC3;
    const int k_device = gUnpacker.get_device_id("SAC3");
    const int k_adc = gUnpacker.get_data_id("SAC3", "adc");
    const int k_tdc = gUnpacker.get_data_id("SAC3", "leading");
    const UInt_t tdc_min = gUser.GetParameter("SAC3_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("SAC3_TDC", 1);
    const Int_t n_seg = 1;
    const Int_t n_ch = 1;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	{ // ADC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_adc);
          for(Int_t m=0; m<n; ++m){
            UInt_t adc = gUnpacker.get(k_device, 0, seg, ch, k_adc, m);
            hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
            hptr_array[hid]->Fill(adc);
            if(has_hit){
              hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(adc);
            } 
          } // nhit
	}
	if(has_hit){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  { // SFV
    DetectorType kDET=kSFV;
    const int k_device = gUnpacker.get_device_id("SFV");
    const int k_adc = gUnpacker.get_data_id("SFV", "adc");
    const int k_tdc = gUnpacker.get_data_id("SFV", "leading");
    const UInt_t tdc_min = gUser.GetParameter("SFV_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("SFV_TDC", 1);
    const Int_t n_seg = 5;
    const Int_t n_ch = 1;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	Bool_t has_hit = false;
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	    if (tdc_min < tdc && tdc < tdc_max)
	      has_hit = true;
	  }
	}
	if(has_hit){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo


  { // COBO
    DetectorType kDET=kCOBO;
    const int k_device = gUnpacker.get_device_id("COBO");
    const int k_adc = gUnpacker.get_data_id("COBO", "adc");
    const int k_tdc = gUnpacker.get_data_id("COBO", "leading");
    const Int_t n_seg = 8;
    const Int_t n_ch = 1;
    Int_t multiplicity = 0;
    for(int seg = 0; seg<n_seg; ++seg){
      for(int ch=0; ch<n_ch; ++ch){
	{ // TDC
	  Int_t n = gUnpacker.get_entries(k_device, 0, seg, ch, k_tdc);
          for(Int_t m=0; m<n; ++m){
            UInt_t tdc = gUnpacker.get(k_device, 0, seg, ch, k_tdc, m);
            hid = gHist.getSequentialID(kDET, ch, kTDC, seg+1);
            hptr_array[hid]->Fill(tdc);
	  }
	}
	hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	hptr_array[hid]->Fill(seg);
	multiplicity++;
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo


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


  // TPC -----------------------------------------------------------
  /*
  //UInt_t cobo_data_size = 0;
  {
    //if(cobo_data_size > 0){
      static const Int_t k_device   = gUnpacker.get_device_id( "TPC" );
      static const Int_t k_adc      = gUnpacker.get_data_id( "TPC", "adc" );
      // static const Int_t k_tdc_high = gUnpacker.get_data_id( "TPC", "tdc_high" );
      // static const Int_t k_tdc_low  = gUnpacker.get_data_id( "TPC", "tdc_low" );
      // sequential id
      static const Int_t tpca_id   = gHist.getSequentialID( kTPC, 0, kADC );
      static const Int_t tpct_id   = gHist.getSequentialID( kTPC, 0, kTDC );
      static const Int_t rms_id    = gHist.getSequentialID( kTPC, 0, kPede );
      static const Int_t tpca2d_id = gHist.getSequentialID( kTPC, 0, kADC2D );
      static const Int_t tpcmul_id = gHist.getSequentialID( kTPC, 0, kMulti );
      static const Int_t agetmul_id = gHist.getSequentialID( kTPC, 3, kMulti );
      static const Int_t amulmax_id = gHist.getSequentialID( kTPC, 4, kMulti );
      
      hptr_array[agetmul_id]->Reset();
      
      // FADC
      Int_t n_active_pad = 0;
      for( Int_t layer=0; layer<32; ++layer ){
	for( Int_t ch=0; ch<272; ++ch ){
	  Int_t pad = gTpcPad.GetPadId( layer, ch );
	  if( pad < 0 ) continue;
	  Int_t nhit = gUnpacker.get_entries( k_device, layer, 0, ch, k_adc );
	  if( nhit == 0 ){
	    hptr_array[tpca2d_id]->SetBinContent( pad+1, 0 );
	    hptr_array[tpca2d_id+1]->SetBinContent( pad+1, 0 );
	    hptr_array[tpca2d_id+2]->SetBinContent( pad+1, 0 );
	    continue;
	  }
	  // if( nhit != 200 ){
	  //   hddaq::cerr << "#W NumOfTimeBucket is wrong " << nhit << "/200 "
	  //               << "(layer=" << layer << ", row=" << ch << ", pad="
	  //               << pad << ")" << std::endl;
	  // }
	  std::vector<Double_t> fadc( nhit );
	  for( Int_t i=0; i<nhit; ++i ){
	    Int_t adc = gUnpacker.get( k_device, layer, 0, ch, k_adc, i );
	    fadc[i] = adc;
	  }
	  Double_t mean = TMath::Mean( nhit, fadc.data() );
	  Double_t rms = TMath::RMS( nhit, fadc.data() );
	  Double_t max_adc = TMath::MaxElement( nhit, fadc.data() );
	  Int_t loc_max = TMath::LocMax( nhit, fadc.data() );

	  if( max_adc - mean <= 0 ) continue;
	  Int_t aget = gTpcPad.GetParam( layer, ch )->AGetId();
	  Int_t asad = gTpcPad.GetParam( layer, ch )->AsAdId();
	  hptr_array[agetmul_id]->Fill( asad*4+aget );

	  hptr_array[tpca_id]->Fill( max_adc - mean );
	  hptr_array[tpct_id]->Fill( loc_max );
	  hptr_array[rms_id]->Fill( rms );
	  hptr_array[tpca2d_id]->SetBinContent( pad+1, max_adc - mean );
	  hptr_array[tpca2d_id+1]->SetBinContent( pad+1, rms );
	  hptr_array[tpca2d_id+2]->SetBinContent( pad+1, loc_max );
	  Double_t pad_z = gTpcPad.GetPoint( pad ).Z();
	  Double_t pad_x = gTpcPad.GetPoint( pad ).X();
	  if( max_adc - mean > 0 ){
	    ++n_active_pad;
	    hptr_array[tpca2d_id+3]->Fill( pad_z, pad_x );
	  }
	}
      }
      // std::cout << "active pad = " << n_active_pad << std::endl;
      hptr_array[tpcmul_id]->Fill( n_active_pad );
      hptr_array[amulmax_id]->Fill( hptr_array[agetmul_id]->GetMaximum() );

      // TDC (Time Stamp)
      // UInt_t tdc_h = gUnpacker.get( k_device, 0, 0, 0, k_tdc_high );
      // UInt_t tdc_l = gUnpacker.get( k_device, 0, 0, 0, k_tdc_low );
      // std::cout << tdc_h << ", " << tdc_l << std::endl;

#if 0
      // Debug, dump data relating this detector
      gUnpacker.dump_data_device(k_device);
#endif
    }
    { ///// TPC-CLOCK
      static const auto device_id = gUnpacker.get_device_id("HTOF");
      static const auto leading_id = gUnpacker.get_data_id("HTOF", "fpga_leading");
      static const auto tdc_hid = gHist.getSequentialID(kTPC, 1, kTDC);
      static const auto mul_hid = gHist.getSequentialID(kTPC, 1, kMulti);
      static const Int_t seg = 34;
      UInt_t multiplicity = 0;
      for(Int_t m=0, n=gUnpacker.get_entries(device_id, 0, seg, 0, leading_id);
          m<n; ++m){
        auto tdc = gUnpacker.get(device_id, 0, seg, 0, leading_id, m);
        if(tdc != 0){
          hptr_array[tdc_hid]->Fill(tdc);
          ++multiplicity;
	}
      }
      hptr_array[mul_hid]->Fill(multiplicity);
      //}
  } // TPC
  */

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

