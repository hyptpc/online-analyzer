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
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

#include "HistMaker.hh"

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
HistMaker&   gHist   = HistMaker::getInstance();
HttpServer&    gHttp = HttpServer::GetInstance();
const auto& gUser    = UserParamMan::GetInstance();
const auto& gGeom    = DCGeomMan::GetInstance();
const auto& gSize    = DetSizeMan::GetInstance();

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
const int nsegs[nhodo]={5,15,5,34,9,1,8,1,1};
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
  ConfMan::getInstance().initializeHodoParamMan();
  ConfMan::getInstance().initializeDCGeomMan();
  ConfMan::getInstance().initializeDetSizeMan();
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
  int port=8083;
  
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
  gHttp.Register(gHist.createHodo(kSFV, "SFV", 1, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kCOBO, "COBO", 8, 1, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,1.45e6,1.7e6));
  
  gHttp.Register(gHist.createTriggerFlag());
  // Chambers
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));

  // BTOF
  gHttp.Register(gHist.createBTOF(true));
  gHttp.Register(gHist.createCorrelation());
  gHttp.Register(gHist.createEventDisplay());

  
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
  // gHttp.Register(http::HTOFThreshold());
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

  std::map<int, std::vector<Int_t>> hit_seg_map;


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
    //std::cout << trigger_flag << std::endl;
  }

  // if(trigger_flag[trigger::kSpillOnEnd] || trigger_flag[trigger::kSpillOffEnd])
  //   return 0;

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
    for(int seg = 0; seg<NumOfSegBHT; ++seg){
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
	hit_seg_map[k_device].push_back(seg);
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
		hit_seg_map[k_device].push_back(seg);
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

  { // BAC
    DetectorType kDET=kBAC;
    const int k_device = gUnpacker.get_device_id("BAC");
    const int k_adc = gUnpacker.get_data_id("BAC", "adc");
    const int k_tdc = gUnpacker.get_data_id("BAC", "leading");
    const UInt_t tdc_min = gUser.GetParameter("BAC_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("BAC_TDC", 1);
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
	if (has_hit && seg == 4){
	  hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	  hptr_array[hid]->Fill(seg);
	  hit_seg_map[k_device].push_back(0);
	  multiplicity++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
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
	hit_seg_map[k_device].push_back(seg);
	multiplicity++;
      }
      if (hit_seg_map[k_device].size() >= 2){
	hid = gHist.getSequentialID(kCorrelation, 0, kHitPat2D, 2);
	hptr_array[hid]->Fill(hit_seg_map[k_device][0], hit_seg_map[k_device][1]);
      }
      // int hid_num = gHist.getSequentialID(kDET, 2, kADCwTDC, seg+1);
      // int hid_den = gHist.getSequentialID(kDET, 2, kADC, seg+1);
      // hid = gHist.getSequentialID(kDET, 0, kThreshold, 0);
      // // hptr_array[hid]->Reset("ICES");
      // hptr_array[hid]->Divide(hptr_array[hid_num], hptr_array[hid_num], 1.0, 1.0);
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  Bool_t has_hit_T1 = false;
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
	    if (tdc_min < tdc && tdc < tdc_max){
	      has_hit    = true;
	      has_hit_T1 = true;
	    }
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
  

  { // KVC
    DetectorType kDET=kKVC;
    const int k_device = gUnpacker.get_device_id("KVC");
    const int k_adc = gUnpacker.get_data_id("KVC", "adc");
    const int k_tdc = gUnpacker.get_data_id("KVC", "leading");
    const UInt_t tdc_min = gUser.GetParameter("KVC_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("KVC_TDC", 1);
    const Int_t n_seg = 8;
    const Int_t n_ch = 5;
    Int_t multiplicity[2] = {0, 0}; // 0: just multi., 1: for eff check
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
	  hit_seg_map[k_device].push_back(seg);
	  multiplicity[0]++;
	  // if ( has_hit_T1 && 2 <= seg && seg <= 5 ) multiplicity[1]++;
	  if ( has_hit_T1 ) multiplicity[1]++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity[0]);
    if(has_hit_T1){
      hid = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(multiplicity[1]);
    }
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
      Bool_t has_hit_ud[3] = {false, false, false};
      for(int ch=0; ch<n_ch; ++ch){
        Bool_t has_hit       = false;
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
	  if (ch != 2) hid = gHist.getSequentialID(kDET, 0, kHitPat, ch);
	  hptr_array[hid]->Fill(seg);
	}
      }
      if (has_hit_ud[0] && has_hit_ud[1]) {
	hid = gHist.getSequentialID(kDET, 0, kHitPat, 2);
	hptr_array[hid]->Fill(seg);
	multiplicity++;
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

  Bool_t has_hit_SFV = false;
  { // SFV
    DetectorType kDET=kSFV;
    const int k_device = gUnpacker.get_device_id("SFV");
    const int k_adc = gUnpacker.get_data_id("SFV", "adc");
    const int k_tdc = gUnpacker.get_data_id("SFV", "leading");
    const UInt_t tdc_min = gUser.GetParameter("SFV_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("SFV_TDC", 1);
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
	    if (tdc_min < tdc && tdc < tdc_max) {
	      has_hit = true;
	      has_hit_SFV = true;
	    }
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

  { // SAC3
    DetectorType kDET=kSAC3;
    const int k_device = gUnpacker.get_device_id("SAC3");
    const int k_adc = gUnpacker.get_data_id("SAC3", "adc");
    const int k_tdc = gUnpacker.get_data_id("SAC3", "leading");
    const UInt_t tdc_min = gUser.GetParameter("SAC3_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("SAC3_TDC", 1);
    const Int_t n_seg = 1;
    const Int_t n_ch = 1;
    Int_t multiplicity[2] = {0, 0};
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
	  multiplicity[0]++;
	  if (has_hit_SFV) multiplicity[1]++;
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity[0]);

    if (has_hit_SFV) {
      hid = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(multiplicity[1]);
    } 
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
	    // std::cout << seg << " " << tdc << std::endl;
            hptr_array[hid]->Fill(tdc);
	    break;
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


      //------------------------------------------------------------------
    // BTOF
    //------------------------------------------------------------------
  
    {
      // Unpacker
      static const Int_t k_d_bht  = gUnpacker.get_device_id("BHT");
      static const Int_t k_d_bh2  = gUnpacker.get_device_id("BH2");
      static const Int_t k_bht_tdc    = gUnpacker.get_data_id("BHT", "leading");
      static const Int_t k_bh2_tdc    = gUnpacker.get_data_id("BH2", "leading");

      // HodoParam
      static const Int_t cid_bht  = 1;
      static const Int_t cid_bh2  = 6;
      static const Int_t plid     = 0;
      
      // Sequential ID
      static const Int_t btof_id  = gHist.getSequentialID(kMisc, 0, kTDC);
      static const auto& hodoMan = HodoParamMan::GetInstance();
      // TDC gate range
      //static const UInt_t tdc_min_bht = gUser.GetParameter("BHT_TDC", 0);
      //static const UInt_t tdc_max_bht = gUser.GetParameter("BHT_TDC", 1);
      static const UInt_t tdc_min_bht = 1400000;
      static const UInt_t tdc_max_bht = 1500000;
      
      
      static const UInt_t tdc_min_bh2 = gUser.GetParameter("BH2_TDC", 0);
      static const UInt_t tdc_max_bh2 = gUser.GetParameter("BH2_TDC", 1);


      //BAC
      
      DetectorType kDET=kBAC;
      const int k_device = gUnpacker.get_device_id("BAC");
      const int k_adc = gUnpacker.get_data_id("BAC","adc");
      const Int_t n_seg = 5;
      const Int_t n_ch = 1;
      Int_t adc_bac = -9999;
      { // ADC
	int ch = 0;
	Int_t n = gUnpacker.get_entries(k_device, 0, 4, ch, k_adc);
	for(Int_t m=0; m<n; ++m){
	  UInt_t adc = gUnpacker.get(k_device, 0, 4, ch, k_adc, m);
	  adc_bac = adc;
	  //hid = gHist.getSequentialID(kDET, ch, kADC, seg+1);
	  //hptr_array[hid]->Fill(adc);
	  /*
            if(has_hit){
	    hid = gHist.getSequentialID(kDET, ch, kADCwTDC, seg+1);
	    //hptr_array[hid]->Fill(adc);
            } 
	  */
	} // nhit
      }
      
      
      // BH2
      Double_t t0  = 1e10;
      Double_t ofs = 10;
      for(Int_t seg=0; seg<NumOfSegBH2; ++seg) {
	Int_t nhitu = gUnpacker.get_entries(k_d_bh2, 0, seg, kU, k_bh2_tdc);
	Int_t nhitd = gUnpacker.get_entries(k_d_bh2, 0, seg, kD, k_bh2_tdc);
	for(Int_t mu=0; mu<nhitu; ++mu) {
	  auto tdcu = gUnpacker.get(k_d_bh2, 0, seg, kU, k_bh2_tdc, mu);
	  if (tdcu < tdc_min_bh2 || tdc_max_bh2 < tdcu) continue;
	  for(Int_t md=0; md<nhitd; ++md) {
	    auto tdcd = gUnpacker.get(k_d_bh2, 0, seg, kD, k_bh2_tdc, md);
	    if (tdcd < tdc_min_bh2 || tdc_max_bh2 < tdcd) continue;
	    Double_t bh2ut=-9999., bh2dt=-9999.;
	    hodoMan.GetTime(cid_bh2, plid, seg, kU, tdcu, bh2ut);
	    hodoMan.GetTime(cid_bh2, plid, seg, kD, tdcd, bh2dt);
	    Double_t bh2mt = (bh2ut + bh2dt)/2.;
	    t0 = bh2mt;
	    if (TMath::Abs(t0) > TMath::Abs(bh2mt)) {
	      hodoMan.GetTime(cid_bh2, plid, seg, 2, 0, ofs);
	      t0 = bh2ut;
	    }
	  }
	}
      }

      // BHT
      for(Int_t seg=0; seg<NumOfSegBHT; ++seg) {
	Int_t nhitu = gUnpacker.get_entries(k_d_bht, 0, seg, kU, k_bht_tdc);
	Int_t nhitd = gUnpacker.get_entries(k_d_bht, 0, seg, kD, k_bht_tdc);
	for(Int_t mu=0; mu<nhitu; ++mu) {
	  auto tdcu = gUnpacker.get(k_d_bht, 0, seg, kU, k_bht_tdc, mu);
	  if (tdcu < tdc_min_bht || tdc_max_bht < tdcu) continue;
	  for(Int_t md=0; md<nhitd; ++md) {
	    auto tdcd = gUnpacker.get(k_d_bht, 0, seg, kD, k_bht_tdc, md);
	    if (tdcd < tdc_min_bht || tdc_max_bht < tdcd) continue;
	    Double_t bhttu = -9999., bhttd = -9999.;
	    hodoMan.GetTime(cid_bht, plid, seg, kU, tdcu, bhttu);
	    hodoMan.GetTime(cid_bht, plid, seg, kD, tdcd, bhttd);
	    Double_t mt = (bhttu+bhttd)/2.;
	    Double_t btof = mt-(t0+ofs);
	    hptr_array[btof_id]->Fill(btof);
	    hptr_array[btof_id+1]->Fill(btof,adc_bac);
	    
	  }// if (tdc)
	}// if (nhit)
      }// for(seg)
    }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  
  {
    const int bht_id = gUnpacker.get_device_id("BHT");
    const int bh2_id = gUnpacker.get_device_id("BH2");
    const int htof_id = gUnpacker.get_device_id("HTOF");
    const int kvc_id = gUnpacker.get_device_id("KVC");
    
    hid = gHist.getSequentialID(kCorrelation, 0, kHitPat2D, 1);
    if (hit_seg_map[bh2_id].size() > 0 && hit_seg_map[bht_id].size() > 0){
      hptr_array[hid]->Fill(hit_seg_map[bht_id][0], hit_seg_map[bh2_id][0]);
    }
    hid = gHist.getSequentialID(kCorrelation, 0, kHitPat2D, 3);
    if (hit_seg_map[bh2_id].size() > 0 && hit_seg_map[htof_id].size() > 0){
      hptr_array[hid]->Fill(hit_seg_map[htof_id][0], hit_seg_map[bh2_id][0]);
    }
    hid = gHist.getSequentialID(kCorrelation, 0, kHitPat2D, 4);
    if (hit_seg_map[kvc_id].size() > 0 && hit_seg_map[htof_id].size() > 0){
      hptr_array[hid]->Fill(hit_seg_map[htof_id][0], hit_seg_map[kvc_id][0]);
    }
  }

  //EventDispaly
  {

    int det_hist_pat_id = 0;
    int det_hist_count_id = 0;
    int hist_id = 0;
    for(int i=0;i<sizeof(EvtDis_Det_name)/sizeof(EvtDis_Det_name[0]);i++){
      const Int_t k_device = gUnpacker.get_device_id(EvtDis_Det_name[i].Data());
      det_hist_pat_id = gHist.getSequentialID(kEventDisplay, 0, kHitPoly, ++hist_id);
      hptr_array[det_hist_pat_id]->Reset();
      det_hist_count_id = gHist.getSequentialID(kEventDisplay, 0, kHitPoly, ++hist_id);
      
      const auto& hit_seg = hit_seg_map[k_device];
      if(!hit_seg.empty()){
	for(int i=0;i<hit_seg.size();i++){
	  int seg = hit_seg[i];
	  int binid = 0;
	  if(k_device == gUnpacker.get_device_id("HTOF")){
	    if(seg == 0) binid = 1;
	    else if(seg == 1 || seg == 2) binid = 2;
	    else if(seg == 3 || seg == 4) binid = 3;
	    else binid = seg-1;
	  }
	  else{
	    binid = seg + 1;
	  }
	  hptr_array[det_hist_pat_id]->SetBinContent(binid,100.);
	  double bin_cont = hptr_array[det_hist_count_id]->GetBinContent(binid);
	  hptr_array[det_hist_count_id]->SetBinContent(binid, bin_cont+1);
	}
      }
    }
  }
  

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

  //std::cout << "\a" << std::endl;

  return 0;
}

}

