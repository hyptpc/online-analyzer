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
#include <TGraph.h>

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
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "DCAnalyzer.hh"
#include "DCLocalTrack.hh"
#include "EventAnalyzer.hh"
#include "RawData.hh"

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
const auto& gTpcPad  = TpcPadHelper::GetInstance();
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
  ConfMan::getInstance().initializeTpcPadHelper();
  ConfMan::getInstance().initializeDCGeomMan();
  ConfMan::getInstance().initializeDCTdcCalibMan();
  ConfMan::getInstance().initializeDCDriftParamMan();
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
  int port=9090;
  
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

  gHttp.Register(gHist.createDAQ(true));
  // TPC
  gHttp.Register(gHist.createTPC(true));

  // BTOF
  gHttp.Register(gHist.createBTOF(true));
  gHttp.Register(gHist.createCorrelation());
  gHttp.Register(gHist.createEventDisplay());

  gHttp.Register(gHist.createBcOutTracking());

  
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  /*
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
  */

  gHttp.Register(http::TPC2D());
  gHttp.Register(http::TPCADCPAD());
  gHttp.Register(http::TPCAGETCond());
  gHttp.Register(http::TPCFADCAGET());
  gHttp.Register(http::SHS2D());


  //=== set directory ===//
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  //=== set directory ===//


  gHttp.Begin();  
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
  // gSystem->ProcessEvents();
  // return 0;
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  static HistMaker&     gHist     = HistMaker::getInstance();
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  static Int_t run_number = -1;
  if(run_number != gUnpacker.get_run_number()){
    for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
      hptr_array[i]->Reset();
    }
    run_number = gUnpacker.get_run_number();
  }
  auto event_number = gUnpacker.get_event_number();

  hptr_array[gHist.getSequentialID(kEventDisplay, 0, kHitPoly, 1)]->SetName(Form("EventDisplay evt%d",event_number));
  hptr_array[gHist.getSequentialID(kTPC, 0, kADC2D)]->SetName(Form("TPC_ADC2D evt%d",event_number));
  
  for (auto& h : hptr_array){
    h->SetTitle(h->GetName() + TString(Form(" run%05d", run_number)));
  }

  UInt_t cobo_data_size = 0;
#if FLAG_DAQ
  { ///// DAQ
    //___ node id
    static const Int_t k_eb = gUnpacker.get_fe_id("k18breb");
    std::vector<Int_t> vme_fe_id;
    std::vector<Int_t> hul_fe_id;
    std::vector<Int_t> ea0c_fe_id;
    std::vector<Int_t> cobo_fe_id;
    for(auto&& c : gUnpacker.get_root()->get_child_list()){
      if(!c.second) continue;
      TString n = c.second->get_name();
      auto id = c.second->get_id();
      if(n.Contains("vme"))
	vme_fe_id.push_back(id);
      if(n.Contains("hul"))
	hul_fe_id.push_back(id);
      if(n.Contains("easiroc"))
	ea0c_fe_id.push_back(id);
      if(n.Contains("cobo"))
	cobo_fe_id.push_back(id);
    }

    //___ sequential id
    static const Int_t eb_hid = gHist.getSequentialID(kDAQ, kEB, kHitPat);
    static const Int_t vme_hid = gHist.getSequentialID(kDAQ, kVME, kHitPat2D);
    static const Int_t hul_hid = gHist.getSequentialID(kDAQ, kHUL, kHitPat2D);
    static const Int_t ea0c_hid = gHist.getSequentialID(kDAQ, kEASIROC, kHitPat2D);
    static const Int_t cobo_hid = gHist.getSequentialID(kDAQ, kCoBo, kHitPat2D);

    { //___ EB
      auto data_size = gUnpacker.get_node_header(k_eb, DAQNode::k_data_size);
      hptr_array[eb_hid]->Fill(data_size);
    }

    { //___ VME
      for(Int_t i=0, n=vme_fe_id.size(); i<n; ++i){
	auto data_size = gUnpacker.get_node_header(vme_fe_id[i], DAQNode::k_data_size);
        hptr_array[vme_hid]->Fill(i, data_size);
      }
    }

    { // EASIROC & VMEEASIROC node
      for(Int_t i=0, n=ea0c_fe_id.size(); i<n; ++i){
        auto data_size = gUnpacker.get_node_header(ea0c_fe_id[i], DAQNode::k_data_size);
        hptr_array[ea0c_hid]->Fill(i, data_size);
      }
    }

    { //___ HUL node
      for(Int_t i=0, n=hul_fe_id.size(); i<n; ++i){
        auto data_size = gUnpacker.get_node_header(hul_fe_id[i], DAQNode::k_data_size);
        hptr_array[hul_hid]->Fill(i, data_size);
      }
    }

    { //___ CoBo node
      for(Int_t i=0, n=cobo_fe_id.size(); i<n; ++i){
        auto data_size = gUnpacker.get_node_header(cobo_fe_id[i], DAQNode::k_data_size);
        hptr_array[cobo_hid]->Fill(i, data_size);
	cobo_data_size += data_size;
      }
    }
  }

#endif

  std::map<int, std::vector<Int_t>> hit_seg_map; 


  Int_t hid;
  bool COSMIC=false;
  bool CLOCK=false;
  
  std::bitset<NumOfSegTFlag> trigger_flag;
  

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
      bool bh2_u_pass = false;
      bool bh2_d_pass = false;
      for(int ud=0; ud<nud[i]; ++ud){	  
        for(int idata=0; idata<ndata; ++idata){
          int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
          for( int m=0; m<nhit; ++m ){
            int val = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
            hid = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
            hptr_array[hid]->Fill(val);	    
            if(data[idata]==k_adc &&
               gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
              hid = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
              hptr_array[hid]->Fill(val);
            }
	    if(ud==0 && data[idata]==k_leading){
	      if(tdc_min < val && val < tdc_max){
		bh2_u_pass = true;
	      }
	    }
	    else if(ud==1 && data[idata]==k_leading){
	      if(tdc_min < val && val < tdc_max){
		bh2_d_pass = true;
	      }
	    }
          } // nhit
        } // data
      }//ud
      if(bh2_u_pass && bh2_d_pass){
	hid = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	hptr_array[hid]->Fill(seg);
	hit_seg_map[k_device].push_back(seg);
	mul++;
      }
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
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);
  } //hodo

    Bool_t has_hit_T2 = false;
  { // T2
    DetectorType kDET=kT2;
    const int k_device = gUnpacker.get_device_id("T2");
    const int k_adc = gUnpacker.get_data_id("T2", "adc");
    const int k_tdc = gUnpacker.get_data_id("T2", "leading");
    const UInt_t tdc_min = gUser.GetParameter("T2_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("T2_TDC", 1);
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
	      hit_seg_map[k_device].push_back(0);
	      has_hit_T2 = true;
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
	}
      }
    }//seg
    hid = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity[0]);
  } //hodo

      // SCH ------------------------------------------------------------
  {
    // data type    
    DetectorType kDET=kSCH;
    const int k_l=0;
    const int k_t=1;
    const int k_device = gUnpacker.get_device_id("SCH");
    const UInt_t tdc_min = gUser.GetParameter("SCH_TDC", 0);
    const UInt_t tdc_max = gUser.GetParameter("SCH_TDC", 1);
    int multiplicity=0;
    for(int seg = 0; seg<NumOfSegSCH; ++seg){
      int ntdc = 0;
      int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, k_l);
      if( nhit==0 ) continue;
      ++multiplicity;
      ntdc = 0;
      std::vector< int > leading_array;
      int leading_size = nhit;
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, seg, 0, k_l, m);
	hid=gHist.getSequentialID(kDET,0,kTDC,seg+1);
	hptr_array[hid]->Fill(tdc);
	leading_array.push_back(tdc);
	if (tdc_min < tdc && tdc < tdc_max) {
	  hid = gHist.getSequentialID(kDET,0,kHitPat,1);
	  hptr_array[hid]->Fill(seg);
	  ntdc++;
	}
      }	
      // traling
      nhit = gUnpacker.get_entries(k_device, 0, seg, 0, k_t);
      if( nhit==0 ) continue;
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, seg, 0, k_t , m);
	hid=gHist.getSequentialID(kDET,0,kTDC2D,seg+1);
	hptr_array[hid]->Fill(tdc);
	// for TOT
	if(m<leading_size){
	  hid=gHist.getSequentialID(kDET,0,kTOT,seg+1);
	  hptr_array[hid]->Fill(leading_array.at(m) - tdc);
	  hid=gHist.getSequentialID(kDET,0,kADC2D,seg+1);
	  hptr_array[hid]->Fill(leading_array.at(m), leading_array.at(m) - tdc);
	}
      }
      if(ntdc>0){
	hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	hptr_array[hid]->Fill(seg);
	hit_seg_map[k_device].push_back(seg);
	multiplicity++;
      }
    }//seg
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);	    
  } //SCH
  

  // TPC -----------------------------------------------------------

  {
    //if(cobo_data_size > 0){
    static const Int_t k_device   = gUnpacker.get_device_id( "TPC" );
    static const Int_t k_adc      = gUnpacker.get_data_id( "TPC", "adc" );
    // static const Int_t k_tdc_high = gUnpacker.get_data_id( "TPC", "tdc_high" );
    // static const Int_t k_tdc_low  = gUnpacker.get_data_id( "TPC", "tdc_low" );
    // sequential id
    static const Int_t tpca_id   = gHist.getSequentialID( kTPC, 0, kADC, 1 );
    static const Int_t tpcaget_a_id   = gHist.getSequentialID( kTPC, 1, kADC);
    static const Int_t tpcaget_rms_id   = gHist.getSequentialID( kTPC, 2, kADC);
    static const Int_t tpct_id   = gHist.getSequentialID( kTPC, 0, kTDC );
    static const Int_t rms_id    = gHist.getSequentialID( kTPC, 0, kPede );
    static const Int_t tpcfa_id   = gHist.getSequentialID( kTPC, 0, kFADC, 1 );
    static const Int_t tpca2d_id = gHist.getSequentialID( kTPC, 0, kADC2D );
    static const Int_t tpcmul_id = gHist.getSequentialID( kTPC, 0, kMulti );
    static const Int_t tpcbp_id   = gHist.getSequentialID( kTPC, 2, kTDC );
    static const Int_t agetmul_id = gHist.getSequentialID( kTPC, 3, kMulti );
    static const Int_t amulmax_id = gHist.getSequentialID( kTPC, 4, kMulti );
      
    hptr_array[tpca2d_id]->Reset();
    hptr_array[tpca2d_id+1]->Reset();
    hptr_array[tpca2d_id+2]->Reset();
    hptr_array[tpca2d_id+4]->Reset();
    hptr_array[agetmul_id]->Reset();
      
    // FADC before
    Int_t n_active_pad = 0;
    std::vector<Double_t> max_fadc( NumOfTimeBucket );
    
    std::vector<std::vector<std::vector<double>>> max_rms_fadc(NumOfAsadTPC,std::vector<std::vector<double>>(4,std::vector<double>(NumOfTimeBucket)));
    double aget_max_rms[NumOfAsadTPC][4]={0.};
    double aget_max_adc[NumOfAsadTPC][4]={0.};
    double aget_rms_mean[NumOfAsadTPC][4]={0.};
    int aget_count[NumOfAsadTPC][4] = {0};
    
    Int_t max_adc = -1;
    Int_t max_tb  = -1;
    Int_t max_pad = -1;
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
	
	std::vector<Double_t> fadc( nhit );
	  
	for( Int_t i=0; i<nhit; ++i ){
	  Int_t adc = gUnpacker.get( k_device, layer, 0, ch, k_adc, i );
	  fadc[i] = adc;
	  if( max_adc < adc && nhit == NumOfTimeBucket ){
	    max_adc = adc;
	    max_tb  = i;
	    max_pad = pad;
	  }
	}
	if( max_pad == pad ){
	  max_fadc = fadc;
	}
	
	Double_t mean = TMath::Mean( nhit, fadc.data() );
	Double_t rms = TMath::RMS( nhit, fadc.data() );
	Double_t max_adc = TMath::MaxElement( nhit, fadc.data() );
	Int_t loc_max = TMath::LocMax( nhit, fadc.data() );

	// Maximum RMS per AGET
	Int_t aget = gTpcPad.GetParam( layer, ch )->AGetId();
	Int_t asad = gTpcPad.GetParam( layer, ch )->AsAdId();
	aget_count[asad][aget]++;
	aget_max_adc[asad][aget]+=max_adc;
	aget_max_adc[asad][aget]/=(double)aget_count[asad][aget];
	aget_rms_mean[asad][aget]+=rms;
	aget_rms_mean[asad][aget]/=(double)aget_count[asad][aget];
	if(aget_max_rms[asad][aget] < rms){
	  aget_max_rms[asad][aget] = rms;
	  max_rms_fadc[asad][aget] = fadc;
	}
	
	

	if( max_adc - mean <= 0 ) continue;
	
	hptr_array[agetmul_id]->Fill( asad*4+aget );

	hptr_array[tpca_id]->Fill( max_adc - mean );
	
	hptr_array[tpcaget_a_id]->SetBinContent(asad*4+aget,aget_max_adc[asad][aget]);
	hptr_array[tpcaget_rms_id]->SetBinContent(asad*4+aget,aget_rms_mean[asad][aget]);

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
	  if( max_tb >= 60 && max_tb < 100
	      && pad_z >= -250 && pad_z <= -200
	      && max_adc - mean > 300
	      ){
	    hptr_array[tpcbp_id]->Fill( pad_x );
	  }
	}
      }
    }

    for(int n_asad = 0;n_asad<NumOfAsadTPC;n_asad++){
      for(int n_aget = 0;n_aget<4;n_aget++){
	hptr_array[tpcfa_id+n_asad*4+n_aget+1]->Reset();
	for(int i=0;i<NumOfTimeBucket;++i){
	  hptr_array[tpcfa_id+n_asad*4+n_aget+1]->SetBinContent(i+1,max_rms_fadc[n_asad][n_aget][i]+(3-n_aget)*700);
	}
      }
    }
    
    
    std::cout << "run# " << run_number << "  ev# " << event_number
	      << "  active pad = " << n_active_pad << std::endl;
    hptr_array[tpcmul_id]->Fill( n_active_pad );
    hptr_array[amulmax_id]->Fill( hptr_array[agetmul_id]->GetMaximum() );
    if( max_adc > 0 ){
      for( Int_t i=0; i<NumOfTimeBucket; ++i ){
	hptr_array[tpcfa_id]->Fill( i, max_fadc[i] );
      }
    }
    // TDC (Time Stamp)
    // UInt_t tdc_h = gUnpacker.get( k_device, 0, 0, 0, k_tdc_high );
    // UInt_t tdc_l = gUnpacker.get( k_device, 0, 0, 0, k_tdc_low );
    // std::cout << tdc_h << ", " << tdc_l << std::endl;

#if 0
    // Debug, dump data relating this detector
    gUnpacker.dump_data_device(k_device);
#endif
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

  //BcOut
  {
    
    RawData rawData;
    for (const auto& name : DCNameList.at("BcOut")) rawData.DecodeHits(name);
    EventAnalyzer evAna;
    evAna.DCRawHit("BcOut", rawData);

    DCAnalyzer dcAna(rawData);
    dcAna.DecodeBcOutHits();
    dcAna.TotCut("BLC2a");
    dcAna.TotCut("BLC2b");
    dcAna.DriftTimeCut("BLC2a");
    dcAna.DriftTimeCut("BLC2b");
    evAna.DCHit("BcOut", dcAna);
    dcAna.TrackSearchBcOut();
    evAna.BcOutTracking(dcAna);

    double ff_offset = 150.; //TPC center : 0, BcOut FF(z=0) : 150.
    double z_start = -900.;
    double z_end = -150.-143.;
    int bcout_count = 5000;
    double step = (z_end-z_start) / (double)bcout_count;
      

    int hist_id = gHist.getSequentialID(kEventDisplay, 0, kHitPoly, 13);
    auto h2 = dynamic_cast<TH2*>(hptr_array[hist_id]);
    h2->Reset();

    
    for(const auto& track : dcAna.GetBcOutTrackContainer()){
      double x0 = track->GetX0();
      double u0 = track->GetU0();
      for(int n_bcout = 0;n_bcout<bcout_count;n_bcout++){
	double real_z = z_start+step*n_bcout;
	h2->Fill(real_z+ff_offset,x0+real_z*u0,100);
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

