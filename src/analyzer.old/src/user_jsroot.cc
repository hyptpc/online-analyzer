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
#include <TVector3.h>

#include <TEveManager.h>
#include <TEveBox.h>

#include <user_analyzer.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <DAQNode.hh>
#include <filesystem_util.hh>

#include "AftHelper.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "EMCParamMan.hh"
#include "EventAnalyzer.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "HistMaker.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "HodoRawHit.hh"
#include "HttpServer.hh"
#include "MacroBuilder.hh"
#include "MatrixParamMan.hh"
#include "MsTParamMan.hh"
#include "TpcPadHelper.hh"
#include "UserParamMan.hh"

#define DEBUG    0
#define FLAG_DAQ 1

namespace
{
using hddaq::unpacker::GConfig;
using hddaq::unpacker::GUnpacker;
using hddaq::unpacker::DAQNode;
std::vector<TH1*> hptr_array;
const auto& gUnpacker = GUnpacker::get_instance();
auto&       gHist     = HistMaker::getInstance();
auto&       gHttp     = HttpServer::GetInstance();
auto&       gMatrix   = MatrixParamMan::GetInstance();
const auto& gAftHelper = AftHelper::GetInstance();
auto&       gMsT      = MsTParamMan::GetInstance();
const auto& gUser     = UserParamMan::GetInstance();
}

namespace analyzer
{

//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  // gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1110);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  // gStyle->SetStatW(.42);
  // gStyle->SetStatH(.35);
  gStyle->SetStatW(.32);
  gStyle->SetStatH(.25);
  // gStyle->SetPalette(55);

  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  // gConfMan.InitializeParameter<HodoParamMan>("HDPRM");
  // gConfMan.InitializeParameter<HodoPHCMan>("HDPHC");
  // gConfMan.InitializeParameter<DCGeomMan>("DCGEO");
  // gConfMan.InitializeParameter<DCTdcCalibMan>("DCTDC");
  // gConfMan.InitializeParameter<DCDriftParamMan>("DCDRFT");
  // gConfMan.InitializeParameter<TpcPadHelper>("TPCPAD");
  gConfMan.InitializeParameter<UserParamMan>("USER");
  if(!gConfMan.IsGood()) return -1;

  Int_t port = 9090;
  if(argv.size()==4)
    port = TString(argv[3]).Atoi();

  gHttp.SetPort(port);
  gHttp.Open();
  gHttp.CreateItem("/", "Online Analyzer");
  gHttp.CreateItem("/Tag", "Tag Check");
  gHttp.SetItemField("/Tag", "_kind", "Text");
  std::stringstream ss;
  ss << "<div style='color: white; background-color: black;"
     << "width: 100%; height: 100%;'>"
     << "Tag cheker is running" << "</div>";
  gHttp.SetItemField("/Tag", "value", ss.str().c_str());
  gHttp.Register(gHist.createHODO());
  // gHttp.Register(gHist.createMatrix());

  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }

  //___ Macro for HttpServer
  gHttp.Register(http::CaenV792_00_15());
  gHttp.Register(http::CaenV792_16_31());
  gHttp.Register(http::HulHrTdc_00_15());
  gHttp.Register(http::HulHrTdc_16_31());

  for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
    hptr_array[i]->SetDirectory(0);
  }

  return 0;
}

//____________________________________________________________________________
int
process_end(void)
{
  hptr_array.clear();
  return 0;
}

//____________________________________________________________________________
int
process_event(void)
{
  static Int_t run_number = -1;
  {
    if(run_number != gUnpacker.get_root()->get_run_number()){
      for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
	hptr_array[i]->Reset();
      }
      run_number = gUnpacker.get_root()->get_run_number();
    }
  }
  auto event_number = gUnpacker.get_event_number();

  { ///// Tag Checker
    static const auto& gConfig = GConfig::get_instance();
    static const TString tout(gConfig.get_control_param("tout"));
    std::stringstream ss;
    ss << "<div style='color: white; background-color: black;"
       << "width: 100%; height: 100%;'>";
    ss << "RUN " << run_number << "   Event " << event_number
       << "<br>";
    if(!gUnpacker.is_good()){
      std::ifstream ifs(tout);
      if(ifs.good()){
	TString buf;
	while(!ifs.eof()){
	  buf.ReadLine(ifs, false);
	  if(buf.Contains("!") && !buf.Contains("............!"))
	    buf = "<font color='yellow'>" + buf + "</font>";
	  // if(buf.Contains("bUuSELselthdg")){
	  //   ss << TString(' ', 24);
	  ss << buf << "<br>";
	}
	ifs.close();
	hddaq::tag_summary.seekp(0, std::ios_base::beg);
      }else{
	ss << Form("Failed to read %s", tout.Data());
      }
      ss << "</div>";
      gHttp.SetItemField("/Tag", "value", ss.str().c_str());
    }
  }

  { ///// HODO
    static const auto device_id = gUnpacker.get_device_id("HODO");
    static const auto adc_id = gUnpacker.get_data_id("HODO", "adc");
    static const auto tdc_id = gUnpacker.get_data_id("HODO", "tdc");
    static const auto adc_hid = gHist.getSequentialID(kHODO, 0, kADC);
    static const auto tdc_hid = gHist.getSequentialID(kHODO, 0, kTDC);
    static const auto awt_hid = gHist.getSequentialID(kHODO, 0, kADCwTDC);
    static const Int_t n_seg = 32;
    for(Int_t seg=0; seg<n_seg; ++seg){
      auto nhit = gUnpacker.get_entries(device_id, 0, seg, 0, adc_id);
      UInt_t adc = 0;
      if (nhit != 0) {
	adc = gUnpacker.get(device_id, 0, seg, 0, adc_id);
	hptr_array[adc_hid + seg]->Fill(adc);
      }
      Bool_t hit_flag = false;
      for(Int_t m=0, n=gUnpacker.get_entries(device_id, 0, seg, 0, tdc_id);
	  m<n; ++m) {
	auto tdc = gUnpacker.get(device_id, 0, seg, 0, tdc_id, m);
	if (tdc != 0) {
	  hptr_array[tdc_hid + seg]->Fill(tdc);
	  // ADC wTDC
	  hit_flag = true;
	}
	if (hit_flag) {
	  hptr_array[awt_hid + seg]->Fill(adc);
	}
      }
    }
#if 0
    gUnpacker.dump_data_device(k_device);
#endif
  }

  gSystem->ProcessEvents();
  return 0;
}

}
