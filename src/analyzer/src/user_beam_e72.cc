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
#include "UnpackerConfig.hh"
#include "UnpackerManager.hh"
#include "UnpackerXMLReadDigit.hh"
#include "DAQNode.hh"
#include "filesystem_util.hh"
#include "RawData.hh"

#include "ConfMan.hh"
#include "HistMaker.hh"
#include "DetectorID.hh"
#include "PsMaker.hh"
#include "MacroBuilder.hh"
#include "UserParamMan.hh"
#include "HodoParamMan.hh"
#include "DCGeomMan.hh"
#include "DetSizeMan.hh"
#include "DCRawHit.hh"
#include "DCAnalyzer.hh"
#include "DCLocalTrack.hh"
#include "EventAnalyzer.hh"

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
const auto& gUConf   = hddaq::unpacker::GConfig::get_instance();

TText text;
TText end;
TString outputname="tmp.root";

}

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{
  ConfMan::getInstance().initialize(argv);
  ConfMan::getInstance().initializeUserParamMan();
  ConfMan::getInstance().initializeHodoParamMan();
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
  int port=8085;
  
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

  gHttp.Register(gHist.createBcInTracking(true));
  gHttp.Register(gHist.createBcOutTracking(true));

  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }

  //=== set directory ===//
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }

  gHttp.Begin();
  return 0;
  
}

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
  static Int_t run_number = -1;
  if(run_number != gUnpacker.get_run_number()){
    for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
      hptr_array[i]->Reset();
    }
    run_number = gUnpacker.get_run_number();
  }
  auto event_number = gUnpacker.get_event_number();

  for (auto& h : hptr_array){
    h->SetTitle(h->GetName() + TString(Form(" run%05d", run_number)));
  }

  RawData rawData;
  for (const auto& name : DCNameList.at("BcIn")) rawData.DecodeHits(name);
  for (const auto& name : DCNameList.at("BcOut")) rawData.DecodeHits(name);

  EventAnalyzer evAna;
  evAna.DCRawHit("BcIn", rawData);
  evAna.DCRawHit("BcOut", rawData);
  
  DCAnalyzer dcAna(rawData);
  
  //BcInTracking
  dcAna.DecodeBcInHits();
  dcAna.TotCut("BLC1a");
  dcAna.TotCut("BLC1b");
  dcAna.DriftTimeCut("BLC1a");
  dcAna.DriftTimeCut("BLC1b");
  evAna.DCHit("BcIn", dcAna);
  dcAna.TrackSearchBcIn();
  evAna.BcInTracking(dcAna);

  //BcOutTracking
  dcAna.DecodeBcOutHits();
  dcAna.TotCut("BLC2a");
  dcAna.TotCut("BLC2b");
  dcAna.DriftTimeCut("BLC2a");
  dcAna.DriftTimeCut("BLC2b");
  evAna.DCHit("BcOut", dcAna);
  dcAna.TrackSearchBcOut();
  evAna.BcOutTracking(dcAna);

  int hist_bcin_x_id = gHist.getSequentialID(kBcInTracking,0,0,1);
  int hist_bcin_y_id = gHist.getSequentialID(kBcInTracking,0,0,2);
  int hist_bcin_2d_id = gHist.getSequentialID(kBcInTracking,0,0,3);
  for(const auto& track : dcAna.GetBcInTrackContainer()){
    hptr_array[hist_bcin_x_id]->Fill(track->GetX0());
    hptr_array[hist_bcin_y_id]->Fill(track->GetY0());
    hptr_array[hist_bcin_2d_id]->Fill(track->GetX0(),track->GetY0());
  }
  
  int hist_bcout_x_id = gHist.getSequentialID(kBcOutTracking,0,0,1);
  int hist_bcout_y_id = gHist.getSequentialID(kBcOutTracking,0,0,2);
  int hist_bcout_2d_id = gHist.getSequentialID(kBcOutTracking,0,0,3);
  for(const auto& track : dcAna.GetBcOutTrackContainer()){
    hptr_array[hist_bcout_x_id]->Fill(track->GetX0());
    hptr_array[hist_bcout_y_id]->Fill(track->GetY0());
    hptr_array[hist_bcout_2d_id]->Fill(track->GetX0(),track->GetY0());
  }
  
  gSystem->ProcessEvents();
  
  return 0;
}
}

