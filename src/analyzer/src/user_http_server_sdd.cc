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
  for(int i=0;(int)i<argv.size();i++){
    std::cout<<i<<"  "<<argv[i]<<std::endl;
  }
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
  
  gHttp.SetPort(8082);
  gHttp.Open();
  gHttp.Register(gHist.createQDC(kSDD,"SDD",32,4096,0,4096));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::QDC(kSDD,"SDD",0, 32,8,4,0,4096));
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  
  return 0;
}
  
  //____________________________________________________________________________
  int
  process_end( void )
  {
    TFile *f=new TFile("hist.root","recreate"); 
    for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
      hptr_array[i]->Write();
    }
    f->Close();
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

  //  gUnpacker.dump_data_fe(513);  
  // Counter ------------------------------------------------------------
  {
    // data type
    static const int k_device = gUnpacker.get_device_id("SDD");
    static const int NumOfSegments=32;
    static const int qdc1_id  = gHist.getSequentialID(kSDD, 0, kADC, 1);
    for(int l = 0; l<NumOfSegments; ++l){
      int nhit = gUnpacker.get_entries(k_device, 0, l, 0, 0);
      if( nhit ){
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, l, 0, 0 , m);
	  hptr_array[qdc1_id + l]->Fill(tdc);
	}
      }
    }    
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  
  return 0;
}
  
}

