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
#include "MacroBuilder.hh"

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

  gHttp.Register(gHist.createBLDC(kBPC,"BPC",8,32,true));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer (TCanvas)
  // for BPC
  gHttp.Register(http::BLDCHitPat(kBPC,"_BPC",32));
  gHttp.Register(http::BLDCMulti(kBPC,"_BPC",32));
  for(int l=0;l<8;l++){ 
    gHttp.Register(http::BLDCWIRE(kBPC,"BPC",l,32,6,6), "BPC");
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
  // Chamber ------------------------------------------------------------
  const int nchm=1;
  DetectorType chm[nchm]={kBPC};
  TString chm_name[nchm]={"BPC"};
  const int nwires[nchm]={32};
  const int nlayers[nchm]={8};
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
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  return 0;
}
  
}

