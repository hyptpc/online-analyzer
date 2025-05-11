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
  gHttp.Register(gHist.createSDD(kSDD ,"SDD",2,4));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }  
  // Macro for HttpServer
  for(int i=0;i<2;i++){
    gHttp.Register(http::SDDADC(kSDD,Form("_SDD_%d",i),i,0,4000,1),"");
  }
  gHttp.Register(http::SDDMHTDC(kSDD,"_SDDgate" ,0,2,4,2),"");
  gHttp.Register(http::SDDMHTDC(kSDD,"_SDDreset",1,2,4,2),"");
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  gHttp.Begin();
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

  int hid;
  //  gUnpacker.dump_data_fe(513);  
  // SDD
  {
    const int k_adc=0;
    const int k_leading=8;
    const int k_trailing=9;
    const int k_resetl=11;
    const int k_resett=12;
    const int ndata=5;
    const int data[ndata]={k_leading,k_trailing,k_resetl,k_resett,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kResetL,kResetT,kADC};    
    DetectorType kDET=kSDD;
    TString tmpname="SDD";
    const int k_device = gUnpacker.get_device_id(tmpname.Data());
    int multiplicity=0;
    int nports=2;
    int nunits=4;
    for(int iport=0;iport<nports;iport++){
      for(int iunit=0;iunit<nunits;iunit++){
	int index=iunit+iport*nunits;
	for(int idata=0; idata<ndata; ++idata){
	  int nhit = gUnpacker.get_entries(k_device, iport, 0, iunit, data[idata]);
	  if(data[idata]==k_leading&&nhit>0){
	    hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
	    hptr_array[hid]->Fill(index);	    
	    multiplicity++;
	  }
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, iport, 0, iunit, data[idata] , m);
	    if(data[idata]==k_adc){
	      for(int i=0;i<8;i++){
		tdc = gUnpacker.get(k_device, iport, 0, iunit, data[idata]+i , m);
		hid  = gHist.getSequentialID(kDET, index, type[idata], i+1);     
		hptr_array[hid]->Fill(tdc);
		if(gUnpacker.get_entries(k_device, iport, 0, iunit, k_leading) ){
		  hid  = gHist.getSequentialID(kDET, index, kADCwTDC, i+1);
		  hptr_array[hid]->Fill(tdc);	  				
		}
	      }
	    }else{
	      hid  = gHist.getSequentialID(kDET, index, type[idata], 1);
	      hptr_array[hid]->Fill(tdc);	    
	    }
	  } // nhit
	} // data
      }// unit
    }// port
    hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
    hptr_array[hid]->Fill(multiplicity);	    
  } // hodox
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  
  return 0;
}
  
}

