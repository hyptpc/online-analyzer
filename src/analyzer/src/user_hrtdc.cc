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

#define DEBUG    0
#define FLAG_DAQ 1
#define WIRE 1 // for Canvas only

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
    TString outputname="tmp.root";

    int nbinsqdc=512;
    int nbinshrtdc=int(1e4);
    double hrtdcmin=0.;
    double hrtdcmax=2e6;    
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
  int port=8085;
  if(argv.size()==4){
    outputname=argv.at(3);
    port=8089;
  }
  gHttp.SetPort(port);
  gHttp.Open();

  // Hodoscopes
  //  gHttp.Register(gHist.createHodo(kBHD		,"BHT"  ,16, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC1	,"HRTDC1"   ,64, 1, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC2	,"HRTDC2"   ,64, 1, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC3	,"HRTDC3"   ,64, 1, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC4	,"HRTDC4"   ,64, 1, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC5	,"HRTDC5"   ,64, 1, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kHRTDC6	,"HRTDC6"   ,64, 1, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::MHTDCTDC(kHRTDC1,"_HRTDC1",0,64,8,8),"");
  gHttp.Register(http::MHTDCTDC(kHRTDC2,"_HRTDC2",0,64,8,8),"");
  gHttp.Register(http::MHTDCTDC(kHRTDC3,"_HRTDC3",0,64,8,8),"");
  gHttp.Register(http::MHTDCTDC(kHRTDC4,"_HRTDC4",0,64,8,8),"");
  gHttp.Register(http::MHTDCTDC(kHRTDC5,"_HRTDC5",0,64,8,8),"");
  gHttp.Register(http::MHTDCTDC(kHRTDC6,"_HRTDC6",0,64,8,8),"");

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
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  static HistMaker&     gHist     = HistMaker::getInstance();
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  int hid;
  double time0=-1;
  // Hodoscope ------------------------------------------------------------
  {
    // data type    
    DetectorType kDET[6]={kHRTDC1,kHRTDC2,kHRTDC3,kHRTDC4,kHRTDC5,kHRTDC6};
    const int ndata=2;
    const int k_leading=0;
    const int k_trailing=1;    
    const int data[ndata]={k_leading,k_trailing};
    const int type[ndata]={kTDC,kTDC2D};

    for(int i=0;i<6;i++){
      TString tmpname=Form("HRTDC%d",i+1);
      const int k_device = gUnpacker.get_device_id(tmpname.Data());
      for(int seg = 0; seg<64; ++seg){
	for(int idata=0; idata<ndata; ++idata){
	  int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, 0, seg, 0, data[idata] , m);
	    hid  = gHist.getSequentialID(kDET[i], 0, type[idata], seg+1);
	    hptr_array[hid]->Fill(tdc);	    
	  } // nhit
	} // nhit
      }//seg
    } //hodo
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  gSystem->ProcessEvents();
  return 0;
}

}

