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
    TTree* tree;
    int v_adc[16];
    int v_tdc[16];    
    int v_leading[32][4];
    int v_trailing[32][4];
    TString outputname="tmp.root";
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
  gHttp.SetPort(port);
  gHttp.Open();
  std::cout<<"QDC1:"<<kQDC1<<std::endl;
  gHttp.Register(gHist.createQDC(kQDC1,"QDC1",32,1024,0,4096));
  std::cout<<"QDC2:"<<kQDC2<<std::endl;
  gHttp.Register(gHist.createQDC(kQDC2,"TDC1",32,1024,0,4096));
  std::cout<<"HRTDC:"<<kHRTDC<<std::endl;
  gHttp.Register(gHist.createMHTDC(kHRTDC,"HRTDC",64,6553600/4.));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  std::cout<<"gHttp.Register: QDC1"<<std::endl;
  gHttp.Register(http::QDC(kQDC1,"QDC1",0,0,16,4,4,0,4096));
  std::cout<<"gHttp.Register: QDC2"<<std::endl;
  gHttp.Register(http::QDC(kQDC2,"TDC1",0,0,16,4,4,0,4096));
  std::cout<<"gHttp.Register: HRTDC"<<std::endl;
  gHttp.Register(http::MHTDCTDC(kHRTDC,"HRTDC",0,32,8,4));
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }

  tree=new TTree("raw","rawdata");
  tree->Branch("adc",v_adc,"adc[16]/I");
  tree->Branch("tdc",v_tdc,"tdc[16]/I");
  tree->Branch("leading",v_leading,"leading[32][4]/I");
  tree->Branch("trailing",v_trailing,"trailing[32][4]/I");
  tree->SetDirectory(0);
  return 0;
}
  
  //____________________________________________________________________________
  int
  process_end( void )
  {
    TFile *root_file=new TFile(outputname,"recreate");
    for(int i=0;i<hptr_array.size();i++)
      if(hptr_array[i]) hptr_array[i]->Write();
    tree->Print();
    tree->Write();
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
  for(int i=0;i<16;i++){
    v_adc[i]=-1;
    v_tdc[i]=-1;
  }
  for(int i=0;i<32;i++){
    for(int j=0;j<4;j++){
      v_leading[i][j]=-1;
      v_trailing[i][j]=-1;
    }
  }

  // Counter EMC ------------------------------------------------------------
  if(0){
    const int k_device = gUnpacker.get_device_id("EMC");
    for(int ch=0;ch<6;ch++){	
      int nhit = gUnpacker.get_entries(k_device, 0, 0, 0, ch);
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0, 0, ch, m);
	//	std::cout<<"EMC:  "<<ch<<"  "<<m<<"  "<<tdc<<std::endl;
      }
    }
  }
  // Counter TDC ------------------------------------------------------------
  {
    const int nmtdc=1;
    DetectorType mtdc[nmtdc]={kHRTDC};
    TString mtdc_name[nmtdc]={"HRTDC1"};
    const int nsegs[nmtdc]={64};
    for(int i=0;i<nmtdc;i++)
      {
	// data type
	DetectorType kDET=mtdc[i];
	const int k_device = gUnpacker.get_device_id(mtdc_name[i].Data());
	const int k_l    = 0;
	const int k_t    = 1;      
	// TDC & HitPat & Multi
	int hid=-1;
	int multiplicity    = 0;
	for(int l = 0; l<nsegs[i]; ++l){
	  int nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_l);
	  if( nhit ){
	    ++multiplicity;
	    hid=gHist.getSequentialID(kDET,0,kHitPat,1);
	    hptr_array[hid]->Fill(l, nhit);
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, l, 0, k_l, m);
	      hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	      hptr_array[hid]->Fill(tdc);
	      if(m<4&&l<32) v_leading[l][m]=tdc;
	    }	
	  }
	  nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_t);
	  if( nhit ){
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, l, 0, k_t, m);
	      hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	      hptr_array[hid]->Fill(tdc);
	      if(m<4&&l<32) v_trailing[l][m]=tdc;
	    } 	
	  }
	}      
	hid=gHist.getSequentialID(kDET,0,kMulti,1);
	hptr_array[hid]->Fill(multiplicity);
      }
  }

  // QDC&PADC ------------------------------------------------------------
  {
    // data type
    const int nadc=2;
    DetectorType adc[nadc]={kQDC2,kQDC1};
    TString adc_name[nadc]={"TDC1","QDC1"};
    const int nsegs[nadc]={32,32};
    for(int i=0;i<nadc;i++){
      DetectorType kDET=adc[i];
      const int k_device = gUnpacker.get_device_id(adc_name[i].Data());
      for(int l = 0; l<nsegs[i]; ++l){
	int nhit = gUnpacker.get_entries(k_device, 0, l, 0, 0);
	//	std::cout<<adc_name[i]<<"  "<<nsegs[i]<<"  "<<nhit<<std::endl;
	if( nhit ){
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, 0, l, 0, 0 , m);
	    if(i==1&&l<16) v_adc[l]=tdc;
	    if(i==0&&l<16) v_tdc[l]=tdc;
	    int hid  = gHist.getSequentialID(kDET, 0, kADC, l+1);
	    hptr_array[hid]->Fill(tdc);
	    if(v_tdc[l]>0 && i==1){
	      hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, l+1);
	      hptr_array[hid]->Fill(tdc);	  		
	    }
	  }
	}
      }    
    }
  }
  tree->Fill();
  //Scaler
  // const int k_device = gUnpacker.get_device_id("Scaler");	
  // int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
  // for(int module=0;module<1;module
  // for( int m=0; m<nhit; ++m ){
  //   if
  //   int tdc = gUnpacker.get(k_device, module, 0, channel);
  // }	    
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  //  sleep(1);
  //  gUnpacker.dump_data_fe(513);
  // gUnpacker.dump_data_fe(1677);
  gSystem->ProcessEvents();
  return 0;
}
  
}
