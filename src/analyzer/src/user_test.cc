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

    // hodoscopes
    TString hodo_name="Counter";
    const int nsegs=64;
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};

    int nbinsqdc=1024*2;
    int nbinshrtdc=int(1e4);
    double hrtdcmin=0.;
    double hrtdcmax=2e6;    

    TTree* tree;
    int v_adc[16];    
    int v_leading[16][4];
    int v_trailing[16][4];
    int v_leading2[16][4];
    int v_trailing2[16][4];

    TString outputname="test.root";
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
    // Hodoscopes
    DetectorType kDET=kDetectorZero;
    gHttp.Register(gHist.createHodo(kDET ,"Counter"  ,64, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
    
    if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
    
    // Macro for HttpServer
    // Hodoscopes
    gHttp.Register(http::QDC(kDET,"_Trigger",0,0,8,4,2,0,1000),"");
    gHttp.Register(http::QDC(kDET,"_barrelTOF",0,8,8,4,2,0,1000),"");
    gHttp.Register(http::QDC(kDET,"_aaaa",1,2,9,3,3,100,1100),"");
    gHttp.Register(http::MHTDCTDC(kDET,"_Counter",0,16,4,4)     ,"");
    gHttp.Register(http::MHTDCHitPatMulti(kDET,"_Counter")      ,"");

    //=== set directory ===//
    for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
      hptr_array[i]->SetDirectory(0);
    }
    //=== set directory ===//
    
    tree=new TTree("raw","rawdata");
    tree->Branch("adc",v_adc,"adc[16]/I");
    tree->Branch("leading",v_leading,"leading[16][4]/I");
    tree->Branch("trailing",v_trailing,"trailing[16][4]/I");
    tree->Branch("leading2",v_leading2,"leading2[16][4]/I");
    tree->Branch("trailing2",v_trailing2,"trailing2[16][4]/I");
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
    for(int i=0;i<16;i++){
      v_adc[i]=-1;
      for(int j=0;j<4;j++){
	v_leading[i][j]=-1;
	v_trailing[i][j]=-1;
	v_leading2[i][j]=-1;
	v_trailing2[i][j]=-1;
      }
    }
    
    // Hodoscope ------------------------------------------------------------
    int hid;
    {
    // data type    
      DetectorType kDET=kDetectorZero;
      const int k_device = gUnpacker.get_device_id(hodo_name.Data());
      int multiplicity=0;
      int mul[2]={0,0};
      for(int seg = 0; seg<nsegs; ++seg){
      int ntdc[2]={0,0};
      int ud=0;
      for(int idata=0; idata<ndata; ++idata){
	int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	if(seg>15){
	  if(data[idata]==k_adc) continue;
	  //	  std::cout<< hodo_name.Data() << " seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	  //	  std::cout<< hodo_name.Data() << " seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" , nhit: "<<nhit<<std::endl;
	}
	if(data[idata]==k_leading&&nhit>0){
	  hid  = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
	  hptr_array[hid]->Fill(seg);	    
	  mul[ud]++;
	  ntdc[ud]=nhit;	  
	}
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	  
	  hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	  hptr_array[hid]->Fill(tdc);	    	 
	  
	  if(data[idata]==k_adc){
	    v_adc[seg]=tdc;
	    if(gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
	      hid  = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
	      hptr_array[hid]->Fill(tdc);	  	      
	    } 
	  }else if(data[idata]==k_leading){
	    if(seg<16){
	      v_leading[seg][m]=tdc;
	    }
	    if(seg>47){
	      v_leading2[seg-48][m]=tdc;
	    }
	  }else if(data[idata]==k_trailing){
	    if(seg<16){
	      v_trailing[seg][m]=tdc;
	    }
	    if(seg>47){
	      v_trailing2[seg-48][m]=tdc;
	    }
	  }
	} // nhit
      } // data
      }//seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(mul[0]);
    }
    tree->Fill();
    // Scaler
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
    
    return 0;
  }
  
}

