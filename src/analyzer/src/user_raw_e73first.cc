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

    // chambers
    const int nchm=6;
    DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b, kBPC,kCDC};
    TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b", "BPC", "CDC"};
    const int nwires[nchm]={32,32,32,32,32,16};
    const int nlayers[nchm]={8,8,8,8,8,118};
    // for CDC analysis
    std::map< int, int > CDCmap;
   
    // hodoscopes
    const int nhodo=9;
    DetectorType hodo[nhodo]={kBHD,kT0,kT0new,kDEF,kCDH,kVeto,kPbF2,kLeak,kRC};
    TString hodo_name[nhodo]={"BHD","T0","T0new","DEF","CDH","Veto","PbF2","Leak","RC"};
    const int nsegs[nhodo]={16,5,1,5,36,4,40,6,8};
    const int nud[nhodo]={2,2,2,2,2,2,1,1,2};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};

    int nbinsqdc=512;
    int nbinshrtdc=int(1e4);
    double hrtdcmin=0.;
    double hrtdcmax=2e6;    

    TString flagnames[16]={"SpillStart","SpillEnd",
			   "Beam","Pion","Kaon2","Kaon3",
			   "KxCDH1","KxCDH2","KxCDH3","KxCDH1xG","KaonxG",
			   "PixCDH","PixPbF2","ExPbF2",
			   "CDH cosmic","clock(10s)"};
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
  gHttp.Register(gHist.createHodo(kBHD		,"BHD"  ,16, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kT0		,"T0"   , 5, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kT0new	,"T0new", 1, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kAC		,"AC"   , 5, 1, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kDEF		,"DEF"  , 5, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kCDH		,"CDH"  ,36, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kVeto		,"Veto" , 4, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kPbF2		,"PbF2" ,40, 1, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax)); 
  gHttp.Register(gHist.createHodo(kLeak		,"Leak" , 6, 1, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kRC		,"RC"   , 8, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createMHTDC(kTriggerFlag ,"TriggerFlag",16,2000));
  // Chambers
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBPC,"BPC",8,32,true));
  gHttp.Register(gHist.createBLDC(kCDC,"CDC",118,16,true));
  gHttp.Register(gHist.createCDC(true));

  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  // Hodoscopes
  gHttp.Register(http::QDC(kBHD,"_BHD_U",0,0,16,4,4,100,500),"BHD");
  gHttp.Register(http::QDC(kBHD,"_BHD_D",1,0,16,4,4,100,500),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_U",0,16,4,4),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_D",1,16,4,4),"BHD");
  gHttp.Register(http::MHTDCHitPatMulti(kBHD,"_BHD"),"BHD");

  gHttp.Register(http::QDC(kT0,"_T0_U",0,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");
  gHttp.Register(http::MHTDCHitPatMulti(kT0,"_T0"),"T0");

  gHttp.Register(http::QDC(kT0new,"_T0new_U",0,0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::QDC(kT0new,"_T0new_D",1,0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::MHTDCTDC(kT0new,"_T0new_U",0,1,1,1),"T0new");
  gHttp.Register(http::MHTDCTDC(kT0new,"_T0new_D",1,1,1,1),"T0new");
  gHttp.Register(http::MHTDCHitPatMulti(kT0new,"_T0new"),"T0new");

  gHttp.Register(http::QDC(kAC,"_AC",0,0,5,3,2,0,2000),"AC");
  gHttp.Register(http::MHTDCTDC(kAC,"_AC",0,1,1,1),"AC");

  gHttp.Register(http::QDC(kDEF,"_DEF_U",0,0,5,3,2,0,2000),"DEF");
  gHttp.Register(http::QDC(kDEF,"_DEF_D",1,0,5,3,2,0,2000),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_U",0,5,3,2),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_D",1,5,3,2),"DEF");
  gHttp.Register(http::MHTDCHitPatMulti(kDEF,"_DEF"),"DEF");

  gHttp.Register(http::QDC(kCDH,"_CDH_U",0,0,36,6,6,0,1000),"CDH");
  gHttp.Register(http::QDC(kCDH,"_CDH_D",1,0,36,6,6,0,1000),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_U",0,36,6,6),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_D",1,36,6,6),"CDH");
  gHttp.Register(http::MHTDCHitPatMulti(kCDH,"_CDH"),"CDH");

  gHttp.Register(http::QDC(kVeto,"_Veto_L",0,0,4,2,2,0,1000),"Veto");
  gHttp.Register(http::QDC(kVeto,"_Veto_R",1,0,2,2,2,0,1000),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_U",0,4,2,2),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_D",1,2,2,2),"Veto");
  gHttp.Register(http::MHTDCHitPatMulti(kVeto,"_Veto"),"Veto");

  gHttp.Register(http::QDC(kPbF2,"_PbF2",0,0,40,8,5,0,500),"PbF2");
  gHttp.Register(http::MHTDCTDC(kPbF2,"_PbF2",0,40,8,5),"PbF2");
  gHttp.Register(http::MHTDCHitPatMulti(kPbF2,"_PbF2"),"PbF2");

  gHttp.Register(http::QDC(kLeak,"_Leak",0,0,6,3,2,0,500),"Leak");
  gHttp.Register(http::MHTDCTDC(kLeak,"_Leak",0,6,3,2),"Leak");
  gHttp.Register(http::MHTDCHitPatMulti(kLeak,"_Leak"),"Leak");

  gHttp.Register(http::QDC(kRC,"_RC_L",0,0,8,4,2,0,1000),"RC");
  gHttp.Register(http::QDC(kRC,"_RC_R",1,0,8,4,2,0,1000),"RC");
  gHttp.Register(http::MHTDCTDC(kRC,"_RC_U",0,8,4,2),"RC");
  gHttp.Register(http::MHTDCTDC(kRC,"_RC_D",1,8,4,2),"RC");
  gHttp.Register(http::MHTDCHitPatMulti(kRC,"_RC"),"RC");

  // Chambers
  for(int i=0;i<nchm-1;i++){
    const std::string tmpstr=Form("_%s",chm_name[i].Data());
    gHttp.Register(http::BLDCHitPat(chm[i],tmpstr,8),chm_name[i]); 
    gHttp.Register(http::BLDCMulti(chm[i],tmpstr,8),chm_name[i]);
    gHttp.Register(http::BLDCTDC(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTOT(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTDCvsTOT(chm[i],tmpstr,8), chm_name[i]);
    for(int l=0;l<8;l++){ 
      gHttp.Register(http::BLDCWIRE(chm[i],tmpstr,l,32,8,4), chm_name[i]);
    }
  }
  // TriggerFlag
  gHttp.Register(http::MHTDCTDC(kTriggerFlag,"_TriggerFlag",0,16,4,4),"TriggerFlag");
  gHttp.Register(http::MHTDCHitPatMulti(kTriggerFlag,"_TriggerFlag"),"TriggerFlag");
  {
    int hid1 = gHist.getSequentialID(kTriggerFlag, 0, kHitPat, 1);  
    hptr_array[hid1]->GetXaxis()->SetTitle("");
    for( Int_t i=0; i<16; ++i ){
      int hid2 = gHist.getSequentialID(kTriggerFlag, 0, kTDC, i+1);
      hptr_array[hid2]->SetTitle(Form("%s_%s",hptr_array[hid2]->GetTitle(),flagnames[i].Data()));
      hptr_array[hid1]->GetXaxis()->SetBinLabel(i+1,flagnames[i]);
    }
  }
		   
  // for CDC analysis
  gHttp.Register(http::BLDCHitPat(kCDC2,"CDC2",15,0,5,3), "CDC2");
  gHttp.Register(http::BLDCMulti(kCDC2,"CDC2",15,0,5,3), "CDC2");
  gHttp.Register(http::BLDCTDC(kCDC2,"CDC2",15,0,5,3), "CDC2");
  gHttp.Register(http::BLDCTOT(kCDC2,"CDC2",15,0,5,3), "CDC2");
  gHttp.Register(http::BLDCTDCvsTOT(kCDC2,"CDC2",15,0,5,3), "CDC2");
  for(int l=0; l<8; l++){
    int n=15;
    if(l==7) n=13;
    std::string tmpname=Form("CDC_%d_%d",l*15,l*15+n);
    gHttp.Register(http::BLDCHitPat(kCDC,tmpname,n,l*15,5,3), "CDC");
    gHttp.Register(http::BLDCMulti(kCDC,tmpname,n,l*15,5,3), "CDC");
  }

#if WIRE
  for(int l=0;l<118;l++){ 
    gHttp.Register(http::BLDCWIRE(kCDC,"CDC",l,16,4,4), "CDC");
  }
#endif

  //=== set directory ===//
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }
  //=== set directory ===//

  std::ifstream data("/home/oper/online_br/CDCmap.txt");
  std::string str;
  if (data.fail()) {
    std::cerr << "Failed to open CDC-map file." << std::endl;
    return -1;
  }
  int layer, wire, asdnum; // 1 origin
  int asdch; // 0 origin
  while(!data.eof()){
    data >> layer >> wire >> asdnum >> asdch;
    CDCmap[1000*(asdnum-1)+asdch] = 1000*(layer-1)+(wire-1);
  }
  
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
  // Hodoscope ------------------------------------------------------------
  int hid;
  {
    // data type    
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
      const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
      int multiplicity=0;
      int mul[2]={0,0};
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	for(int ud=0; ud<nud[i]; ++ud){	  
	  for(int idata=0; idata<ndata; ++idata){
	    int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	    if(data[idata]==k_leading&&nhit>0){
	      hid  = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
	      hptr_array[hid]->Fill(seg);	    
	      mul[ud]++;
	      ntdc[ud]=nhit;	  
	    }
	    for( int m=0; m<nhit; ++m ){
	      int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	      hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	      //std::cout<< hodo_name[i].Data() << " seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	      hptr_array[hid]->Fill(tdc);	    
	      if(data[idata]==k_adc&&gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
		      hid  = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
		      hptr_array[hid]->Fill(tdc);	  		
	      } 
	    } // nhit
	  } // data
	}//ud
	if(ntdc[0]>0&&ntdc[1]>0){
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
    } //hodo
  }

  // Aerogel 
  {
    double ac_sum=0;
    DetectorType kDET=kAC;
    const int k_device = gUnpacker.get_device_id("AC");
    const int nch = 4;
    const int k_ac_leading=4;
    bool ACHIT=false;
    // tdc leading
    int nhit=gUnpacker.get_entries(k_device, 0, 0, 0, k_ac_leading);
    if(nhit>0){
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0, 0, k_ac_leading , m);
	hid  = gHist.getSequentialID(kDET, 0, kTDC, 1);
	if(tdc<6.7e5&&tdc>6.3e5) ACHIT=true;
	hptr_array[hid]->Fill(tdc);	    
      }
    }
    // adc
    for(int ch=0; ch<nch; ++ch){
      // for(int idata=0; idata<ndata; ++idata){
      int nhit = gUnpacker.get_entries(k_device, 0, 0, 0, ch);
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0, 0 ,ch, m);
	hid  = gHist.getSequentialID(kDET, 0, kADC, ch+1);
	//  std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	ac_sum+=tdc;
	hptr_array[hid]->Fill(tdc);	    
	if(ACHIT){	   
	  hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, ch+1);
	  hptr_array[hid]->Fill(tdc);	    
	} 
      } // nhit
    }//ch
    hid = gHist.getSequentialID(kDET, 0, kADC, 5);
    hptr_array[hid]->Fill(ac_sum);	    
    if(ACHIT){
      hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, 5);
      hptr_array[hid]->Fill(ac_sum);	    
    }
  }
  // Chamber ------------------------------------------------------------
  // BLDCs
  for(int i=0;i<nchm-1;i++)
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

	// traling
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

  // CDC
  {
    int i = nchm-1;
    // data type
    DetectorType kDET=chm[i];
    const int k_device = gUnpacker.get_device_id(chm_name[i].Data());
    const int k_l    = 0;
    const int k_t    = 1;
    int mulCDC[15] = {0,0,0,0,0,
		      0,0,0,0,0,
		      0,0,0,0,0};
    DetectorType kDET2 = kCDC2;    
    // TDC & HitPat & Multi
    TString hname;
    int hid=-1;
    for(int l = 0; l<nlayers[i]; ++l){ // ASD number
      int multiplicity    = 0;
      for(int w = 0; w<nwires[i]; ++w){ // ASD channel
	if( !CDCmap.count(1000*l+w) ) continue;
	int val = CDCmap[1000*l+w];
	int layer = val/1000;
	int wire  = val%1000;
	int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	// std::cerr<<l<<" "<<w<<" "<<layer<<" "<<wire<<" "<<nhit<<std::endl;
	if( nhit==0 ) continue;
	// This wire fired at least one times.
	++multiplicity;
	mulCDC[layer]++;
	hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	hptr_array[hid]->Fill(w, nhit);
	hid=gHist.getSequentialID(kDET2,0,kHitPat,layer+1);
	hptr_array[hid]->Fill(wire);
	std::vector< int > leading_array(nhit);
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	  hid=gHist.getSequentialID(kDET,l+1,kTDC,w+1);
	  hptr_array[hid]->Fill(tdc);
	  hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	  hptr_array[hid]->Fill(tdc);
	  hid=gHist.getSequentialID(kDET2,0,kTDC,layer+1);
	  hptr_array[hid]->Fill(tdc);
	  leading_array[m] = tdc;
	}	

	nhit = gUnpacker.get_entries(k_device, l, 0, w, k_t);
	if( nhit==0 ) continue;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_t , m);
	  hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	  hptr_array[hid]->Fill(tdc);
	  hid=gHist.getSequentialID(kDET2,0,kTDC2D,layer+1);
	  hptr_array[hid]->Fill(tdc);
		if(m<leading_size){
			hid=gHist.getSequentialID(kDET,l+1,kTOT,l+1);
			hptr_array[hid]->Fill(tdc - leading_array[m]);
			hid=gHist.getSequentialID(kDET2,0,kTOT,layer+1);
			hptr_array[hid]->Fill(leading_array[m] - tdc);

			hid=gHist.getSequentialID(kDET2,0,kADC2D,layer+1);
			hptr_array[hid]->Fill(leading_array[m],leading_array[m] - tdc);
		}
	} 
      } // for(int w = 0; w<nwires[i]; ++w){
      hid=gHist.getSequentialID(kDET,0,kMulti,l+1);
      hptr_array[hid]->Fill(multiplicity);
    }

    for( int i=0; i<15; i++ ){
      hid=gHist.getSequentialID(kDET2,0,kMulti,i+1);
      hptr_array[hid]->Fill(mulCDC[i]);
    }
  }
  // TriggerFlag
  {
    // data type
    int hid=-1;
    const int nmhtdc=1;
    DetectorType mhtdc[nmhtdc]={kTriggerFlag};
    TString mhtdc_name[nmhtdc]={"TriggerFlag"};
    const int nsegs[nmhtdc]={16};
    const int k_leading=0;
    const int k_trailing=1;
    const int ndata=2;
    const int data[ndata]={k_leading,k_trailing};
    const int type[ndata]={kTDC,kTDC2D};
    const int ud=0;
    for(int i=0;i<nmhtdc;i++){
      DetectorType kDET=mhtdc[i];
      const int k_device = gUnpacker.get_device_id(mhtdc_name[i].Data());
      int multiplicity=0;
      int mul[2]={0,0};
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	for(int idata=0; idata<ndata; ++idata){
	  int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	  if(data[idata]==k_leading&&nhit>0){
	    hid  = gHist.getSequentialID(kDET, 0, kHitPat, 1);
	    hptr_array[hid]->Fill(seg);	    
	    mul[ud]++;
	    ntdc[ud]=nhit;	  
	  }
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	    hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	    hptr_array[hid]->Fill(tdc);	    
	  } // nhit
	} // data
      }// seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(multiplicity);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(mul[1]);	    
    } // hodo
  }
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

