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
#include <TRandom3.h>

#include "BLDCWireMapMan.hh"
#include "XTMapMan.hh"       
#include "HodoAnalyzer.hh"
#include "DCAnalyzer.hh"
#include "Controller.hh"
#include "HttpServer.hh"
#include "Updater.hh"
#include "LocalTrack.hh"

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
    //    HodoParamMan&   gHodo      = HodoParamMan::GetInstance();
    UserParamMan&  gUser= UserParamMan::GetInstance();
    const XTMapMan& gXt  = XTMapMan::GetInstance();
    const BLDCWireMapMan& gWire = BLDCWireMapMan::GetInstance();
    DCAnalyzer   *DCAna;
    HodoAnalyzer *hodoAna;
    TString flagnames[16]={"SpillStart","SpillEnd",
			   "Beam","Pion","Kaon2","Kaon3",
			   "KxCDH1","KxCDH2","KxCDH3","KxCDH1xG","KaonxG",
			   "PixCDH","BxPbF2",
			   "CDH cosmic","none","clock(10s)"};
  }

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{  
  ConfMan& gConfMan = ConfMan::getInstance();
  gConfMan.initialize(argv);
  gConfMan.initializeHodoParamMan();
  gConfMan.initializeUserParamMan();
  gConfMan.initializeXTMapMan();
  gConfMan.initializeBLDCWireMapMan();
  gConfMan.initializeDCTdcCalibMan();
  gUser.Print();
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
  gHttp.SetPort(8080);
  gHttp.Open();

  int tdcnbins=2000;
  double tdcmin=0;
  double tdcmax=2e6;
  gHttp.Register(gHist.createHodo(kBHD,"BHD",16,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT0, "T0",5,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kAC ,"AC",5,1,1024,0,2048,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT0new,"T0new",1,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kDEF,"DEF",5,2,4096,0,4096,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kCDH,"CDH"  ,36, 2, 4096,0,4096,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kVeto,"Veto",2,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kPbF2,"PbF2",40,1,2048,0,2048,tdcnbins,tdcmin,tdcmax));
  //gHttp.Register(gHist.createHodo(kRC,"RC",8,2,2048,0,2048,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createMHTDC(kTriggerFlag ,"Triggerflag",16,4000));
  gHttp.Register(gHist.createBLDC(kBPC,"BPC",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  gHttp.Register(gHist.createAna(1));

  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::Check(1,"QDC"),"Analysis");
  gHttp.Register(http::Ana(1),"Analysis");
  gHttp.Register(http::Ana(2),"Analysis");
  gHttp.Register(http::TOF(0,"BHD"),"Analysis");
  gHttp.Register(http::AC(5,"All"),"Analysis");
  gHttp.Register(http::AC(6,"Beam"),"Analysis");
  gHttp.Register(http::BLDCProf(kBLC1a,"BLC1a"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC1b,"BLC1b"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC2a,"BLC2a"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC2b,"BLC2b"),"Profile2");
  gHttp.Register(http::BLDCProf(kBPC,"BPC"),"Profile2");
  gHttp.Register(http::BLDCXYProf(kBLC1a,"BLC1a",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC1b,"BLC1b",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a_FF",10),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b_FF",10),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC,"BPC",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC,"BPC_FF",10),"Profile");

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

  gHttp.Register(http::QDC(kBHD,"_BHD_U",0,0,16,4,4,0,500),"BHD");
  gHttp.Register(http::QDC(kBHD,"_BHD_D",1,0,16,4,4,0,500),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_U",0,16,4,4),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_D",1,16,4,4),"BHD");
  gHttp.Register(http::MHTDCMeanTime(kBHD,"_BHD",0,16,4,4),"BHD");

  gHttp.Register(http::QDC(kT0,"_T0_U",0,0,5,3,2,0,500),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,0,5,3,2,0,500),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");
  gHttp.Register(http::MHTDCMeanTime(kT0,"_T0",0,5,3,2),"T0");

  gHttp.Register(http::QDC(kT0new,"_T0new_U",0,0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::QDC(kT0new,"_T0new_D",1,0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::MHTDCMeanTime(kT0new,"_T0new",0,1,1,1),"T0new");

  gHttp.Register(http::QDC(kAC,"AC",0,0,5,3,2,0,2000),"AC");

  gHttp.Register(http::QDC(kDEF,"_DEF_U",0,0,5,3,2,0,2000),"DEF");
  gHttp.Register(http::QDC(kDEF,"_DEF_D",1,0,5,3,2,0,2000),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_U",0,5,3,2),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_D",1,5,3,2),"DEF");
  gHttp.Register(http::MHTDCMeanTime(kDEF,"_DEF",0,5,3,2),"DEF");

  gHttp.Register(http::QDC(kCDH,"_CDH_U",0,0,36,6,6,0,2000),"CDH");
  gHttp.Register(http::QDC(kCDH,"_CDH_D",1,0,36,6,6,0,2000),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_U",0,36,6,6),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_D",1,36,6,6),"CDH");
  gHttp.Register(http::MHTDCHitPatMulti(kCDH,"_CDH"),"CDH");

  gHttp.Register(http::QDC(kVeto,"_Veto_L",0,0,2,2,2,0,4000),"Veto");
  gHttp.Register(http::QDC(kVeto,"_Veto_R",1,0,2,2,2,0,4000),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_L",0,2,2,2),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_R",1,2,2,2),"Veto");
  gHttp.Register(http::MHTDCMeanTime(kVeto,"_Veto",0,2,2,2),"Veto");

  gHttp.Register(http::QDC(kPbF2,"_PbF2",0,0,40,8,5,0,1200),"PbF2");
  gHttp.Register(http::MHTDCTDC(kPbF2,"_PbF2",0,40,8,5),"PbF2");

  gHttp.Register(http::BLDCHitPat(kBLC1a,"_BLC1a",8));
  gHttp.Register(http::BLDCMulti(kBLC1a,"_BLC1a",8));
  gHttp.Register(http::BLDCHitPat(kBLC1b,"_BLC1b",8));
  gHttp.Register(http::BLDCMulti(kBLC1b,"_BLC1b",8));
  gHttp.Register(http::BLDCHitPat(kBLC2a,"_BLC2a",8));
  gHttp.Register(http::BLDCMulti(kBLC2a,"_BLC2a",8));
  gHttp.Register(http::BLDCHitPat(kBLC2b,"_BLC2b",8));
  gHttp.Register(http::BLDCMulti(kBLC2b,"_BLC2b",8));
  gHttp.Register(http::BLDCHitPat(kBPC,"_BPC",8));
  gHttp.Register(http::BLDCMulti(kBPC,"_BPC",8));

  gHttp.Register(http::BLDCTDC(kBPC,"BPC",8), "BPC");
  gHttp.Register(http::BLDCTDC(kBLC1a,"BLC1a",8), "BLC1a");
  gHttp.Register(http::BLDCTDC(kBLC1b,"BLC1b",8), "BLC1b");
  gHttp.Register(http::BLDCTDC(kBLC2a,"BLC2a",8), "BLC2a");
  gHttp.Register(http::BLDCTDC(kBLC2b,"BLC2b",8), "BLC2b");
  for(int l=0;l<8;l++){ 
    gHttp.Register(http::BLDCWIRE(kBPC,"BPC",l,32,8,4), "BPC");
    gHttp.Register(http::BLDCWIRE(kBLC1a,"BLC1a",l,32,8,4), "BLC1a");
    gHttp.Register(http::BLDCWIRE(kBLC1b,"BLC1b",l,32,8,4), "BLC1b");
    gHttp.Register(http::BLDCWIRE(kBLC2a,"BLC2a",l,32,8,4), "BLC2a");
    gHttp.Register(http::BLDCWIRE(kBLC2b,"BLC2b",l,32,8,4), "BLC2b");
  }
  for( Int_t i=0, n=hptr_array.size(); i<n; ++i ){
    hptr_array[i]->SetDirectory(0);
  }

  gHttp.Begin();  
  std::cout<<__func__<<std::endl;
  return 0;
}
  
  //____________________________________________________________________________
  int
  process_end( void )
  {
    TFile *root_file;
    root_file=new TFile("tmp.root","recreate");
    root_file->cd();
    for(int i=0;i<hptr_array.size();i++)
      if(hptr_array[i]) hptr_array[i]->Write();
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

  RawData  *rawData=new RawData;
  rawData->DecodeHits();
  DCAna=new DCAnalyzer;
  hodoAna=new HodoAnalyzer;
  hodoAna->DecodeRawHits( rawData);


#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  // MHTDC only
  bool TRIG[16]={0};
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
	TRIG[seg]=false;
	int ntdc[2]={0,0};
	bool beamtiming[2]={false,false};
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
	    if(tdc>gUser.GetParameter("TRIGTDC",0)
	       &&tdc<gUser.GetParameter("TRIGTDC",1)){
	      TRIG[seg]=true;
	    }
	  } // nhit
	} // data	
      }//seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(multiplicity);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(mul[1]);	    
    } //hodo
  }
  // Hodoscope ------------------------------------------------------------
  int hid;
  double time0=-9999;
  bool KAON=false;
  bool PION=false;
  bool ACHIT=false;  
  int t0segment=-1;
  int DEFsegment=-1;
  double DEFADC[2]={-1,-1};
  {
    // data type
    const int nhodo=7;
    int Cid[nhodo]={DetIdT0new,DetIdBHD,DetIdT0,DetIdDEF,DetIdVeto1,DetIdPbF2,DetIdCDH}; 
    DetectorType hodo[nhodo]={kT0new, kBHD, kT0, kDEF, kVeto,kPbF2,kCDH};
    TString hodo_name[nhodo]={"T0new","BHD","T0","DEF","Veto","PbF2","CDH"};
    const int nsegs[nhodo]={1,16,5,5,2,40,36};
    const int nud[nhodo]={2,2,2,2,2,1,2};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};
    // hodo rawdata    
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
      int segud[2]={-1,-1};
      const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
      int multiplicity=0;
      int mul[2]={0,0};
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	bool beamtiming[2]={false,false};
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
	      //		std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	      hptr_array[hid]->Fill(tdc);	    
	      if(data[idata]==k_leading){
		switch(hodo[i]){
		case kBHD:
		  if(tdc>gUser.GetParameter("BHDTDC",0)&&tdc<gUser.GetParameter("BHDTDC",1)){
		    segud[ud]=seg;
		    beamtiming[ud]=true;
		  }
		  break;
		default:
		  if(tdc>gUser.GetParameter("HODOTDC",0)&&tdc<gUser.GetParameter("HODOTDC",1)){
		    segud[ud]=seg;
		    beamtiming[ud]=true;
		  }
		}
	      }
	      if(data[idata]==k_adc&&gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
		if((hodo[i]!=kT0new)&&!beamtiming[0]) continue;
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
	  DEFADC[0] = gUnpacker.get(k_device, 0, seg, 0, k_adc , 0);
	  DEFADC[1] = gUnpacker.get(k_device, 0, seg, 1, k_adc , 0);
	}
      }//seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);      hptr_array[hid]->Fill(multiplicity);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);      hptr_array[hid]->Fill(mul[1]);	    
      if(segud[0]>-1&&segud[1]>-1){
	hid  = gHist.getSequentialID(kDET, 0, kHitPat2D, 0); hptr_array[hid]->Fill(segud[0],segud[1]);	    
      }
    } //hodo

    // AC Hit?
    ACHIT=false;
    {
      DetectorType kDET=kAC;
      const int k_device = gUnpacker.get_device_id("AC");
      const int k_leading = 4; 
      int nhit=gUnpacker.get_entries(k_device, 0, 0, 0, k_leading);
      if(nhit>0){
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, 0,0, k_leading , m);
	  hid  = gHist.getSequentialID(kDET, 0, kTDC,  1);
	  if(tdc>gUser.GetParameter("ACTDC",0)
	     &&tdc<gUser.GetParameter("ACTDC",1)){
	    ACHIT=true;
	  }
	  hptr_array[hid]->Fill(tdc);	    
	  //	  std::cout<<hptr_array[hid]->GetTitle()<<"  "<<tdc<<std::endl;
	}
      }
    }
    // hodo tof
    for(int ihodo=0;ihodo<4;++ihodo){
      DetectorType kDET=hodo[ihodo];
      int nh = hodoAna->GetNHits(Cid[ihodo]);
      //    std::cout<<name[ihodo]<<"  "<<nh<<std::endl;
      //    HF1( Hid[ihodo]+10, double(nh) );
      for( int i=0; i<nh; ++i ){
	Hodo2Hit *hit = hodoAna->GetHit(Cid[ihodo],i);
	HodoRawHit *raw = hit->GetRawHit();
	if(!hit) continue;
	int seg = hit->SegmentId()+1;
	for(int it=0;it<hit->GetIndex();it++){
	  double mt  = hit->MeanTime(it); 
	  hid  = gHist.getSequentialID(kDET, 0, kCTime,  seg);
	  hptr_array[hid]->Fill(mt);
	  if(kDET==kT0new&&mt>-20&&mt<20) time0=mt;
	  double tof=mt-time0;
	  if(kDET==kT0&&mt>-50&&mt<50){
	    t0segment=seg;	    
	  }
	  if(kDET==kDEF&&mt>-50&&mt<50){
	    DEFsegment=seg;	    
	  }
	  if(kDET==kBHD){
	    if(tof>gUser.GetParameter("TOFK",0)
	       &&tof<gUser.GetParameter("TOFK",1)){
	      KAON=true;
	    }
	    if(tof>gUser.GetParameter("TOFPi",0)
	       &&tof<gUser.GetParameter("TOFPi",1)){
	      PION=true;
	    }
	  }	    
	  if(ihodo!=0&&time0>-9000){
	    hid  = gHist.getSequentialID(kAna, 0, 1, 1)+(ihodo-1)*4; 
	    hptr_array[hid]->Fill(tof);
	    if(ACHIT)  hptr_array[hid+1]->Fill(tof);
	    if(KAON)   hptr_array[hid+2]->Fill(tof);
	    if(PION)   hptr_array[hid+3]->Fill(tof);
	  }
	}
      }//for(i)
    }    
  }
  
  // Aerogel 
  {
    double ac_sum=0;
    int nch=4;
    DetectorType kDET=kAC;
    const int k_device = gUnpacker.get_device_id("AC");
    for(int ch=0; ch<nch; ++ch){
      int nhit = gUnpacker.get_entries(k_device, 0, 0, 0, ch);
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0, 0, ch , m);
	hid  = gHist.getSequentialID(kDET, 0, kADC, ch+1);
	ac_sum+=tdc;
	hptr_array[hid]->Fill(tdc);	    
	if(ACHIT){	   
	  hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, ch+1);
	  hptr_array[hid]->Fill(tdc);	    
	} 
      } // nhit
    }//ch
    hid  = gHist.getSequentialID(kDET, 0, kADC, 5);
    hptr_array[hid]->Fill(ac_sum);	    
    if(ACHIT){
      hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, 5);
      hptr_array[hid]->Fill(ac_sum);	    
    }
    hid  = gHist.getSequentialID(kAna, 0, 5, 1);
    hptr_array[hid]->Fill(ac_sum);	    
    if(ACHIT) hptr_array[hid+1]->Fill(ac_sum);	    
    if(KAON&&!PION)  hptr_array[hid+2]->Fill(ac_sum);	        
    if(PION&&!KAON)  hptr_array[hid+3]->Fill(ac_sum);	        
    if(TRIG[2]){
      hid  = gHist.getSequentialID(kAna, 0, 6, 1);
      hptr_array[hid]->Fill(ac_sum);	    
      if(ACHIT) hptr_array[hid+1]->Fill(ac_sum);	    
      if(KAON&&!PION)  hptr_array[hid+2]->Fill(ac_sum);	        
      if(PION&&!KAON)  hptr_array[hid+3]->Fill(ac_sum);	        
    }
  }


  // Chamber ------------------------------------------------------------
  const int nchm=5;
  int ChmCid[nchm]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdBPC};
  DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b, kBPC};
  TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC"};
  const int nwires[nchm]={32,32,32,32,32};
  const int nlayers[nchm]={8,8,8,8,8};
  double rotz[nchm]={45,45,135,135,0};
  double ffz[nchm]={0,0,130,130,20};
  for(int i=0;i<nchm;i++)
    {
      // data type
      DetectorType kDET=chm[i];
      const int k_device = gUnpacker.get_device_id(chm_name[i].Data());
      const int k_l    = 0;
      const int k_t    = 1;      
      // TDC & HitPat & Multi
      TString hname;
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
	    hid=gHist.getSequentialID(kDET,l+1,kTDC,w+1);   hptr_array[hid]->Fill(tdc);
	    hid=gHist.getSequentialID(kDET,0,kTDC,l+1);	    hptr_array[hid]->Fill(tdc);
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
	// if( kDET==kBPC){
	//   std::cout<<"fill: "<<k_device<<"  "<<l<<"  "<<multiplicity<<std::endl;
	//	}
      }
    }
  if(1)
  {    
    DCAna->DecodeRawHits(rawData);
    for(int ichm=0;ichm<5;ichm++){
      DetectorType kDET=chm[ichm];
      int cid=ChmCid[ichm];
      int nl=nlayers[ichm];
      int nw=nwires[ichm];    
      int xy[8]={0,0,1,1,0,0,1,1};
      if(kDET==kBLC2b){
	xy[0]=1;
	xy[1]=1;
	xy[2]=0;
	xy[3]=0;
	xy[4]=1;
	xy[5]=1;
	xy[6]=0;
	xy[7]=0;
      }  
      bool Single=true;
      for( int layer=1; layer<=nl; ++layer ){
	const DCHitContainer &cont = DCAna->GetDCHC(cid, layer);
	int mul=cont.size();
	//	std::cout<<cid<<"  "<<layer<<"  "<<mul<<std::endl;
	if(mul!=1) Single=false;
      }
      LocalTrack *track=new LocalTrack();
      if(Single){
	int wire[8];
	for( int layer=1; layer<=nl; ++layer ){
	  const DCHitContainer &cont = DCAna->GetDCHC(cid, layer);
	  int mul=cont.size();
	  for(int ihit=0;ihit<mul;ihit++){
	    DCHit* hit=cont[ihit];
	    if(mul==1){
	      wire[layer-1]=hit->GetWire();
	      track->AddHit(hit,xy[layer-1]);
	    }
	  }
	}
	for(int i=0;i<4;i++){
	  hid = gHist.getSequentialID(kDET, 0, kWireCorr, i);
	  hptr_array[hid]->Fill(wire[i],wire[i+4]);
	}
	track->DoFit();
	track->ConvLocalToGlobal(rotz[ichm]);
	double x,y;
	track->XYPosatZ(0,x,y);
	double dxdz=track->dx();
	double dydz=track->dy();
	
	hid=gHist.getSequentialID(kDET,0,kProf,1);	hptr_array[hid]->Fill(x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,2);	hptr_array[hid]->Fill(x);
	hid=gHist.getSequentialID(kDET,0,kProf,3);	hptr_array[hid]->Fill(y);
	hid=gHist.getSequentialID(kDET,0,kProf,31);	hptr_array[hid]->Fill(dxdz,dydz);
	hid=gHist.getSequentialID(kDET,0,kProf,32);	hptr_array[hid]->Fill(x,dxdz);
	hid=gHist.getSequentialID(kDET,0,kProf,33);	hptr_array[hid]->Fill(y,dydz);
	if(kDET==kBLC2a){
	  if(t0segment>0){
	    hid=gHist.getSequentialID(kAna,0,2,t0segment);	    hptr_array[hid]->Fill(x,y);
	  }
	}else if(kDET==kBLC2b){
	  if(t0segment>0){
	    hid=gHist.getSequentialID(kAna,0,3,t0segment);	    hptr_array[hid]->Fill(x,y);
	  }
	}else if(kDET==kBPC){
	  if(DEFsegment>0){
	    hid=gHist.getSequentialID(kAna,0,4,DEFsegment);	    hptr_array[hid]->Fill(x,y);
	  }
	}
	if(KAON){
	  hid=gHist.getSequentialID(kDET,0,kProf,4);	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,5);	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,6);	  hptr_array[hid]->Fill(y);
	  if(DEFsegment>0){
	    hid=gHist.getSequentialID(kDET,0,kProf,21);	    hptr_array[hid]->Fill(x,y);
	    hid=gHist.getSequentialID(kDET,0,kProf,22);	    hptr_array[hid]->Fill(x);
	    hid=gHist.getSequentialID(kDET,0,kProf,23);	    hptr_array[hid]->Fill(y);
	  }
	}else if(PION){
	  hid=gHist.getSequentialID(kDET,0,kProf,7);	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,8);	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,9);	  hptr_array[hid]->Fill(y);
	}
	track->XYPosatZ(ffz[ichm],x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,11);	hptr_array[hid]->Fill(x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,12);	hptr_array[hid]->Fill(x);
	hid=gHist.getSequentialID(kDET,0,kProf,13);	hptr_array[hid]->Fill(y);
	if(KAON){
	  hid=gHist.getSequentialID(kDET,0,kProf,14);	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,15);	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,16);	  hptr_array[hid]->Fill(y);
	}else if(PION){
	  hid=gHist.getSequentialID(kDET,0,kProf,17);	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,18);	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,19);	  hptr_array[hid]->Fill(y);
	}
      }
      delete track;
    }
  }

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  if(hodoAna) delete hodoAna;  
  if(DCAna) delete DCAna;  
  if(rawData) delete rawData;
    return 0;
}
  
}

