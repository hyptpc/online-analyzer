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
    TString flagnames[32]={"SpillStart","SpillEnd", //1
			   "Beam","Pion","Kaon2","Kaon3", //5
			   "KxCDH1","KxCDH2","KxCDH3","KxCDH1xG","KaonxG", //10
			   "PixCDH","PixPbF2", // 12
			   "ExPbF2","CDH cosmic","Clock", //15
                           "dummy16","dummy17","dummy18","dummy19",
                           "dummy20","dummy21","dummy22","dummy23",
                           "dummy24","dummy25","dummy26","dummy27",
                           "dummy28","dummy29","dummy30","dummy31"};
    TString outputname="test.root";
    
    const int nchm=6;
    int ChmCid[nchm]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdBPC1,DetIdBPC2};
    DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b, kBPC1,kBPC2};
    TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC1","BPC2"};
    const int nwires[nchm]={32,32,32,32,15,32};
    const int nlayers[nchm]={8,8,8,8,8,8};
    double rotz[nchm]={45,45,135,135,-45,0};
    double ffz[nchm]={0,0,130,130,50,37};
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
  int port=8080;
  if(argv.size()==4){
    outputname=argv.at(3);
    port=8088;
  }
  gHttp.SetPort(port);
  gHttp.Open();

  int tdcnbins=2000;
  double tdcmin=0.9e6;
  double tdcmax=1.4e6;
  gHttp.Register(gHist.createBHT(kBHT,"BHT",63,2,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT0,"T0",5,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kT1,"T1",1,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kDEF,"DEF",5,2,4096,0,4096,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kBTC,"BTC",2,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createBLDC(kBPC1 ,"BPC1" ,8,15,false,true));
  gHttp.Register(gHist.createBLDC(kBPC2 ,"BPC2" ,8,32,false,true));
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,false,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,false,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,false,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,false,true));
  gHttp.Register(gHist.createHodo(kAC ,"AC",5,1,1024,0,2048,2000,0,2000));
  gHttp.Register(gHist.createMHTDC(kTriggerFlag ,"Triggerflag",32,2000));
  gHttp.Register(gHist.createAnaTOF(1));
  gHttp.Register(gHist.createAnaTrack(1));

  std::cout<<"histogram created"<<std::endl;
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::Check(1,"QDC"),"Analysis");
  gHttp.Register(http::Ana(1),"Analysis");
  gHttp.Register(http::Ana(2),"Analysis");
  gHttp.Register(http::Ana(3),"Analysis");
  gHttp.Register(http::Ana(4),"Analysis");
  gHttp.Register(http::Ana(10),"Analysis");
  gHttp.Register(http::TOF(0,"BHT"),"Analysis");
  gHttp.Register(http::TOF(4,"BTC"),"Analysis");
  gHttp.Register(http::TOF2(0,"BHT"),"Analysis");
  gHttp.Register(http::TOF2(4,"BTC"),"Analysis");
  gHttp.Register(http::TOF_Btrg(0,"BHT"),"Analysis");
  gHttp.Register(http::TOF_Btrg(4,"BTC"),"Analysis");
  // gHttp.Register(http::TOF2D(7,"VETO"),"Analysis");
  // gHttp.Register(http::TOF2D(8,"BTC"),"Analysis");
  // gHttp.Register(http::TOF2D(9,"VetodE"),"Analysis");
  // gHttp.Register(http::TOF2D(10,"BTCdE"),"Analysis");
  gHttp.Register(http::AC(10,"All"),"Analysis");
  gHttp.Register(http::AC(11,"Beam"),"Analysis");
  gHttp.Register(http::BLDCProf(kBLC1a,"BLC1a"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC1b,"BLC1b"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC2a,"BLC2a"),"Profile2");
  gHttp.Register(http::BLDCProf(kBLC2b,"BLC2b"),"Profile2");
  gHttp.Register(http::BLDCProf(kBPC1,"BPC1"),"Profile2");
  gHttp.Register(http::BLDCProf(kBPC2,"BPC2"),"Profile2");
  gHttp.Register(http::BLDCXYProf(kBLC1a,"BLC1a",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC1b,"BLC1b",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a_FF",10),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b_FF",10),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC1,"BPC1",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC1,"BPC1_FF",10),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC2,"BPC2",0),"Profile");
  gHttp.Register(http::BLDCXYProf(kBPC2,"BPC2_FF",10),"Profile");
  gHttp.Register(http::BeamAxis(),"Profile");
  gHttp.Register(http::BLDCCorr(0,"BLC1"),"BLDCCorr");
  gHttp.Register(http::BLDCCorr(1,"BLC2"),"BLDCCorr");
  gHttp.Register(http::BLDCCorr(2,"BPC"),"BLDCCorr");

  for(int i=0;i<nchm;i++){
    const std::string tmpstr=Form("_%s",chm_name[i].Data());
    gHttp.Register(http::BLDCResid(chm[i],tmpstr,8),"Residual"); 
    gHttp.Register(http::BLDCMulti(chm[i],tmpstr,8),"Multiplicity");
  }
  gHttp.Register(http::MHTDCTDC(kTriggerFlag,"_TriggerFlag",0,32,8,4),"TriggerFlag");
  gHttp.Register(http::MHTDCHitPatMulti(kTriggerFlag,"_TriggerFlag"),"TriggerFlag");

  std::cout<<"canvas created"<<std::endl;
  {
    int hid1 = gHist.getSequentialID(kTriggerFlag, 0, kHitPat, 1);  
    hptr_array[hid1]->GetXaxis()->SetTitle("");
    for( Int_t i=0; i<32; ++i ){
      int hid2 = gHist.getSequentialID(kTriggerFlag, 0, kTDC, i+1);
      hptr_array[hid2]->SetTitle(Form("%s_%s",hptr_array[hid2]->GetTitle(),flagnames[i].Data()));
      hptr_array[hid1]->GetXaxis()->SetBinLabel(i+1,flagnames[i]);
    }
  }
  std::cout<<"set attributes for some histograms"<<std::endl;
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
    root_file=new TFile(outputname,"recreate");
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
  bool TRIG[32]={0};
  {
    // data type
    int hid=-1;
    const int nmhtdc=1;
    DetectorType mhtdc[nmhtdc]={kTriggerFlag};
    TString mhtdc_name[nmhtdc]={"TriggerFlag"};
    const int nsegs[nmhtdc]={32};
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
	for(int idata=0; idata<ndata; ++idata){
	  int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	  if(data[idata]==k_leading&&nhit>0){
	    hid  = gHist.getSequentialID(kDET, 0, kHitPat, 1);
	    hptr_array[hid]->Fill(seg);	    
	    mul[ud]++;
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
  if(TRIG[14]&&!TRIG[15]) return 0;
  // Hodoscope ------------------------------------------------------------
  int hid;
  double time0=-9999;
  bool KAON=false;
  bool PION=false;
  bool PROTON=false;
  bool DEUTERON=false;
  bool ACHIT=false;  
  bool TOF=false;
  int t0segment=-1;
  int DEFsegment=-1;

  {
    // data type
    const int nhodo=5;
    int Cid[nhodo]={DetIdT1,DetIdBHT,DetIdT0,DetIdDEF,DetIdBTC}; 
    DetectorType hodo[nhodo]={kT1, kBHT, kT0, kDEF,kBTC};
    for(int ihodo=0;ihodo<nhodo;ihodo++){
      DetectorType kDET=hodo[ihodo];
      int detid=Cid[ihodo];
      int mul=0,mulu=0,muld=0;
      if(detid==DetIdBHT) continue;
      const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
      int nh = cont.size();
      int tdc_ll=gUser.GetParameter("HODOTDC",0);
      int tdc_ul=gUser.GetParameter("HODOTDC",1);
      for( int i=0; i<nh; ++i ){
        HodoRawHit *raw = cont[i];
        if(!raw) continue;
	int seg = raw->SegmentId();
        TString segstr=Form("_seg%d",seg);
        double au  = 0;
        double ad  = 0;
        if(detid!=DetIdBHT){
          au  = raw->GetAdcUp();
          ad  = raw->GetAdcDown();
          hid  = gHist.getSequentialID(kDET, 0, kADC, seg+1);
          hptr_array[hid]->Fill(au);	    
          hid  = gHist.getSequentialID(kDET, 1, kADC, seg+1);
          hptr_array[hid]->Fill(ad);	    
        }
        bool ngateu=0;
        bool ngated=0;
        int ntu=raw->GetSizeTdcUp();
        for(int it=0;it<ntu;it++){
          double tu  = raw->GetTdcUp(it);
          hid  = gHist.getSequentialID(kDET, 0, kTDC, seg+1);     
          hptr_array[hid]->Fill(tu);	    
          if(tu>tdc_ll&&tu<tdc_ul){
            ngateu++;
          }
        }
        if(ngateu>0){
          mulu++;
	  hid  = gHist.getSequentialID(kDET, 0, kHitPat, 1);
          hptr_array[hid]->Fill(seg);	  		
          if(detid!=DetIdBHT){            
            hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, seg+1);
            hptr_array[hid]->Fill(au);	  		
          }
        }
        int ntd=raw->GetSizeTdcDown();
        for(int it=0;it<ntd;it++){
          double td  = raw->GetTdcDown(it);
          hid  = gHist.getSequentialID(kDET, 1, kTDC, seg+1);     
          hptr_array[hid]->Fill(td);	    
        }
        if(ngated>0){
          muld++;
	  hid  = gHist.getSequentialID(kDET, 0, kHitPat, 2);
          hptr_array[hid]->Fill(seg);	  
          if(detid!=DetIdBHT){            		
            hid  = gHist.getSequentialID(kDET, 1, kADCwTDC, seg+1);
            hptr_array[hid]->Fill(ad);	  
          }		
        }
        if(ngateu>0 && ngated>0){
          mul++;
          hid  = gHist.getSequentialID(kDET, 0, kHitPat, 0);
          hptr_array[hid]->Fill(seg);	  		
        }
      }
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(mul);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mulu);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(muld);	    
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

    std::vector<double> tof_bht;
    std::vector<double> tof_btc;
    // hodo tof

    //    
    for(int ihodo=0;ihodo<5;++ihodo){
      DetectorType kDET=hodo[ihodo];
      int nh = hodoAna->GetNHits(Cid[ihodo]);
      //    HF1( Hid[ihodo]+10, double(nh) );
      for( int i=0; i<nh; ++i ){
	Hodo2Hit *hit = hodoAna->GetHit(Cid[ihodo],i);
	if(!hit) continue;
	int seg = hit->SegmentId()+1;
	for(int it=0;it<hit->GetIndex();it++){
	  double mt  = hit->MeanTime(it); 
	  if(kDET==kBHT){
	    if(hit->GetAUp(it)  <8000) continue;
	    if(hit->GetADown(it)<8000) continue;
	  }
	  hid  = gHist.getSequentialID(kDET, 0, kCTime,  seg);
	  hptr_array[hid]->Fill(mt);
	  if(kDET==kT1&&mt>-20&&mt<20) time0=mt;
	  if(kDET==kT0&&mt>-50&&mt<50){
	    t0segment=seg;	    
	  }
	  if(kDET==kDEF&&mt>-50&&mt<50){
	    DEFsegment=seg;	    
	  }
	  if(ihodo!=0&&time0>-9000){
	    double tof=mt-time0;
	    if(kDET==kBHD){
	      if(tof>gUser.GetParameter("TOFK",0)
		 &&tof<gUser.GetParameter("TOFK",1)){
		KAON=true;
	      }
	      if(tof>gUser.GetParameter("TOFPi",0)
		 &&tof<gUser.GetParameter("TOFPi",1)){
		PION=true;
	      }
	      if(tof>gUser.GetParameter("TOFP",0)
		 &&tof<gUser.GetParameter("TOFP",1)){
		PROTON=true;
	      }
	      if(tof>gUser.GetParameter("TOFD",0)
		 &&tof<gUser.GetParameter("TOFD",1)){
		DEUTERON=true;
	      }
	    }	    	   
	    if(kDET==kBHT)  tof_bht.push_back(tof);
	    if(kDET==kBTC)  tof_btc.push_back(tof);
	    hid  = gHist.getSequentialID(kAna, 0, 1, 1)+(ihodo-1)*13; 
	    hptr_array[hid]->Fill(tof);
	    if(ACHIT)  hptr_array[hid+1]->Fill(tof);
	    if(KAON)   hptr_array[hid+2]->Fill(tof);
	    if(PION)   hptr_array[hid+3]->Fill(tof);
	    if(PROTON)   hptr_array[hid+4]->Fill(tof);
	    if(DEUTERON)   hptr_array[hid+5]->Fill(tof);
	    if(TRIG[2]){//btrg
	      hptr_array[hid+6]->Fill(tof); //all
	      if(TRIG[23]) hptr_array[hid+10]->Fill(tof); //pion
	      if(TRIG[20]) hptr_array[hid+11]->Fill(tof); //kaon
	      if(TRIG[24]) hptr_array[hid+12]->Fill(tof); //proton
	    }
	    if(TRIG[3])          hptr_array[hid+7]->Fill(tof); //pitrg
	    if(TRIG[4]||TRIG[5]) hptr_array[hid+8]->Fill(tof); //ktrg
	    if(TRIG[15])         hptr_array[hid+9]->Fill(tof); //ptrg
	    if(tof>-10&&tof<10) TOF=true;
	  }
	}
      }//for(i)
    } //hodo  
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
    hid  = gHist.getSequentialID(kAna, 0, 10, 1);
    hptr_array[hid]->Fill(ac_sum);	    
    if(ACHIT) hptr_array[hid+1]->Fill(ac_sum);	    
    if(KAON&&!PION)  hptr_array[hid+2]->Fill(ac_sum);	        
    if(PION&&!KAON)  hptr_array[hid+3]->Fill(ac_sum);	        
    if(TRIG[2]){ // beam trigger
      hid  = gHist.getSequentialID(kAna, 0, 11, 1);
      hptr_array[hid]->Fill(ac_sum);	    
      if(ACHIT) hptr_array[hid+1]->Fill(ac_sum);	    
      if(KAON&&!PION)  hptr_array[hid+2]->Fill(ac_sum);	        
      if(PION&&!KAON)  hptr_array[hid+3]->Fill(ac_sum);	        
    }
  }

  // Chamber ------------------------------------------------------------
  for(int i=0;i<nchm-1;i++)
    {
      // data type
      DetectorType kDET=chm[i];
      const int k_device = gUnpacker.get_device_id(chm_name[i].Data());
      const int k_l    = 0;
      const int k_t    = 1;      
      // TDC & HitPat & Multi
      int t_ll=1300;
      int t_ul=1450;
      int detid=ChmCid[i];
      if(detid==DetIdBLC1a||detid==DetIdBLC1b) t_ll=1200;
      if(detid==DetIdBPC1){ t_ll=1200; t_ul=1400;}
      if(detid==DetIdBPC2){ t_ll=1250; t_ul=1400;}
      TString hname;
      for(int l = 0; l<nlayers[i]; ++l){
	int multiplicity    = 0;
	int mul_gate    = 0;
	for(int w = 0; w<nwires[i]; ++w){
	  int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	  if( nhit==0 ) continue;
	  // This wire fired at least one times.
	  ++multiplicity;
	  hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	  hptr_array[hid]->Fill(w, nhit);
          bool tmp=false;
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	    hid=gHist.getSequentialID(kDET,l+1,kTDC,w+1);   hptr_array[hid]->Fill(tdc);
	    hid=gHist.getSequentialID(kDET,0,kTDC,l+1);	    hptr_array[hid]->Fill(tdc);
            if(tdc>t_ll && tdc<t_ul) tmp=true;
	  }	
          if(tmp) mul_gate++;
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
        if(TRIG[7]){//KxCDH
          hid=gHist.getSequentialID(kDET,0,kMulti2D,l+1);
          hptr_array[hid]->Fill(mul_gate);
        }
	// if( kDET==kBPC){
	//   std::cout<<"fill: "<<k_device<<"  "<<l<<"  "<<multiplicity<<std::endl;
	//	}
      }
    }

  // if(TOF)
    {    
      double xpos[6];
      double ypos[6];
      double xdir[6];
      double ydir[6];
      DCAna->DecodeRawHits(rawData);
      for(int ichm=0;ichm<6;ichm++){
        xpos[ichm]=-999999;
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
	  int tmpmul=cont.size();
	  //	std::cout<<cid<<"  "<<layer<<"  "<<mul<<std::endl;
          int mul=0;
          for(int ihit=0;ihit<tmpmul;ihit++){
            DCHit* hit=cont[ihit];
            for(int itdc=0; itdc<hit->GetTdcSize();itdc++){
              if(hit->IsWithinRange(itdc)){
                mul++;
              }
            }
          }
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
              for(int itdc=0; itdc<hit->GetTdcSize();itdc++){
                if(hit->IsWithinRange(itdc)){
                  wire[layer-1]=hit->GetWire();                
                  track->AddHit(hit,xy[layer-1]);
                  break;
                }
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
	  for(int xy=0;xy<2;xy++){
	    for(int i=0;i<track->nhit(xy);i++){
	      double resid=track->resid(xy,i);
	      int layer=track->layer(xy,i);
	      hid=gHist.getSequentialID(kDET,0,kResid,layer+1);
	      hptr_array[hid]->Fill(resid);
	    }
	  }

	  double dxdz=track->dx();
	  double dydz=track->dy();
          xpos[ichm]=x;	  
          ypos[ichm]=y;	  
          xdir[ichm]=dxdz;	  
          ydir[ichm]=dydz;	  

	  hid=gHist.getSequentialID(kDET,0,kProf,1);	hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,2);	hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,3);	hptr_array[hid]->Fill(y);
	  hid=gHist.getSequentialID(kDET,0,kProf,10);	hptr_array[hid]->Fill(dxdz);
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
	  }else if(kDET==kBPC1){
	    if(DEFsegment>0){
              double tmpx,tmpy;
              track->XYPosatZ(180.,tmpx,tmpy);
	      hid=gHist.getSequentialID(kAna,0,4,DEFsegment);	    hptr_array[hid]->Fill(tmpx,tmpy);
	    }
	  }else if(kDET==kBPC2){
	    if(DEFsegment>0){
              double tmpx,tmpy;
              track->XYPosatZ(50.,tmpx,tmpy);
	      hid=gHist.getSequentialID(kAna,0,5,DEFsegment);	    hptr_array[hid]->Fill(tmpx,tmpy);
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
	  track->XYPosatZ(ffz[ichm]*10,x,y);
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
	  if(kDET==kBPC2){
	    track->XYPosatZ(370,x,y);
	    hid=gHist.getSequentialID(kDET,0,kProf,41);	hptr_array[hid]->Fill(x,y);
	    hid=gHist.getSequentialID(kDET,0,kProf,42);	hptr_array[hid]->Fill(x);
	    hid=gHist.getSequentialID(kDET,0,kProf,43);	hptr_array[hid]->Fill(y);
	    if(KAON){
	      hid=gHist.getSequentialID(kDET,0,kProf,44);	  hptr_array[hid]->Fill(x,y);
	      hid=gHist.getSequentialID(kDET,0,kProf,45);	  hptr_array[hid]->Fill(x);
	      hid=gHist.getSequentialID(kDET,0,kProf,46);	  hptr_array[hid]->Fill(y);
	    }else if(PION){
	      hid=gHist.getSequentialID(kDET,0,kProf,47);	  hptr_array[hid]->Fill(x,y);
	      hid=gHist.getSequentialID(kDET,0,kProf,48);	  hptr_array[hid]->Fill(x);
	      hid=gHist.getSequentialID(kDET,0,kProf,49);	  hptr_array[hid]->Fill(y);
	    }
	  }
	}
	delete track;
      }
      if(xpos[0]>-9999&&xpos[1]>-9999){
        hid=gHist.getSequentialID(kAna,0,40,1);
        hptr_array[hid+0]->Fill(xpos[0],xpos[1]); 
        hptr_array[hid+1]->Fill(ypos[0],ypos[1]); 
        hptr_array[hid+2]->Fill(xdir[0],xdir[1]); 
        hptr_array[hid+3]->Fill(ydir[0],ydir[1]); 
      }
      if(xpos[2]>-9999&&xpos[3]>-9999){
        hid=gHist.getSequentialID(kAna,0,40,11);
        hptr_array[hid+0]->Fill(xpos[2],xpos[3]); 
        hptr_array[hid+1]->Fill(ypos[2],ypos[3]); 
        hptr_array[hid+2]->Fill(xdir[2],xdir[3]); 
        hptr_array[hid+3]->Fill(ydir[2],ydir[3]); 
      }
      if(xpos[4]>-9999&&xpos[5]>-9999){
        hid=gHist.getSequentialID(kAna,0,40,21);
        hptr_array[hid]  ->Fill(xpos[4],xpos[5]); 
        hptr_array[hid+1]->Fill(ypos[4],ypos[5]); 
        hptr_array[hid+2]->Fill(xdir[4],xdir[5]); 
        hptr_array[hid+3]->Fill(ydir[4],ydir[5]); 
      }  
    }

    
    if( gUnpacker.get_counter()%100 == 0 ){
      auto prev_level = gErrorIgnoreLevel;
      gErrorIgnoreLevel = kError;
      http::UpdateBLDCEfficiency();
      gErrorIgnoreLevel = prev_level;
    }
    
#if DEBUG
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
    if(hodoAna) delete hodoAna;  
    if(DCAna) delete DCAna;  
    if(rawData) delete rawData;
    gSystem->ProcessEvents();     
    return 0;
}
  
}

