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
    const XTMapMan& gXt = XTMapMan::GetInstance();
    const BLDCWireMapMan& gWire = BLDCWireMapMan::GetInstance();
    DCAnalyzer   *DCAna;
    HodoAnalyzer *hodoAna;
    TText text;
    TText end;
    TFile *root_file;
    // double sdd_ch2e[32];
    // double sdd_offset[32];
  }

//____________________________________________________________________________
int
process_begin( const std::vector<std::string>& argv )
{
  root_file=new TFile("tmp.root","recreate");
  // std::ifstream ifs("/home/oper/sddcalib.txt");
  // for(int i=0;i<32;i++){
  //   sdd_ch2e[i]=sdd_offset[i]=0;
  // }
  // if(ifs.good()){
  //   int ch;
  //   double ch2e, offs;
  //   while(ifs>>ch>>offs>>ch2e){
  //     sdd_ch2e[ch]=ch2e;
  //     sdd_offset[ch]=offs;
  //   }
  // }
  // for(int i=0;i<32;i++){
  //   std::cout<<i<<"  "<<sdd_ch2e[i]<<"  "<<sdd_offset[i]<<std::endl;
  // }
  
  ConfMan& gConfMan = ConfMan::getInstance();
  gConfMan.initialize(argv);
  gConfMan.initializeHodoParamMan();
  gConfMan.initializeXTMapMan();
  gConfMan.initializeBLDCWireMapMan();
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

  int tdcnbins=3000;
  double tdcmin=0;
  double tdcmax=5e6;
  gHttp.Register(gHist.createAna(1));
  //  gHttp.Register(gHist.createSDD(kSDD1,"SDD"));
  gHttp.Register(gHist.createHodo(kBHD,"BHD",16,2,1024,0,1024,4096,0,4096));
  gHttp.Register(gHist.createHodo(kT0, "T0",5,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kAC ,"AC",5,1,1024,0,2048,1024,0,4096));
  gHttp.Register(gHist.createHodo(kT0new,"T0new",1,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kE0 ,"E0",3,2,2048,0,2048,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kDEF,"DEF",3,2,4096,0,4096,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kSTART ,"START",8,2,1024,0,1024,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kSTOP ,"STOP",5,2,1024,0,2048,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createHodo(kSDD1 ,"SDD",32,1,1024,0,4096,1024,0,10000));
  gHttp.Register(gHist.createHodo(kSDDGate ,"SDDGate",4,1,1024,0,4096,1024,0,8192));
  gHttp.Register(gHist.createHodo(kSDDReset ,"SDDReset",4,1,1024,0,4096,1024,0,8192));
  gHttp.Register(gHist.createHodo(kSDDVeto ,"SDDVeto",4,1,1024,0, 2048,tdcnbins,tdcmin,tdcmax));
  gHttp.Register(gHist.createMHTDC(kTriggerFlag ,"Triggerflag",16,2000));
  gHttp.Register(gHist.createBLDC(kBPC,"BPC",8,16,true));
  gHttp.Register(gHist.createBLDC(kSDC,"SDC",8,16,true));
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  gHttp.Register(gHist.createBLDC(kFDC,"FDC",6,64,true));
  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  gHttp.Register(http::BLDCXYProf(kBLC1a,"BLC1a",0),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC1b,"BLC1b",0),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a",0),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b",0),"Profile");;
  gHttp.Register(http::BLDCXYProf(kSDC,"SDC",0),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC2a,"BLC2a_FF",10),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b_FF",10),"Profile");;
  gHttp.Register(http::BLDCXYProf(kSDC,"SDC_FF",10),"Profile");;
  gHttp.Register(http::BLDCXYProf(kBLC2b,"BLC2b_DEF",20),"Profile");;

  gHttp.Register(http::MHTDCTDC(kTriggerFlag,"_TriggerFlag",0,16,4,4),"TriggerFlag");
  // gHttp.Register(http::SDD(kSDD1,"_Energy1",0,kEnergy1,32,8,4,-2,18),"SDD");
  // gHttp.Register(http::SDD(kSDD1,"_Energy2",0,kEnergy2,32,8,4,-2,18),"SDD");
  // gHttp.Register(http::SDD(kSDD1,"_Energy3",0,kEnergy3,32,8,4,-2,18),"SDD");
  // gHttp.Register(http::SDD(kSDD1,"_Energy4",0,kEnergy4,32,8,4,-2,18),"SDD");
  // gHttp.Register(http::SDD(kSDD1,"_Energy5",0,kEnergy5,32,8,4,-2,18),"SDD");
  gHttp.Register(http::QDC(kSDD1,"_SDD",0,32,8,4,0,4000),"SDD");
  gHttp.Register(http::MHTDCTDC(kSDD1,"_SDD",0,32,8,4),"SDD");
  gHttp.Register(http::MHTDCTDC(kSDDGate,"_SDDgate",0,4,2,2),"SDD");
  gHttp.Register(http::MHTDCTDC(kSDDReset,"_SDDreset",0,4,2,2),"SDD");

  gHttp.Register(http::QDC(kBHD,"_BHD_U",0,16,4,4,0,1000),"BHD");
  gHttp.Register(http::QDC(kBHD,"_BHD_D",1,16,4,4,0,1000),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_U",0,16,4,4),"BHD");
  gHttp.Register(http::MHTDCTDC(kBHD,"_BHD_D",1,16,4,4),"BHD");
  gHttp.Register(http::MHTDCMeanTime(kBHD,"_BHD",0,16,4,4),"BHD");

  gHttp.Register(http::QDC(kT0,"_T0_U",0,5,3,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,5,3,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");
  gHttp.Register(http::MHTDCMeanTime(kT0,"_T0",0,5,3,2),"T0");

  gHttp.Register(http::QDC(kT0new,"_T0new_U",0,1,1,1,0,1000),"T0new");
  gHttp.Register(http::QDC(kT0new,"_T0new_D",1,1,1,1,0,1000),"T0new");
  gHttp.Register(http::MHTDCMeanTime(kT0new,"_T0new",0,1,1,1),"T0new");

  gHttp.Register(http::QDC(kAC,"AC",0,5,3,2,0,2000),"AC");

  gHttp.Register(http::QDC(kE0,"_E0_U",0,3,2,2,0,1000),"E0");
  gHttp.Register(http::QDC(kE0,"_E0_D",1,3,2,2,0,1000),"E0");
  gHttp.Register(http::MHTDCTDC(kE0,"_E0_U",0,3,2,2),"E0");
  gHttp.Register(http::MHTDCTDC(kE0,"_E0_D",1,3,2,2),"E0");

  gHttp.Register(http::QDC(kDEF,"_DEF_U",0,3,2,2,0,4000),"DEF");
  gHttp.Register(http::QDC(kDEF,"_DEF_D",1,3,2,2,0,4000),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_U",0,3,2,2),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_D",1,3,2,2),"DEF");

  gHttp.Register(http::QDC(kSTART,"_START_U",0,8,4,2,0,1000),"START");
  gHttp.Register(http::QDC(kSTART,"_START_D",1,8,4,2,0,1000),"START");
  gHttp.Register(http::MHTDCTDC(kSTART,"_START_U",0,8,4,2),"START");
  gHttp.Register(http::MHTDCTDC(kSTART,"_START_D",1,8,4,2),"START");

  gHttp.Register(http::QDC(kSTOP,"_STOP_U",0,5,3,2,0,2000),"STOP");
  gHttp.Register(http::QDC(kSTOP,"_STOP_D",1,5,3,2,0,2000),"STOP");
  gHttp.Register(http::MHTDCTDC(kSTOP,"_STOP_U",0,5,3,2),"STOP");
  gHttp.Register(http::MHTDCTDC(kSTOP,"_STOP_D",1,5,3,2),"STOP");

  gHttp.Register(http::QDC(kSDDVeto,"_SDDVeto",0,4,2,2,0,2000),"SDD");
  gHttp.Register(http::MHTDCTDC(kSDDVeto,"_SDDVeto",0,4,2,2),"SDD");
  
  gHttp.Register(http::BLDCHitPat(kBLC1a,"_BLC1a",8));
  gHttp.Register(http::BLDCMulti(kBLC1a,"_BLC1a",8));
  gHttp.Register(http::BLDCHitPat(kBLC1b,"_BLC1b",8));
  gHttp.Register(http::BLDCMulti(kBLC1b,"_BLC1b",8));
  gHttp.Register(http::BLDCHitPat(kBLC2a,"_BLC2a",8));
  gHttp.Register(http::BLDCMulti(kBLC2a,"_BLC2a",8));
  gHttp.Register(http::BLDCHitPat(kBLC2b,"_BLC2b",8));
  gHttp.Register(http::BLDCMulti(kBLC2b,"_BLC2b",8));
  gHttp.Register(http::BLDCHitPat(kFDC,"_FDC",6));
  gHttp.Register(http::BLDCMulti(kFDC,"_FDC",6));
  gHttp.Register(http::BLDCHitPat(kSDC,"_SDC",8));
  gHttp.Register(http::BLDCMulti(kSDC,"_SDC",8));
  gHttp.Register(http::BLDCHitPat(kBPC,"_BPC",8));
  gHttp.Register(http::BLDCMulti(kBPC,"_BPC",8));

  gHttp.Register(http::BLDCTDC(kBPC,"BPC",8), "BPC");
  gHttp.Register(http::BLDCTDC(kSDC,"SDC",8), "SDC");
  gHttp.Register(http::BLDCTDC(kBLC1a,"BLC1a",8), "BLC1a");
  gHttp.Register(http::BLDCTDC(kBLC1b,"BLC1b",8), "BLC1b");
  gHttp.Register(http::BLDCTDC(kBLC2a,"BLC2a",8), "BLC2a");
  gHttp.Register(http::BLDCTDC(kBLC2b,"BLC2b",8), "BLC2b");
  gHttp.Register(http::BLDCTDC(kFDC,"FDC",6), "FDC");
  for(int l=0;l<8;l++){ 
    gHttp.Register(http::BLDCWIRE(kBPC,"BPC",l,16,4,4), "BPC");
    gHttp.Register(http::BLDCWIRE(kSDC,"SDC",l,16,4,4), "SDC");
    gHttp.Register(http::BLDCWIRE(kBLC1a,"BLC1a",l,32,8,4), "BLC1a");
    gHttp.Register(http::BLDCWIRE(kBLC1b,"BLC1b",l,32,8,4), "BLC1b");
    gHttp.Register(http::BLDCWIRE(kBLC2a,"BLC2a",l,32,8,4), "BLC2a");
    gHttp.Register(http::BLDCWIRE(kBLC2b,"BLC2b",l,32,8,4), "BLC2b");
  }
  for(int l=0;l<6;l++){ 
    gHttp.Register(http::BLDCWIRE(kFDC,"FDC",l,64,8,8), "FDC");
  }
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
  DCAna->DecodeRawHits( rawData);

#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
  // MHTDC only
  {
    // data type
    int hid=-1;
    const int nmhtdc=3;
    DetectorType mhtdc[nmhtdc]={kSDDGate, kSDDReset,kTriggerFlag};
    TString mhtdc_name[nmhtdc]={"SDDGate","SDDReset","TriggerFlag"};
    const int nsegs[nmhtdc]={4,4,16};
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
  double time0=0;
  bool KAON=false;
  bool PION=false;
  int t0segment=-1;
  int DEFsegment=-1;
  // double sddenergy[32],sddtrailing[32],sddleading[32];
  // for(int i=0;i<32;i++) sddenergy[i]=sddtrailing[i]=sddleading[i]=-999;
  {
    // data type
    const int nhodo=9;
    int Cid[5]={DetIdT0new,DetIdBHD,DetIdT0,DetIdE0,DetIdDEF}; 
    DetectorType hodo[nhodo]={kT0new, kBHD,kT0,kE0,kDEF,kSTART,kSTOP,kSDDVeto, kSDD1};
    TString hodo_name[nhodo]={"T0new","BHD","T0","E0","DEF","START","STOP","SDDVeto","SDD"};
    const int nsegs[nhodo]={1,16,5,3,3,8,5,4,32};
    const int nud[nhodo]={2,2,2,2,2,2,2,1,1};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};

    
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
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
		  if(tdc>1600&&tdc<1700){
		    beamtiming[ud]=true;
		  }
		  break;
		// case kSDD1:
		//   sddleading[seg]=tdc;
		default:
		  if(tdc>2.3e6&&tdc<2.5e6){
		    beamtiming[ud]=true;
		  }
		}
	      }
	      if(data[idata]==k_adc&&gUnpacker.get_entries(k_device, 0, seg, ud, k_leading)){
		if(hodo[i]==kSDD1){
		  hid  = gHist.getSequentialID(kDET, ud, kADCwTDC, seg+1);
		  hptr_array[hid]->Fill(tdc);	  		
		  continue;
		}
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
	}
      }//seg
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 0);
      hptr_array[hid]->Fill(multiplicity);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 1);
      hptr_array[hid]->Fill(mul[0]);	    
      hid  = gHist.getSequentialID(kDET, 0, kMulti, 2);
      hptr_array[hid]->Fill(mul[1]);	    
    } //hodo
    for(int ihodo=0;ihodo<5;++ihodo){
      DetectorType kDET=hodo[ihodo];
      int nh = hodoAna->GetNHits(Cid[ihodo]);
      //    std::cout<<name[ihodo]<<"  "<<nh<<std::endl;
      //    HF1( Hid[ihodo]+10, double(nh) );
      for( int i=0; i<nh; ++i ){
	Hodo2Hit *hit = hodoAna->GetHit(Cid[ihodo],i);
	HodoRawHit *raw = hit->GetRawHit();
	if(!hit) continue;
	int seg = hit->SegmentId()+1;
	hid  = gHist.getSequentialID(kDET, 0, kCTime,  seg);
	for(int it=0;it<hit->GetIndex();it++){
	  double mt  = hit->MeanTime(it); 
	  hptr_array[hid]->Fill(mt);
	  if(kDET==kT0new&&mt>-25&&mt<-15) time0=mt;
	  if(kDET==kT0&&mt>-70&&mt<-60){
	    t0segment=seg;	    
	  }
	  if(kDET==kDEF&&mt>-70&&mt<-60){
	    DEFsegment=seg;	    
	  }
	  if(kDET==kBHD){
	    double tof=time0-mt;
	    hid  = gHist.getSequentialID(kAna, 0, 1,  1); 
	    hptr_array[hid]->Fill(tof);
	    if(tof>509&&tof<514) KAON=true;
	    else if(tof<507&&tof>503) PION=true;
	  }	    
	}
      }//for(i)
    }    
  }
  
  // Aerogel 
  {
    double ac_sum=0;
    DetectorType kDET=kAC;
    const int k_device = gUnpacker.get_device_id("AC");
    const int nch = 4;
    int nhit=gUnpacker.get_entries(k_device, 0, 0, 0, 1);
    bool ACHIT=false;
    if(nhit>0){
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0,0, 1 , m);
	hid  = gHist.getSequentialID(kDET, 0, kTDC,  1);
	if(tdc<420&&tdc>360) ACHIT=true;
	hptr_array[hid]->Fill(tdc);	    
      }
    }
    for(int ch=0; ch<nch; ++ch){
      //      for(int idata=0; idata<ndata; ++idata){
      int nhit = gUnpacker.get_entries(k_device, 0, 0, ch, 0);
      for( int m=0; m<nhit; ++m ){
	int tdc = gUnpacker.get(k_device, 0, 0, ch, 0 , m);
	hid  = gHist.getSequentialID(kDET, 0, kADC, ch+1);
	//  std::cout<< "seg:"<<seg <<" ,ud:"<<ud<<" ,data:"<<idata<<" ,hid:"<<hid<<" ,val: "<<tdc<<std::endl;
	ac_sum+=tdc;
	hptr_array[hid]->Fill(tdc);	    
	if(ACHIT){	   
	  hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, ch+1);
	  hptr_array[hid]->Fill(tdc);	    
	} 
      } // nhit
	//      } // data
    }//ch
    hid  = gHist.getSequentialID(kDET, 0, kADC, 5);
    hptr_array[hid]->Fill(ac_sum);	    
    if(ACHIT){
      hid  = gHist.getSequentialID(kDET, 0, kADCwTDC, 5);
      hptr_array[hid]->Fill(ac_sum);	    
    }
  }


  // Chamber ------------------------------------------------------------
  const int nchm=7;
  int ChmCid[nchm]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdSDC,0,0};
  DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b,kSDC, kBPC,kFDC};
  TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b","SDC", "BPC","FDC"};
  const int nwires[nchm]={32,32,32,32,16,15,64};
  const int nlayers[nchm]={8,8,8,8,8,8,6};
  double rotz[nchm]={45,45,-45,-45,0,0,0};
  double ffz[nchm]={0,0,130,130,20,0,0};
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

  {
    for(int ichm=0;ichm<5;ichm++){
      DetectorType kDET=chm[ichm];
      int cid=ChmCid[ichm];
      int nl=nlayers[ichm];
      int nw=nwires[ichm];    
      int xy[8]={0,0,1,1,0,0,1,1};
      bool Single=true;
      for( int layer=1; layer<=nl; ++layer ){
	const DCHitContainer &cont = DCAna->GetDCHC(cid, layer);
	int mul=cont.size();
	//	std::cout<<ichm<<"  "<<layer<<"  "<<mul<<std::endl;
	if(mul!=1) Single=false;
      }
      LocalTrack *track=new LocalTrack();
      for( int layer=1; layer<=nl; ++layer ){
	const DCHitContainer &cont = DCAna->GetDCHC(cid, layer);
	int mul=cont.size();
	for(int ihit=0;ihit<mul;ihit++){
	  DCHit* hit=cont[ihit];      	
	  if(mul==1){
	    track->AddHit(hit,xy[layer-1]);
	  }
	}
      }
      if(Single){
	track->DoFit();
	track->ConvLocalToGlobal(rotz[ichm]);
	double x,y;
	track->XYPosatZ(0,x,y);

	hid=gHist.getSequentialID(kDET,0,kProf,1);
	hptr_array[hid]->Fill(x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,2);
	hptr_array[hid]->Fill(x);
	hid=gHist.getSequentialID(kDET,0,kProf,3);
	hptr_array[hid]->Fill(y);
	if(KAON){
	  hid=gHist.getSequentialID(kDET,0,kProf,4);
	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,5);
	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,6);
	  hptr_array[hid]->Fill(y);
	  if(DEFsegment>0){
	    hid=gHist.getSequentialID(kDET,0,kProf,21);
	    hptr_array[hid]->Fill(x,y);
	    hid=gHist.getSequentialID(kDET,0,kProf,22);
	    hptr_array[hid]->Fill(x);
	    hid=gHist.getSequentialID(kDET,0,kProf,23);
	    hptr_array[hid]->Fill(y);
	  }
	  // if(t0segment>0){
	  //   hid=gHist.getSequentialID(kAna,0,2,t0segment);
	  //   hptr_array[hid]->Fill(x,y);
	  //   // hid=gHist.getSequentialID(kDET,0,kProf,30+t0segment);
	  //   // hptr_array[hid]->Fill(x);
	  //   // hid=gHist.getSequentialID(kDET,0,kProf,40+t0segment);
	  //   // hptr_array[hid]->Fill(y);
	  // }
	}else if(PION){
	  hid=gHist.getSequentialID(kDET,0,kProf,7);
	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,8);
	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,9);
	  hptr_array[hid]->Fill(y);
	}

	track->XYPosatZ(ffz[ichm],x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,11);
	hptr_array[hid]->Fill(x,y);
	hid=gHist.getSequentialID(kDET,0,kProf,12);
	hptr_array[hid]->Fill(x);
	hid=gHist.getSequentialID(kDET,0,kProf,13);
	hptr_array[hid]->Fill(y);
	if(KAON){
	  hid=gHist.getSequentialID(kDET,0,kProf,14);
	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,15);
	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,16);
	  hptr_array[hid]->Fill(y);
	}else if(PION){
	  hid=gHist.getSequentialID(kDET,0,kProf,17);
	  hptr_array[hid]->Fill(x,y);
	  hid=gHist.getSequentialID(kDET,0,kProf,18);
	  hptr_array[hid]->Fill(x);
	  hid=gHist.getSequentialID(kDET,0,kProf,19);
	  hptr_array[hid]->Fill(y);
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

