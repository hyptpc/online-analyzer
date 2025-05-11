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
    const int nchm=7;
    DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b, kBPC1, kBPC2, kCDC};
    TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b", "BPC1", "BPC2",  "CDC"};
    const int nwires[nchm]={32,32,32,32,16,32,16};
    const int nlayers[nchm]={8,8,8,8,8,8,118};
    // for CDC analysis
    std::map< int, int > CDCmap;
    std::map< int, std::pair<int,int> > VFTmap;   
   
    // hodoscopes
    const int nhodo=10;
    DetectorType hodo[nhodo]={kT0,kT1,kDEF,kCDH,kVeto,kPbF2,kPbG,kBTC,kCVC,kNC};
    TString hodo_name[nhodo]={"T0","T1","DEF","CDH","Veto","PbF2","PbG","BTC","CVC","NC"};
    const int nsegs[nhodo]={5,1,4,36,4,36,40,2,10,6};
    const int nud[nhodo]={2,2,2,2,2,1,1,2,2,2};
    const int k_adc=0;
    const int k_leading=1;
    const int k_trailing=2;
    const int ndata=3;
    const int data[ndata]={k_leading,k_trailing,k_adc};
    const int type[ndata]={kTDC,kTDC2D,kADC};

    //    int nbinsqdc=512;
    int nbinsqdc=1024;
    int nbinshrtdc=int(2e3);
    double hrtdcmin=0.3e6;
    double hrtdcmax=0.9e6;    
    double hrtdcmin2=0.0e6;
    double hrtdcmax2=1.0e6;    

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
  if(argv.size()==5){
    outputname=argv.at(3);
    port=std::strtoull(argv.at(4).c_str(), nullptr, 10);
  }
  
  std::cout << "aa" <<std::endl;
  gHttp.SetPort(port);
  gHttp.Open();

  // Hodoscopes
  //  gHttp.Register(gHist.createHodo(kBHD		,"BHT"  ,16, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createBHT(kBHT   ,"BHT"  ,63, 2, nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kT0	,"T0"   , 5, 2, nbinsqdc/4,-0.5,511.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kT1	,"T1"   , 1, 2, nbinsqdc/4,-0.5,1023.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kAC	,"AC"   , 5, 1, nbinsqdc/2,-0.5,2047.5,2000,0,2000));
  gHttp.Register(gHist.createHodo(kDEF	,"DEF"  , 4, 2, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kCDH	,"CDH"  ,36, 2, nbinsqdc/2,-0.5,2047.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kVeto	,"Veto" , 4, 2, nbinsqdc/2,-0.5,1023.5,nbinshrtdc,hrtdcmin,hrtdcmax));
  gHttp.Register(gHist.createHodo(kPbF2	,"PbF2" ,36, 1, nbinsqdc/4,-0.5,511.5,nbinshrtdc,hrtdcmin,hrtdcmax)); 
  gHttp.Register(gHist.createHodo(kPbG	,"PbG"  ,40, 1, nbinsqdc/2,-0.5,1023.5,nbinshrtdc,hrtdcmin,hrtdcmax)); 
  gHttp.Register(gHist.createHodo(kBTC	,"BTC"  , 2, 2, nbinsqdc/4,-0.5,511.5,nbinshrtdc,hrtdcmin,hrtdcmax));
    gHttp.Register(gHist.createHodo(kCVC	,"CVC"  ,10, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin2,hrtdcmax2));
  gHttp.Register(gHist.createHodo(kNC	,"NC"   , 6, 2, nbinsqdc,-0.5,2047.5,nbinshrtdc,hrtdcmin2,hrtdcmax2)); 
  //  gHttp.Register(gHist.createHodo(kCNCtest	,"CNCtest"   , 4, 2, nbinsqdc,-0.5,4095.5,nbinshrtdc,hrtdcmin-0.2e6,hrtdcmax)); 
  gHttp.Register(gHist.createMHTDC(kTriggerFlag ,"TriggerFlag",32,2000));
  // Chambers
  gHttp.Register(gHist.createBLDC(kBLC1a,"BLC1a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC1b,"BLC1b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2a,"BLC2a",8,32,true));
  gHttp.Register(gHist.createBLDC(kBLC2b,"BLC2b",8,32,true));
  gHttp.Register(gHist.createBLDC(kBPC1,"BPC1",8,16,true));
  gHttp.Register(gHist.createBLDC(kBPC2,"BPC2",8,32,true));
  //  gHttp.Register(gHist.createBLDC(kVFT,"VFT",4,224,true));
  gHttp.Register(gHist.createBLDC(kVFT,"VFT",7,128,true));
  gHttp.Register(gHist.createCDC(true));
  gHttp.Register(gHist.createBLDC(kCDC,"CDC",118,16,true));

  if(0 != gHist.setHistPtr(hptr_array)){ return -1; }
  
  // Macro for HttpServer
  // BHT
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_U",0,63,8,8),"BHT");
  gHttp.Register(http::MHTDCTDC(kBHT,"_BHT_D",1,63,8,8),"BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_U",63,0,8,8), "BHT");
  gHttp.Register(http::BHTTOT(kBHT,"_D",63,1,8,8), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_U",63,0,8,8), "BHT");
  gHttp.Register(http::BHTTDCvsTOT(kBHT,"_D",63,1,8,8), "BHT");
  gHttp.Register(http::MHTDCHitPatMulti(kBHT,"_BHT"),"BHT");

  // Hodoscopes
  gHttp.Register(http::QDC(kT0,"_T0_U",0,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::QDC(kT0,"_T0_D",1,0,5,3,2,0,1000),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_U",0,5,3,2),"T0");
  gHttp.Register(http::MHTDCTDC(kT0,"_T0_D",1,5,3,2),"T0");
  gHttp.Register(http::MHTDCHitPatMulti(kT0,"_T0"),"T0");

  gHttp.Register(http::QDC(kT1,"_T1_U",0,0,1,1,1,0,1000),"T1");
  gHttp.Register(http::QDC(kT1,"_T1_D",1,0,1,1,1,0,1000),"T1");
  gHttp.Register(http::MHTDCTDC(kT1,"_T1_U",0,1,1,1),"T1");
  gHttp.Register(http::MHTDCTDC(kT1,"_T1_D",1,1,1,1),"T1");
  gHttp.Register(http::MHTDCHitPatMulti(kT1,"_T1"),"T1");

  gHttp.Register(http::QDC(kAC,"_AC",0,0,5,3,2,0,2000),"AC");
  gHttp.Register(http::MHTDCTDC(kAC,"_AC",0,1,1,1),"AC");

  gHttp.Register(http::QDC(kDEF,"_DEF_U",0,0,4,2,2,0,2000),"DEF");
  gHttp.Register(http::QDC(kDEF,"_DEF_D",1,0,4,2,2,0,2000),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_U",0,4,2,2),"DEF");
  gHttp.Register(http::MHTDCTDC(kDEF,"_DEF_D",1,4,2,2),"DEF");
  gHttp.Register(http::MHTDCHitPatMulti(kDEF,"_DEF"),"DEF");

  gHttp.Register(http::QDC(kCDH,"_CDH_U",0,0,36,6,6,0,1000),"CDH");
  gHttp.Register(http::QDC(kCDH,"_CDH_D",1,0,36,6,6,0,1000),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_U",0,36,6,6),"CDH");
  gHttp.Register(http::MHTDCTDC(kCDH,"_CDH_D",1,36,6,6),"CDH");
  gHttp.Register(http::MHTDCHitPatMulti(kCDH,"_CDH"),"CDH");

  gHttp.Register(http::QDC(kVeto,"_Veto_L",0,0,4,2,2,0,1000),"Veto");
  gHttp.Register(http::QDC(kVeto,"_Veto_R",1,0,4,2,2,0,1000),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_L",0,4,2,2),"Veto");
  gHttp.Register(http::MHTDCTDC(kVeto,"_Veto_R",1,4,2,2),"Veto");
  gHttp.Register(http::MHTDCHitPatMulti(kVeto,"_Veto"),"Veto");

  gHttp.Register(http::QDC(kBTC,"_BTC_L",0,0,2,2,2,0,1000),"BTC");
  gHttp.Register(http::QDC(kBTC,"_BTC_R",1,0,2,2,2,0,1000),"BTC");
  gHttp.Register(http::MHTDCTDC(kBTC,"_BTC_L",0,2,2,2),"BTC");
  gHttp.Register(http::MHTDCTDC(kBTC,"_BTC_R",1,2,2,2),"BTC");
  gHttp.Register(http::MHTDCHitPatMulti(kBTC,"_BTC"),"BTC");

  gHttp.Register(http::QDC(kCVC,"_CVC_U",0,0,10,4,3,0,1000),"CVC");
  gHttp.Register(http::QDC(kCVC,"_CVC_D",1,0,10,4,3,0,1000),"CVC");
  gHttp.Register(http::MHTDCTDC(kCVC,"_CVC_U",0,10,4,3),"CVC");
  gHttp.Register(http::MHTDCTDC(kCVC,"_CVC_D",1,10,4,3),"CVC");
  gHttp.Register(http::MHTDCHitPatMulti(kCVC,"_CVC"),"CVC");

  gHttp.Register(http::QDC(kNC,"_NC_U",0,0,6,3,2,0,1000),"NC");
  gHttp.Register(http::QDC(kNC,"_NC_D",1,0,6,3,2,0,1000),"NC");
  gHttp.Register(http::MHTDCTDC(kNC,"_NC_U",0,6,3,2),"NC");
  gHttp.Register(http::MHTDCTDC(kNC,"_NC_D",1,6,3,2),"NC");
  gHttp.Register(http::MHTDCHitPatMulti(kNC,"_NC"),"NC");

  // gHttp.Register(http::QDC(kCNCtest,"_CNCtest_U",0,0,4,2,2,0,1000),"CNCtest");
  // gHttp.Register(http::QDC(kCNCtest,"_CNCtest_D",1,0,4,2,2,0,1000),"CNCtest");
  // gHttp.Register(http::MHTDCTDC(kCNCtest,"_CNCtest_U",0,4,2,2),"CNCtest");
  // gHttp.Register(http::MHTDCTDC(kCNCtest,"_CNCtest_D",1,4,2,2),"CNCtest");
  // gHttp.Register(http::MHTDCHitPatMulti(kCNCtest,"_CNCtest"),"CNCtest");

  gHttp.Register(http::QDC(kPbF2,"_PbF2",0,0,36,6,6,0,500),"PbF2");
  gHttp.Register(http::MHTDCTDC(kPbF2,"_PbF2",0,36,6,6),"PbF2");
  gHttp.Register(http::MHTDCHitPatMulti(kPbF2,"_PbF2"),"PbF2");

  gHttp.Register(http::PbG(kPbG,"_PbG",0,0,40,7,7,0,2000),"PbG");
  gHttp.Register(http::MHTDCTDC(kPbG,"_PbG",0,40,8,5),"PbG");
  gHttp.Register(http::MHTDCHitPatMulti(kPbG,"_PbG"),"PbG");

  // Chambers except for CDC
  for(int i=0;i<nchm-1;i++){
    const std::string tmpstr=Form("_%s",chm_name[i].Data());
    gHttp.Register(http::BLDCHitPat(chm[i],tmpstr,8),chm_name[i]); 
    gHttp.Register(http::BLDCMulti(chm[i],tmpstr,8),chm_name[i]);
    gHttp.Register(http::BLDCTDC(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTOT(chm[i],tmpstr,8), chm_name[i]);
    gHttp.Register(http::BLDCTDCvsTOT(chm[i],tmpstr,8), chm_name[i]);
    for(int l=0;l<8;l++){ 
      gHttp.Register(http::BLDCWIRE(chm[i],tmpstr,l,nwires[i],8,4), chm_name[i]);
    }
  }
  // TriggerFlag
  gHttp.Register(http::MHTDCTDC(kTriggerFlag,"_TriggerFlag",0,32,8,4),"TriggerFlag");
  gHttp.Register(http::MHTDCHitPatMulti(kTriggerFlag,"_TriggerFlag"),"TriggerFlag");
  {
    int hid1 = gHist.getSequentialID(kTriggerFlag, 0, kHitPat, 1);  
    hptr_array[hid1]->GetXaxis()->SetTitle("");
    for( Int_t i=0; i<32; ++i ){
      int hid2 = gHist.getSequentialID(kTriggerFlag, 0, kTDC, i+1);
      hptr_array[hid2]->SetTitle(Form("%s_%s",hptr_array[hid2]->GetTitle(),flagnames[i].Data()));
      hptr_array[hid1]->GetXaxis()->SetBinLabel(i+1,flagnames[i]);
    }
  }		   
  // for VFT analysis
  // gHttp.Register(http::BLDCHitPat(  kVFT,"VFT",4,0,2,2), "VFT");
  // gHttp.Register(http::BLDCMulti(   kVFT,"VFT",4,0,2,2), "VFT");
  // gHttp.Register(http::BLDCTDC(     kVFT,"VFT",4,0,2,2), "VFT");
  // gHttp.Register(http::BLDCTOT(     kVFT,"VFT",4,0,2,2), "VFT");
  // gHttp.Register(http::BLDCTDCvsTOT(kVFT,"VFT",4,0,2,2), "VFT");
  // for VFT analysis
  gHttp.Register(http::BLDCHitPat(  kVFT,"VFT",7,0,4,2), "VFT");
  gHttp.Register(http::BLDCMulti(   kVFT,"VFT",7,0,4,2), "VFT");
  gHttp.Register(http::BLDCTDC(     kVFT,"VFT",7,0,4,2), "VFT");
  gHttp.Register(http::BLDCTOT(     kVFT,"VFT",7,0,4,2), "VFT");
  gHttp.Register(http::BLDCTDCvsTOT(kVFT,"VFT",7,0,4,2), "VFT");

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

  std::ifstream data2("/home/oper/online_br/VFTmap.txt");
  if (data2.fail()) {
    std::cerr << "Failed to open VFT-map file." << std::endl;
    return -1;
  }
  int fiber; // 0 origin
  int easirocnum,easirocch; // 0 origin
  while(!data2.eof()){
    data2 >> easirocnum >> easirocch >> layer >> fiber;
    VFTmap[1000*easirocnum+easirocch] = std::make_pair(layer,fiber);
  }

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

  static Int_t run_number = -1;
  if(run_number != gUnpacker.get_run_number()){
    for(Int_t i=0, n=hptr_array.size(); i<n; ++i){
      hptr_array[i]->Reset();
    }
    run_number = gUnpacker.get_run_number();
  }
  auto event_number = gUnpacker.get_event_number();

  int hid;
  bool COSMIC=false;
  bool CLOCK=false;
  // TriggerFlag
  {
    // data type
    hid=-1;
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
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	for(int idata=0; idata<ndata; ++idata){
	  int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, data[idata]);
	  if(data[idata]==k_leading&&nhit>0){
	    hid  = gHist.getSequentialID(kDET, 0, kHitPat, 1);
	    hptr_array[hid]->Fill(seg);	    
            if(seg==14) COSMIC=true;
            if(seg==15) CLOCK=true;
	  }
	  for( int m=0; m<nhit; ++m ){
	    int tdc = gUnpacker.get(k_device, 0, seg, ud, data[idata] , m);
	    hid  = gHist.getSequentialID(kDET, ud, type[idata], seg+1);
	    hptr_array[hid]->Fill(tdc);	    
	  } // nhit
	} // data
      }// seg
    } // triggerflag
  }
  //if(COSMIC&&!CLOCK) return 0;
  // BHT ------------------------------------------------------------
  {
    // data type    
    DetectorType kDET=kBHT;
    const int k_l=0;
    const int k_t=1;
    const int k_device = gUnpacker.get_device_id("BHT");
    int multiplicity=0;
    int mul[2]={0,0};
    int segud[2]={};
    for(int seg = 0; seg<63; ++seg){
      int ntdc[2]={0,0};
      for(int ud=0; ud<2; ++ud){	  
	int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, k_l);
	if( nhit==0 ) continue;
	segud[ud]=seg;
	// This wire fired at least one times.
	++multiplicity;
	mul[ud]++;
	ntdc[ud]=nhit;	  
	std::vector< int > leading_array;
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, k_l, m);
	  hid=gHist.getSequentialID(kDET,ud,kTDC,seg+1);
	  //	  if(ud==0&&seg==0) std::cout<<seg<<"  "<<m<<"  "<<tdc<<"  "<<hid<<std::endl;
	  hptr_array[hid]->Fill(tdc);
	  leading_array.push_back(tdc);
          if(tdc>1.15e6&&tdc<1.2e6){
            hid=gHist.getSequentialID(kDET,0,kHitPat,ud+1);
            hptr_array[hid]->Fill(seg);
          }
	}	
	// traling
	nhit = gUnpacker.get_entries(k_device, 0, seg, ud, k_t);
	if( nhit==0 ) continue;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, 0, seg, ud, k_t , m);
	  hid=gHist.getSequentialID(kDET,ud,kTDC2D,seg+1);
	  hptr_array[hid]->Fill(tdc);
	  // for TOT
	  if(m<leading_size){
	    hid=gHist.getSequentialID(kDET,ud,kTOT,seg+1);
	    hptr_array[hid]->Fill(leading_array.at(m) - tdc);
	    hid=gHist.getSequentialID(kDET,ud,kADC2D,seg+1);
	    hptr_array[hid]->Fill(leading_array.at(m), leading_array.at(m) - tdc);
	  }
	}
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
    if(multiplicity>0){
      hid  = gHist.getSequentialID(kDET, 0, kHitPat2D, 0); 
      hptr_array[hid]->Fill(segud[0],segud[1]);	    
    }
  } //bht

  // Hodoscope ------------------------------------------------------------
  {
    // data type    
    for(int i=0;i<nhodo;i++){
      DetectorType kDET=hodo[i];
      const int k_device = gUnpacker.get_device_id(hodo_name[i].Data());
      int multiplicity=0;
      int mul[2]={0,0};
      int segud[2]={-1,-1};
      for(int seg = 0; seg<nsegs[i]; ++seg){
	int ntdc[2]={0,0};
	for(int ud=0; ud<nud[i]; ++ud){	  
	  for(int idata=0; idata<ndata; ++idata){
	    int nhit = gUnpacker.get_entries(k_device, 0, seg, ud, data[idata]);
	    if(data[idata]==k_leading&&nhit>0){
	      hid  = gHist.getSequentialID(kDET, 0, kHitPat, ud+1);
	      hptr_array[hid]->Fill(seg);	    
	      segud[ud]=seg;
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
      if(multiplicity>0){
	hid  = gHist.getSequentialID(kDET, 0, kHitPat2D, 0); 
	hptr_array[hid]->Fill(segud[0],segud[1]);	    
      }
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
        //	if(tdc<6.7e5&&tdc>6.3e5) ACHIT=true;
	if(tdc<1060&&tdc>1020) ACHIT=true;
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

	// trailing
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
	    hid=gHist.getSequentialID(kDET,0,kTOT,l+1);
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

  // VFT
  {
    // data type
    DetectorType kDET=kVFT;
    const int k_device = gUnpacker.get_device_id("VFT");
    const int k_l    = gUnpacker.get_data_id("VFT" , "leading");
    const int k_t    = gUnpacker.get_data_id("VFT" , "trailing");
    // TDC & HitPat & Multi
    TString hname;
    int hid=-1;
    std::pair<int,int> map;
    for(int l = 0; l<7; ++l){ // hul mhtdc number
      for(int w = 0; w<128; ++w){ // hul chan
	int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	if( nhit==0 ) continue;
	hid=gHist.getSequentialID(kDET,0,kHitPat,l+1);
	hptr_array[hid]->Fill(w);
	std::vector< int > leading_array(nhit);
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	  hid=gHist.getSequentialID(kDET,0,kTDC,l+1);
	  hptr_array[hid]->Fill(tdc);
	  leading_array[m] = tdc;
	}
	nhit = gUnpacker.get_entries(k_device, l, 0, w, k_t);
	if( nhit==0 ) continue;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_t , m);
	  hid=gHist.getSequentialID(kDET,0,kTDC2D,l+1);
	  hptr_array[hid]->Fill(tdc);
	  if(m<leading_size){
	    hid=gHist.getSequentialID(kDET,0,kTOT,l+1);
	    hptr_array[hid]->Fill(leading_array[m] - tdc);	    
	    hid=gHist.getSequentialID(kDET,0,kADC2D,l+1);
	    hptr_array[hid]->Fill(leading_array[m],leading_array[m] - tdc);
	  }
	} 
      }
    }
  }

#if 0
  // VFT
  {
    // data type
    DetectorType kDET=kVFT;
    const int k_device = gUnpacker.get_device_id("VFT");
    const int k_l    = gUnpacker.get_data_id("VFT" , "leading");
    const int k_t    = gUnpacker.get_data_id("VFT" , "trailing");
    int mulVFT[4] = {0,0,0,0};
    // TDC & HitPat & Multi
    TString hname;
    int hid=-1;
    std::pair<int,int> map;
    for(int l = 0; l<14; ++l){ // easiroc number
      for(int w = 0; w<64; ++w){ // easiroc num chan
        if( !VFTmap.count(1000*l+w) ) continue;
        map = VFTmap[1000*l+w];
	int layer = map.first;
	int fiber  = map.second;
	int nhit = gUnpacker.get_entries(k_device, l, 0, w, k_l);
	if( nhit==0 ) continue;
	mulVFT[layer]++;
	hid=gHist.getSequentialID(kDET,0,kHitPat,layer+1);
	hptr_array[hid]->Fill(fiber);
	std::vector< int > leading_array(nhit);
	int leading_size = nhit;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_l, m);
	  hid=gHist.getSequentialID(kDET,0,kTDC,layer+1);
	  hptr_array[hid]->Fill(tdc);
	  leading_array[m] = tdc;
	}
	nhit = gUnpacker.get_entries(k_device, l, 0, w, k_t);
	if( nhit==0 ) continue;
	for( int m=0; m<nhit; ++m ){
	  int tdc = gUnpacker.get(k_device, l, 0, w, k_t , m);
	  hid=gHist.getSequentialID(kDET,0,kTDC2D,layer+1);
	  hptr_array[hid]->Fill(tdc);
	  if(m<leading_size){
	    hid=gHist.getSequentialID(kDET,0,kTOT,layer+1);
	    hptr_array[hid]->Fill(leading_array[m] - tdc);	    
	    hid=gHist.getSequentialID(kDET,0,kADC2D,layer+1);
	    hptr_array[hid]->Fill(leading_array[m],leading_array[m] - tdc);
	  }
	} 
      }
      for( int i=0; i<4; i++ ){
        hid=gHist.getSequentialID(kDET,0,kMulti,i+1);
        hptr_array[hid]->Fill(mulVFT[i]);
      }
    }
  }
#endif
#if DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

  gSystem->ProcessEvents();
  return 0;
}

}

