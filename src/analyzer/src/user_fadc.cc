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
//#include "MatrixParamMan.hh"

#define DEBUG    0
#define FLAG_DAQ 1

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    TFile*outRootFile;
    TTree*rawevt;
    //TTree*Header;
    
    std::vector<std::vector<double>> trvecADC;
    //std::vector<std::vector<unsigned int>> trvecADC;
    unsigned int CLKCounter,TIGArCounter;
    
    //int eventcounter=0;
    const int nch=6;
    TString outputname="t98_tmp.root";
  }

//____________________________________________________________________________
  int
  process_begin( const std::vector<std::string>& argv )
  {
    ConfMan& gConfMan = ConfMan::getInstance();
    ConfMan::getInstance().initialize(argv);
    //gUser.Print();

    if(argv.size()==4){
      //int runnumber = stoi(argv.at(3));
      //outputname=Form("/home/oper/t98_ctrl/data/run_%05d.root",runnumber); 
      outputname=argv.at(3);
    }

    outRootFile = new TFile(outputname,"recreate");  
    outRootFile->SetCompressionAlgorithm(5);  
    rawevt = new TTree("rawevt",  "rawevt"); 
    rawevt -> Branch("trvecADC"   ,    &trvecADC); 
    rawevt -> Branch("CLKCounter" ,    &CLKCounter);  
    rawevt -> Branch("TIGArCounter" ,  &TIGArCounter);  
    rawevt->SetDirectory(0); 
    std::cout<<__func__<<std::endl;

    // unpacker and all the parameter managers are initialized at this stage
    
    // Macro for HttpServer
     std::cout<<"gHttp.Register: QDC1"<<std::endl;
    // gHttp.Register(http::QDC(kQDC1,"QDC1",0,0,16,4,4,0,4096));


    return 0;
  }
 
 int
  process_end( void )
  {
    outRootFile->cd();
    // Header->Print();  
    rawevt->Print();  
    // Header->Write();  
    rawevt->Write();  
    
    outRootFile->Close(); 
    return 0;
  }

  //____________________________________________________________________________
  int
  process_event( void )
  {
    static UnpackerManager& gUnpacker = GUnpacker::get_instance();
    // static HistMaker&     gHist     = HistMaker::getInstance();
#if DEBUG
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif

    // GRAMS40 ------------------------------------------------------------
    //int size=gUnpacker.get_entries(gUnpacker.get_device_id("GRAMS40"), 0, 0, 0, gUnpacker.get_data_id("GRAMS40", "fadc"));
    int size = 40000;
    std::cout << "Data Size : " << size << std::endl;
    trvecADC.resize(6,std::vector<double>(size));
    //trvecADC.resize(6,std::vector<unsigned int>(size));
     std::cout << "FADC Value : " << std::endl;

    static const int k_device   = gUnpacker.get_device_id("GRAMS40");
    static const int k_fadc      = gUnpacker.get_data_id("GRAMS40", "fadc");
    int nbin,nbin0;
    if(gUnpacker.get_counter()!=0){
      for(int ich=0;ich <nch;ich++){
	nbin = gUnpacker.get_entries(k_device, 0, 0, ich, k_fadc);
        //std::cout << "Nbin : " << nbin << std::endl;
 
	if(ich==0){
	  nbin0=nbin;
	  //printf("nbin0=%d\n",nbin0);
	}
	if(nbin==0) continue;
	for (int ibin=0; ibin < nbin; ibin++) {
	  unsigned int fadc = gUnpacker.get(k_device, 0, 0, ich, k_fadc ,ibin);
	  //std::cout << "FADC Value : " << fadc << std::endl;
	  //if(fadc==0xffff)continue;
	  //std::cout << "Event : " <<  gUnpacker.get_counter() << " Ch : " << ich << "  Bin  : "<< ibin << std::endl;
	  trvecADC[ich][ibin] = fadc;
          //if(fadc!=65535)std::cout << "FADC Value : " << fadc << std::endl;
	  //printf("Ev%d, nbin%d trvecADC[%d][%d]=%d\n",gUnpacker.get_counter(),nbin ,ich,ibin,int(trvecADC[ich][ibin]));
	}
      }
      CLKCounter = gUnpacker.get_counter();

      Double_t tmpPedVal=0;
      int tmps=0, StartTimeBit=0;
      for(int ibin=0;ibin<gUnpacker.get_data_id("GRAMS40", "fadc");ibin++){
        if(trvecADC[0][ibin]!=0xFFFF){
          tmps=ibin;
          break;
        } 
      }
      tmpPedVal=trvecADC[0][tmps];
      //for(int ibin=tmps;ibin<tmps+10;ibin++)tmpPedVal+=0.1*trvecADC[0][ibin];
      //printf("tmpPedVal=%f, nbin0-1=%d\n",tmpPedVal,nbin0-1);
      for(Int_t ibin=nbin0-1; ibin>0; --ibin){
        if(trvecADC[0][ibin]!=0&&trvecADC[0][ibin]<-3000+tmpPedVal){
          StartTimeBit = ibin-(21*0.4*4);
          //printf("  StartTimeBit=%d\n",StartTimeBit);
          break;
        }
      }
      if(StartTimeBit<=0) TIGArCounter = 4294967295;
      else{
      	TIGArCounter = 0;
      	for(Int_t ibit=0; ibit<32; ++ibit){
      	  if(trvecADC[0][StartTimeBit-(Int_t)(ibit*12.5*0.4*4)]<-2000+tmpPedVal &&
      	     trvecADC[0][StartTimeBit-(Int_t)(ibit*12.5*0.4*4)]!=0&&
             trvecADC[0][StartTimeBit-(Int_t)(ibit*12.5*0.4*4)]!=0xffff){
      	    TIGArCounter+=pow(2,ibit);
      	  }
      	}
      }
      TIGArCounter--;
      //printf("  TIGArCounter=%d\n",TIGArCounter);
      rawevt -> Fill();
      std::cout << "Event : " <<  gUnpacker.get_counter() << std::endl; 
      //if(gUnpacker.get_counter()==0)std::cout << "Event : " <<  gUnpacker.get_counter() << std::endl; 
        
    }

    //
    // // QDC&PADC ------------------------------------------------------------
    // {
    //   // data type
    //   const int nadc=1;
    //   DetectorType adc[nadc]={kQDC1};
    //   TString adc_name[nadc]={"QDC1"};
    //   const int nsegs[nadc]={32};
    //   for(int i=0;i<nadc;i++){
    //     DetectorType kDET=adc[i];
    //     const int k_device = gUnpacker.get_device_id(adc_name[i].Data());
    //     for(int l = 0; l<nsegs[i]; ++l){
    // 	int nhit = gUnpacker.get_entries(k_device, 0, l, 0, 0);
    // 	//
    // 	if( nhit ){
    // 	  for( int m=0; m<nhit; ++m ){
    // 	    int tdc = gUnpacker.get(k_device, 0, l, 0, 0 , m);
    // 	    int hid  = gHist.getSequentialID(kDET, 0, kADC, l+1);
    // 	    // std::cout<<adc_name[i]<<"  "<<l<<"  "
    // 	    // 	     <<m<<" / "<<nhit<<"  "<<tdc<<std::endl;
    // 	    hptr_array[hid]->Fill(tdc);
    // 	  }
    // 	}
    //     }    
    //   }
    // }
    // FADC ------------------------------------------------------------
    // {
    //   // data type
    //   const int nadc=1;
    //   DetectorType adc[nadc]={kFADC};
    //   TString adc_name[nadc]={"FADC"};
    //   const int nsegs[nadc]={32};
    //   for(int i=0;i<nadc;i++){
    //     DetectorType kDET=adc[i];
    //     const int k_device = gUnpacker.get_device_id(adc_name[i].Data());
    //     for(int l = 0; l<nsegs[i]; ++l){
    // 	int nhit = gUnpacker.get_entries(k_device, 0, l, 0, 0);
    // 	//
    // 	if( nhit ){
    // 	  for( int m=0; m<nhit; ++m ){
    // 	    int tdc = gUnpacker.get(k_device, 0, l, 0, 0 , m);
    // 	    int hid  = gHist.getSequentialID(kDET, 0, kADC, l+1);
    // 	    // std::cout<<adc_name[i]<<"  "<<l<<"  "
    // 	    // 	     <<m<<" / "<<nhit<<"  "<<tdc<<std::endl;
    // 	    hptr_array[hid]->Fill(tdc);
    // 	  }
    // 	}
    //     }    
    //   }
    // }

    //   //------------------------------------------------------------------
    //   // BGO
    //   //------------------------------------------------------------------
    //   {
    //     // data type
    //     static const int k_device   = gUnpacker.get_device_id("BGO");
    //     //    static const int k_adc      = gUnpacker.get_data_id("BGO", "adc");
    //     //    static const int k_tdc      = gUnpacker.get_data_id("BGO", "tdc");
    //     static const int k_fadc      = gUnpacker.get_data_id("BGO", "fadc");
    //     static const int k_leading   = gUnpacker.get_data_id("BGO", "leading");
    // //    static const int k_trailing   = gUnpacker.get_data_id("BGO", "trailing");

    //     // TDC gate range
    //     static const unsigned int tdc_min = gUser.GetParameter("BGO_TDC", 0);
    //     static const unsigned int tdc_max = gUser.GetParameter("BGO_TDC", 1);

    //     int bgo_fa_id = gHist.getSequentialID(kBGO, 0, kFADC);
    //     int bgo_fawt_id = gHist.getSequentialID(kBGO, 0, kADCwTDC);
    //     int bgo_t_id = gHist.getSequentialID(kBGO, 0, kTDC);
    //     int bgo_hit_id  = gHist.getSequentialID(kBGO, 0, kHitPat,  1);
    //     int bgo_chit_id = gHist.getSequentialID(kBGO, 0, kHitPat,  2);
    //     int bgo_mul_id  = gHist.getSequentialID(kBGO, 0, kMulti,   1);
    //     int bgo_cmul_id = gHist.getSequentialID(kBGO, 0, kMulti,   2);
    //     unsigned int multiplicity  = 0;
    //     unsigned int cmultiplicity = 0;

    //     for(int seg=0; seg<NumOfSegBGO; ++seg){
    //       // FADC
    //       int nhit_f = gUnpacker.get_entries(k_device, 0, seg, 0, k_fadc);
    //       if(nhit_f==0) continue;
    //       for (int i=0; i<nhit_f; ++i) {
    // 	unsigned int fadc = gUnpacker.get(k_device, 0, seg, 0, k_fadc ,i);
    // 	hptr_array[bgo_fa_id + seg]->Fill( i+1, fadc);
    //       }


    //       // TDC && Hit pattern && multiplicity
    //       unsigned int nhit_l = gUnpacker.get_entries(k_device, 0, seg, 0, k_leading);
    // //      unsigned int nhit_t = gUnpacker.get_entries(k_device, 0, seg, 0, k_trailing);
    // //      unsigned int hit_l_max = 0;
    // //      unsigned int hit_t_max = 0;
    // //
    // //      if(nhit_l != 0){
    // //        hit_l_max = gUnpacker.get(k_device, 0, seg, 0, k_leading,  nhit_l - 1);
    // //      }
    // //      if(nhit_t != 0){
    // //        hit_t_max = gUnpacker.get(k_device, 0, seg, 0, k_trailing, nhit_t - 1);
    // //      }

    //       if(nhit_l!=0){
    //         for(unsigned int m = 0; m<nhit_l; ++m){
    //       	unsigned int tdc = gUnpacker.get(k_device, 0, seg, 0, k_leading, m);
    //       	  if(tdc!=0){
    //       	    hptr_array[bgo_t_id + seg]->Fill(tdc);
    // 	    hptr_array[bgo_hit_id]->Fill(seg);
    // 	    ++multiplicity;
    // 	    if(true && tdc_min < tdc && tdc < tdc_max){
    //               hptr_array[bgo_chit_id]->Fill(seg);
    // 	      ++cmultiplicity;
    // 	    }

    //       	    // ADC wTDC
    //             for (int i=0; i<nhit_f; ++i) {
    //               unsigned int fadc = gUnpacker.get(k_device, 0, seg, 0, k_fadc ,i);
    //               hptr_array[bgo_fawt_id + seg]->Fill( i+1, fadc);
    //             }
    //       	  }
    //         }
    //       }
    //     }

    //     for(int seg=0; seg<NumOfSegBGO_T; ++seg){
    //       unsigned int nhit_l = gUnpacker.get_entries(k_device, 1, seg, 0, k_leading);
    //       nhit_l = gUnpacker.get_entries(k_device, 1, seg, 0, k_leading);
    //       if(nhit_l!=0){
    //         for(unsigned int m = 0; m<nhit_l; ++m){
    //       	unsigned int tdc = gUnpacker.get(k_device, 1, seg, 0, k_leading, m);
    //       	  if(tdc!=0){
    //       	    hptr_array[bgo_t_id + NumOfSegBGO + seg]->Fill(tdc);
    //       	  }
    //         }
    //       }
    //     }

    //     hptr_array[bgo_mul_id]->Fill(multiplicity);
    //     hptr_array[bgo_cmul_id]->Fill(cmultiplicity); // CMulti

#if 0
    // Debug, dump data relating this detector
    //gUnpacker.dump_data_device(k_device);
#endif
    //} BGO


#if DEBUG
    std::cout << __FILE__ << " " << __LINE__ << std::endl;
#endif
    //  gUnpacker.dump_data_fe(521);
    gSystem->ProcessEvents();
    return 0;
  }
  
}
