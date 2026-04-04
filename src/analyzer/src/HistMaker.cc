// -*- C++ -*-

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <utility>

#include "DetectorID.hh"
#include "HistHelper.hh"
#include "HistMaker.hh"

#include "DCGeomMan.hh"
#include "DetSizeMan.hh"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TDirectory.h>
#include <TString.h>
#include <TH2Poly.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "TpcPadHelper.hh"

#include <Unpacker.hh>
#include <UnpackerManager.hh>

ClassImp(HistMaker)

// getStr_FromEnum ----------------------------------------------------------
// The method to get std::string from enum value
#define CONV_STRING(x) getStr_FromEnum(#x)
std::string getStr_FromEnum(const char* c){
  std::string str = c;
  return str.substr(1);
}

namespace
{
  const std::string& MyName("HistMaker::");
  const DCGeomMan& gGeom  = DCGeomMan::GetInstance();
  const DetSizeMan& gSize = DetSizeMan::GetInstance();
  const auto& gTpcPad  = TpcPadHelper::GetInstance();
  using hddaq::unpacker::GUnpacker;
  const auto& gUnpacker = GUnpacker::get_instance();
}

// Constructor -------------------------------------------------------------
HistMaker::HistMaker( void )
  : current_hist_id_(0)
{
}

// Destructor -------------------------------------------------------------
HistMaker::~HistMaker( void )
{
}

// -------------------------------------------------------------------------
// getListOfDetectors
// -------------------------------------------------------------------------
void HistMaker::getListOfDetectors(std::vector<std::string>& vec)
{
  HistMaker& g = HistMaker::getInstance();
  std::copy(g.name_created_detectors_.begin(),
	    g.name_created_detectors_.end(),
	    back_inserter(vec)
	    );
}

// -------------------------------------------------------------------------
// getListOfPsFiles
// -------------------------------------------------------------------------
void HistMaker::getListOfPsFiles(std::vector<std::string>& vec)
{
  HistMaker& g = HistMaker::getInstance();
  std::copy(g.name_ps_files_.begin(),
	    g.name_ps_files_.end(),
	    back_inserter(vec)
	    );

}

// -------------------------------------------------------------------------
// setHistPtr
// -------------------------------------------------------------------------
int HistMaker::setHistPtr(std::vector<TH1*>& vec)
{
  static const std::string MyFunc = "setHistPtr ";

  vec.resize(current_hist_id_);
  for(int i = 0; i<current_hist_id_; ++i){
    int unique_id = getUniqueID(i);
    vec[i] = GHist::get(unique_id);
    if(vec[i] == NULL){
      std::cerr << "#E: " << MyName << MyFunc
		<< "Pointer is NULL\n"
		<< " Unique ID    : " << unique_id << "\n"
		<< " Sequential ID: " << i << std::endl;
      gDirectory->ls();
      return -1;
    }
  }

  return 0;
}

// -------------------------------------------------------------------------
// CreateTH1
// -------------------------------------------------------------------------
TH1* HistMaker::createTH1(int unique_id, const char* title,
			  int nbinx, double xmin, double xmax,
			  const char* xtitle, const char* ytitle
			  )
{
  static const std::string MyFunc = "createTH1 ";

  // Add information to dictionaries which will be used to find the histogram
  int sequential_id = current_hist_id_++;
  TypeRetInsert ret =
    idmap_seq_from_unique_.insert(std::make_pair(unique_id,     sequential_id));
  idmap_seq_from_name_.insert(  std::make_pair(title,         sequential_id));
  idmap_unique_from_seq_.insert(std::make_pair(sequential_id, unique_id));
  if(!ret.second){
    std::cerr << "#E: " << MyName << MyFunc
	      << "The unique id overlaps with other id"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
  }

  // create histogram using the static method of HistHelper class
  //  std::cout<< "#D: " << MyFunc << " " << title << " " << unique_id << std::endl;
  TH1 *h = GHist::D1(unique_id, title, nbinx, xmin, xmax);
  if(!h){
    std::cerr << "#E: " << MyName << MyFunc
	      << "Fail to create TH1"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
    //    return h;
  }

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  return h;
}

// -------------------------------------------------------------------------
// CreateTH2
// -------------------------------------------------------------------------
TH2* HistMaker::createTH2(int unique_id, const char* title,
			  int nbinx, double xmin, double xmax,
			  int nbiny, double ymin, double ymax,
			  const char* xtitle, const char* ytitle
			  )
{
  static const std::string MyFunc = "createTH2 ";

  // Add information to dictionaries which will be used to find the histogram
  int sequential_id = current_hist_id_++;
  TypeRetInsert ret =
    idmap_seq_from_unique_.insert(std::make_pair(unique_id,     sequential_id));
  idmap_seq_from_name_.insert(  std::make_pair(title,         sequential_id));
  idmap_unique_from_seq_.insert(std::make_pair(sequential_id, unique_id));
  if(!ret.second){
    std::cerr << "#E: " << MyName << MyFunc
	      << "The unique id overlaps with other id"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
  }
  
  // create histogram using the static method of HistHelper class
  TH2 *h = GHist::I2(unique_id, title,
		     nbinx, xmin, xmax,
		     nbiny, ymin, ymax);
  if(!h){
    std::cerr << "#E: " << MyName << MyFunc
	      << "Fail to create TH2"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
    //    return h;
  }

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  return h;
}

//_____________________________________________________________________________
TH2* HistMaker::createTH2Poly( Int_t unique_id, const TString& title,
                               Double_t xmin, Double_t xmax,
                               Double_t ymin, Double_t ymax )
{
  
  static const std::string MyFunc = "createTH2Poly ";
  Int_t sequential_id = current_hist_id_++;
  TypeRetInsert ret =
    idmap_seq_from_unique_.insert( std::make_pair( unique_id, sequential_id ) );
  idmap_seq_from_name_.insert( std::make_pair( title, sequential_id ) );
  idmap_unique_from_seq_.insert( std::make_pair( sequential_id, unique_id ) );
  if( !ret.second ){
    std::cerr << "#E " << MyName << MyFunc
	      << " The unique id overlaps with other id"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
  }

  TH2 *h = GHist::P2( unique_id, title, xmin, xmax, ymin, ymax );
  if( !h ){
    std::cerr << "#E " << MyName << " Fail to create TH2" << std::endl
	      << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
    //    return h;
  } else {
    gDirectory->GetList()->Add( h );
  }
  return h;
}

//_____________________________________________________________________________
TH3* HistMaker::createTH3( Int_t unique_id, const TString& title,
			   int nbinx, Double_t xmin, Double_t xmax,
			   int nbiny, Double_t ymin, Double_t ymax,
			   int nbinz, Double_t zmin, Double_t zmax)
{
  
  static const std::string MyFunc = "createTH3 ";
  Int_t sequential_id = current_hist_id_++;
  TypeRetInsert ret =
    idmap_seq_from_unique_.insert( std::make_pair( unique_id, sequential_id ) );
  idmap_seq_from_name_.insert( std::make_pair( title, sequential_id ) );
  idmap_unique_from_seq_.insert( std::make_pair( sequential_id, unique_id ) );
  if( !ret.second ){
    std::cerr << "#E " << MyName << MyFunc
	      << " The unique id overlaps with other id"
	      << std::endl;
    std::cerr << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
  }

  TH3 *h = GHist::D3( unique_id, title, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax );
  if( !h ){
    std::cerr << "#E " << MyName << " Fail to create TH3" << std::endl
	      << " " << unique_id << " " << title << std::endl;
    std::exit(-1);
    //    return h;
  } else {
    gDirectory->GetList()->Add( h );
  }
  return h;
}


//_____________________________________________________________________
TList*
HistMaker::createTimeStamp( bool flag_ps )
{
  // Determine the detector name
  std::string strDet = CONV_STRING(kTimeStamp);
  // name list of crearted detector
  name_created_detectors_.push_back(strDet);

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  {
    // Make histogram and add it
    // Make unique ID
    int target_id = getUniqueID( kTimeStamp, 0, kTDC, 0);
    for( int i=0; i<NumOfPlaneVmeRm; ++i ){
      int seg = i+1; // 1 origin
      const char* title = Form("%s_%d", nameDetector, seg);
      top_dir->Add( createTH1( ++target_id, title, // 1 origin
			       0x1000, -0x1000, 0x1000,
			       "TimeStamp", ""));
    }
  }

  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createTriggerFlag(Bool_t flag_ps)
{
  std::string strDet = "TriggerFlag";
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  { // TDC
    Int_t target_id = getUniqueID(kTriggerFlag, 0, kTDC, 0);
    for(Int_t i = 0; i<NumOfSegTFlag; ++i){
      TString title = Form("%s #%d", trigger::STriggerFlag[i].Data(), i);
      top_dir->Add(createTH1(++target_id, title,
			     1000, 0, 2000000,
			     "TDC [ch]", ""));
    }
  }
  { // HitPat
    const char* title = "TriggerFlag_HitPat";
    Int_t target_id = getUniqueID(kTriggerFlag, 0, kHitPat);
    auto h = createTH1(target_id, title,
		       NumOfSegTFlag, 0., NumOfSegTFlag,
		       "", "");
    for(Int_t i=0, n=trigger::STriggerFlag.size(); i<n; ++i){
      h->GetXaxis()->SetBinLabel(i+1, trigger::STriggerFlag.at(i));
    }
    h->SetStats(0);
    top_dir->Add(h);
  }
  return top_dir;
}




//_____________________________________________________________________________
TList*
HistMaker::createBcInTracking(Bool_t flag_ps)
{
  const std::string strDet("BcInTracking");
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  TString title;

  std::vector<std::pair<DataType, std::string>> particle_types = {
    {kPi, "Pi"},
    {kKaon, "K"},
    {kAll, "All"}
  };
  {
    for (auto& [type, particle] : particle_types) {
      Int_t target_id = getUniqueID(kBcInTracking, 0, type, 0);
      title = Form("%s_X_%s", nameDetector,particle.c_str());
      top_dir->Add(createTH1(++target_id, title,
			     300, -300, 300,
			     "X [mm]", ""));

      title = Form("%s_Y_%s", nameDetector,particle.c_str());
      top_dir->Add(createTH1(++target_id, title,
			     300, -300, 300,
			     "Y [mm]", ""));
      title = Form("%s_2D_%s", nameDetector,particle.c_str());
      top_dir->Add( createTH2(++target_id, title,
			      100, -300, 300,
			      100, -300, 300,
			      "X [mm]", "Y [mm]" ) );
    }
  }
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createBcOutTracking(Bool_t flag_ps)
{
  const std::string strDet("BcOutTracking");
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  TString title;
  std::vector<std::pair<DataType, std::string>> particle_types = {
    {kPi, "Pi"},
    {kKaon, "K"},
    {kAll, "All"}
  };
  {
    for (auto& [type, particle] : particle_types) {
      Int_t target_id = getUniqueID(kBcOutTracking, 0, type, 0);
      title = Form("%s_X_%s", nameDetector,particle.c_str());
      top_dir->Add(createTH1(++target_id, title,
			     200, -200, 200,
			     "X [mm]", ""));

      title = Form("%s_Y_%s", nameDetector,particle.c_str());
      top_dir->Add(createTH1(++target_id, title,
			     200, -200, 200,
			     "Y [mm]", ""));
      title = Form("%s_BH2_%s", nameDetector,particle.c_str());
      top_dir->Add( createTH2(++target_id, title,
			      100, -200, 200,
			      100, -200, 200,
			      "X [mm]", "Y [mm]" ) );

      title = Form("%s_HTOF_window_%s", nameDetector,particle.c_str());
      top_dir->Add( createTH2(++target_id, title,
			      100, -200, 200,
			      100, -200, 200,
			      "X [mm]", "Y [mm]" ) );

      title = Form("%s_Target_%s", nameDetector,particle.c_str());
      top_dir->Add( createTH2(++target_id, title,
			      100, -200, 200,
			      100, -200, 200,
			      "X [mm]", "Y [mm]" ) );

    }
  }
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createTPC(Bool_t flag_ps)
{

  const std::string strDet = "TPC";
  name_created_detectors_.push_back(strDet);
  if( flag_ps ) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  {
    Int_t target_id = getUniqueID( kTPC, 0, kADC2D );
    auto title = Form( "%s_ADC2D", nameDetector );
    auto h_adc = dynamic_cast<TH2Poly*>
      ( createTH2Poly( target_id++, title, -300, 300, -300, 300 ) );
    title = Form( "%s_RMS2D", nameDetector );
    auto h_rms = dynamic_cast<TH2Poly*>
      ( createTH2Poly( target_id++, title, -300, 300, -300, 300 ) );
    title = Form( "%s_TDC2D", nameDetector );
    auto h_loc = dynamic_cast<TH2Poly*>
      ( createTH2Poly( target_id++, title, -300, 300, -300, 300 ) );
    title = Form( "%s_Hit2D", nameDetector );
    auto h_hit = dynamic_cast<TH2Poly*>
      ( createTH2Poly( target_id++, title, -300, 300, -300, 300 ) );
    title = Form( "%s_Threshold2D", nameDetector );
    auto h_thre = dynamic_cast<TH2Poly*>
      ( createTH2Poly( target_id++, title, -300, 300, -300, 300 ) );
    top_dir->Add( createTH3( getUniqueID(kTPC, 0, kTDC3D), "TPC_TDC3D", 100,-300,300,100,-300,300,250,-250,250));


    Double_t X[5]={};
    Double_t Y[5]={};
    for( Int_t i=0; i<32; i++ ){
      Double_t pLength = TpcPadHelper::PadParameter[i][5];
      Double_t st      = ( 180. -
                           360./TpcPadHelper::PadParameter[i][3] *
                           TpcPadHelper::PadParameter[i][1]/2. );
      Double_t sTheta  = (-1+st/180.)*TMath::Pi();
      Double_t dTheta  = ( (360./TpcPadHelper::PadParameter[i][3]) /
                           180.*TMath::Pi() );
      Double_t cRad    = TpcPadHelper::PadParameter[i][2];
      Int_t    nPad    = TpcPadHelper::PadParameter[i][1];
      for( Int_t j=0; j<nPad; j++ ){
	X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
	X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
	X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
	X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
	X[0] = X[4];
	Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
	Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
	Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
	Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
	Y[0] = Y[4];
	for( Int_t ii=0; ii<5; ii++ ) X[ii] -=143;
        h_adc->AddBin( 5, X, Y );
        h_rms->AddBin( 5, X, Y );
        h_loc->AddBin( 5, X, Y );
        h_hit->AddBin( 5, X, Y );
	h_thre->AddBin( 5, X, Y );
      }
    }
    
    h_adc->SetStats( 0 );
    h_adc->SetMinimum(    0. );
    h_adc->SetMaximum( 4000. );
    h_rms->SetStats( 0 );
    h_rms->SetMinimum(    0. );
    h_rms->SetMaximum(  200. );
    h_loc->SetStats( 0 );
    h_loc->SetMinimum(    0. );
    h_loc->SetMaximum(  NumOfTimeBucket );
    h_thre->SetMinimum(   0. );
    //h_thre->SetMaximum( 2000. );
    top_dir->Add( h_adc );
    top_dir->Add( h_rms );
    top_dir->Add( h_loc );
    top_dir->Add( h_hit );
    top_dir->Add( h_thre );
    
  }
  // ADC
  top_dir->Add( createTH1( getUniqueID( kTPC, 0, kADC, 1 ),
                           "TPC_ADC", 4000, 0, 4000 ) );
  top_dir->Add( createTH1( getUniqueID( kTPC, 0, kPede ),
                           "TPC_RMS", 1000, 0, 1000 ) );
  // TDC
  top_dir->Add( createTH1( getUniqueID( kTPC, 0, kTDC ),
                           "TPC_TDC", NumOfTimeBucket, 0, NumOfTimeBucket ) );
  top_dir->Add( createTH2( getUniqueID( kTPC, 0, kTDC2D ),
                           "TPC_Hit_ZY", 100, -300, 300,
			   NumOfTimeBucket, 0, NumOfTimeBucket,
			   "z [mm]", "tb" ) );
  top_dir->Add( createTH2( getUniqueID( kTPC, 1, kTDC2D ),
                           "TPC_Hit_XY", 100, -300, 300,
			   NumOfTimeBucket, 0, NumOfTimeBucket,
			   "x [mm]", "tb" ) );
  top_dir->Add( createTH2( getUniqueID( kTPC, 2, kTDC2D ),
                           "TPC_ZY", 100, -300, 300,
			   NumOfTimeBucket, 0, NumOfTimeBucket,
			   "z [mm]", "tb" ) );
  top_dir->Add( createTH2( getUniqueID( kTPC, 3, kTDC2D ),
                           "TPC_XY", 100, -300, 300,
			   NumOfTimeBucket, 0, NumOfTimeBucket,
			   "x [mm]", "tb" ) );
  // FADC
  top_dir->Add( createTH2( getUniqueID( kTPC, 0, kFADC,1 ),
                           "TPC_FADC",
                           NumOfTimeBucket, 0, NumOfTimeBucket,
			   200, 0, 0x1000,
                           "Time bucket", "ADC" ) );

  // FADC per AGET
  for(int n_asad=0;n_asad < NumOfAsadTPC;n_asad++){
    for(int n_aget = 0;n_aget<4;n_aget++){
      top_dir->Add( createTH1( getUniqueID( kTPC, 0, kFADC, 2+n_asad*4+n_aget),
			       Form("TPC_FADC_AsAd%d_AGET%d",n_asad+1,n_aget),
			       NumOfTimeBucket, -0.5, NumOfTimeBucket-0.5));
    }
  }

  
  // Multiplicity
  top_dir->Add( createTH1( getUniqueID( kTPC, 0, kMulti ),
                           "TPC_multiplicity", 600, 0, 6000 ) );
  top_dir->Add( createTH1( getUniqueID( kTPC, 3, kMulti ),
                           "TPC_AGET_multiplicity", 124, 0, 124,
			   "AsAdID#times4+AGETID", "Multiplicity/AGET/Event" ) );
  top_dir->Add( createTH1( getUniqueID( kTPC, 4, kMulti ),
                           "TPC_AGET_multiplicity_Max", 64, 0, 64 ) );


  top_dir->Add( createTH1( getUniqueID( kTPC, 1, kADC),
                           "TPC_AGET_Count", 124, 0,124,
			   "AsAdID#times4+AGETID", "Count" ) );
  
  // ClusterSize
  top_dir->Add( createTH2( getUniqueID( kTPC, 2, kMulti ),
                           "TPC_ClusterSize", 42, -10, 32, 10, 0, 10,
			   "Layer ID", "Cluster size") );
  // TPC-CLOCK
  top_dir->Add( createTH1( getUniqueID( kTPC, 1, kTDC ),
			   "TPC_CLOCK", 50000, 0, 2000000,
			   "TDC", "") );
  top_dir->Add( createTH1( getUniqueID( kTPC, 1, kMulti ),
			  "TPC_CLOCK_multiplicity", 10, 0, 10) );
  // TPC-BeamProfile
//  top_dir->Add( createTH1( getUniqueID( kTPC, 2, kTDC ),
//		 	  "TPC_BeamProfile", 34, 0, 34,
//		 	  "Pad", "") );
  top_dir->Add( createTH1( getUniqueID( kTPC, 2, kTDC ),
			   "TPC_BeamProfile", 20, -100, 100,
			   "x [mm]", "") );
  return top_dir;
}


//_____________________________________________________________________________
TList *
HistMaker::createBTOF(Bool_t flag_ps)
{
  TString det_name = "BTOF";
  auto top_dir = new TList;
  top_dir->SetName(det_name);
  const auto hid = getUniqueID(kMisc, 0, kTDC);

  top_dir->Add(createTH1(hid, det_name,
                         300, -20, 20, "[ns]", ""));
  top_dir->Add(createTH2(hid+1,"BAC vs BTOF",300,-20,20,100,0,4000,"BTOF [ns]","ADC"));
    
  /*
  top_dir->Add(createTH1(hid + 1, "BH1-6_BH2-4",
                         600, 50000, 350000, "[ch]", ""));
  
  top_dir->Add(createTH1(hid + 2, det_name + "_wide",
                         1800, -90, 5, "[ns]", ""));
  */
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createCorrelation(Bool_t flag_ps)
{
  std::string strDet = "Correlation";
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  { 
    Int_t target_id = getUniqueID(kCorrelation, 0, kHitPat2D, 0);
    TString title = Form("%s BH2%BHT", nameDetector);
    top_dir->Add(createTH2(++target_id, title,
			   NumOfSegBHT, -0.5, NumOfSegBHT-0.5,
			   NumOfSegBH2, -0.5, NumOfSegBH2-0.5,
			   "BHT segment", "BH2 segment"));
    title = Form("%s HTOF-Mp2", nameDetector);
    top_dir->Add(createTH2(++target_id, title,
			   NumOfSegHTOF, -0.5, NumOfSegHTOF-0.5,
			   NumOfSegHTOF, -0.5, NumOfSegHTOF-0.5,
			   "HTOF segment1", "HTOF segment2"));
    title = Form("%s BH2%%HTOF", nameDetector);
    top_dir->Add(createTH2(++target_id, title,
			   NumOfSegHTOF, -0.5, NumOfSegHTOF-0.5,
			   NumOfSegBH2, -0.5, NumOfSegBH2-0.5,
			   "HTOF segment", "BH2 segment"));
    title = Form("%s KVC%%HTOF", nameDetector);
    top_dir->Add(createTH2(++target_id, title,
			   NumOfSegHTOF, -0.5, NumOfSegHTOF-0.5,
			   NumOfSegKVC, -0.5, NumOfSegKVC-0.5,
			   "HTOF segment", "KVC segment"));
  }
  return top_dir;
}

TList*
HistMaker::createEventDisplay(Bool_t flag_ps)
{
  

  //TargetID order (kEventDisplay, kHitPoly)
  //1  : HTOF Hit pattern
  //2  : HTOF Hit count
  //3  : BH2 Hit pattern
  //4  : BH2 Hit count
  //5  : BAC Hit pattern
  //6  : BAC Hit count
  //7  : KVC Hit pattern
  //8  : KVC Hit count
  //9  : T2 Hit pattern
  //10 : T2 Hit count
  //11 : SCH Hit pattern
  //12 : SCH Hit count
  //13 : BcOut
  
  std::string strDet = "EventDisplay";
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName("EventDisplay");

  top_dir->SetName(nameDetector);
  Int_t target_id = getUniqueID(kEventDisplay, 0, kHitPoly, 0);

  const int tpc_id = gGeom.GetDetectorId("HypTPC");
  const double tpc_z = gGeom.GetGlobalPosition(tpc_id).z();
  
  
  TString title;

  
  //HTOF
  {
    title = Form("%s_HTOF_HitPatternPoly", nameDetector);
    
    auto htof_pattern_poly = dynamic_cast<TH2Poly*>
      ( createTH2Poly( ++target_id, title, c_x_min,c_x_max,c_y_min,c_y_max ) );
    title = Form( "%s_HTOF_HitCountPoly", nameDetector );
    auto htof_count_poly = dynamic_cast<TH2Poly*>
      ( createTH2Poly( ++target_id, title, c_x_min,c_x_max,c_y_min,c_y_max) );

    const auto& htof_size = gSize.GetSize("HTOF");
    
    const Double_t w = htof_size.x();
    const Double_t L = 337;
    const Double_t t = htof_size.z();

    Double_t theta[8];
    Double_t X[5];
    Double_t Y[5];
    Double_t seg_X[5];
    Double_t seg_Y[5];
    for( Int_t i=0; i<8; i++ ){
      theta[i] = (-180+45*i)*acos(-1)/180.;
      for( Int_t j=0; j<4; j++ ){
	seg_X[1] = L-t/2.;
	seg_X[2] = L+t/2.;
	seg_X[3] = L+t/2.;
	seg_X[4] = L-t/2.;
	seg_X[0] = seg_X[4];
	seg_Y[1] = w*j-2*w;
	seg_Y[2] = w*j-2*w;
	seg_Y[3] = w*j-1*w;
	seg_Y[4] = w*j-1*w;
	seg_Y[0] = seg_Y[4];
	for( Int_t k=0; k<5; k++ ){
	  X[k] = cos(theta[i])*seg_X[k]-sin(theta[i])*seg_Y[k];
	  Y[k] = sin(theta[i])*seg_X[k]+cos(theta[i])*seg_Y[k];
	}
	htof_pattern_poly->AddBin( 5, X, Y );
	htof_count_poly->AddBin( 5, X, Y );
      }
    }
    htof_pattern_poly->SetStats( 0 );
    htof_pattern_poly->SetMinimum( 0. );
    htof_pattern_poly->SetMaximum( 170. );

    htof_count_poly->SetStats( 0 );
    htof_count_poly->SetMinimum( 0. );
    
    top_dir->Add( htof_pattern_poly );
    top_dir->Add( htof_count_poly);
  }
  
  for (int i = 1; i < sizeof(EvtDis_Det_name)/sizeof(EvtDis_Det_name[0]); i++) {
    auto name = EvtDis_Det_name[i];
    title = Form("%s_%s_HitPatternPoly",nameDetector,name.Data());
    auto h_pattern_poly = dynamic_cast<TH2Poly*>
      ( createTH2Poly(++target_id, title, c_x_min,c_x_max,c_y_min,c_y_max) );
    title = Form("%s_%s_HitCountPoly", nameDetector, name.Data());
    auto h_count_poly = dynamic_cast<TH2Poly*>
      ( createTH2Poly(++target_id, title, c_x_min,c_x_max,c_y_min,c_y_max) );
    const auto& det_size = gSize.GetSize(name.Data());
    const Double_t w = det_size.x();
    Double_t t = det_size.z();
    
    if(name == "BH2")t*=2;
    if(name == "T2")t*=2;
    
    
    const int det_id = gGeom.GetDetectorId(name.Data());
    auto det_pos = gGeom.GetGlobalPosition(det_id);
    const Double_t det_x = det_pos.x();
    Double_t det_z = det_pos.z() - tpc_z;
    if(name == "SCH"){
      t*=4;
      //det_z = 800;
    }
    
    double zmin = det_z - t/2.;
    double zmax = det_z + t/2.;

    double totalw = NumOfSeg[i]*w;
    double x_start = -totalw / 2.;

    for(Int_t j=0;j<NumOfSeg[i];j++){
      double xmin = x_start + j*w + det_x;
      double xmax = xmin + w;
	
      h_pattern_poly->AddBin(zmin,xmin,zmax,xmax);
      h_count_poly->AddBin(zmin,xmin,zmax,xmax);
      
    }
    
    h_pattern_poly->SetStats( 0 );
    h_pattern_poly->SetMinimum( 0. );
    h_pattern_poly->SetMaximum( 200. );

    h_count_poly->SetStats( 0 );
    h_count_poly->SetMinimum( 0. );

    top_dir->Add( h_pattern_poly );
    top_dir->Add( h_count_poly );
    
  }

  //BcOut
  {
    title = Form("%s_BcOut", nameDetector);
    top_dir->Add(createTH2(++target_id,title,
			   500,-900,0,
			   500,-500,0,
			   "Z [mm]","X [mm]"));

  }
  
  return top_dir;
  
}




// -------------------------------------------------------------------------
// createQDC
// -------------------------------------------------------------------------
TList* HistMaker::createQDC(DetectorType kDET, std::string strDet, const int nsegments, int nbins, double xmin, double xmax, bool flag_ps)
{
  name_created_detectors_.push_back(strDet);
  if(flag_ps){
    // name list which are displayed in Ps tab
    name_ps_files_.push_back(strDet);
  }

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kADC,  "ADC",    "ADC [ch]", nbins,xmin,xmax) );
  list.push_back( HistMakerInfo(kADCwTDC,  "ADCwTDC",  "ADC [ch]", nbins, xmin,xmax) );
  size_t size=list.size();
  for(size_t i=0; i<size;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

// -------------------------------------------------------------------------
// createHodo
// -------------------------------------------------------------------------
TList* HistMaker::createHodo(DetectorType kDET, std::string strDet, const int nsegments, const int nud, int nbins, double xmin, double xmax,int nbins2, double xmin2, double xmax2, bool flag_ps)
{
  name_created_detectors_.push_back(strDet);
  if(flag_ps){
    // name list which are displayed in Ps tab
    name_ps_files_.push_back(strDet);
  }

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kADC,    "ADC",    "ADC [ch]", nbins,xmin,xmax) );
  list.push_back( HistMakerInfo(kADCwTDC,"ADCwTDC",  "ADC [ch]", nbins, xmin,xmax) );
  list.push_back( HistMakerInfo(kTDC,    "TDC",       "TDC [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTDC2D,  "Trailing",  "Trainling [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTOT,    "TOT",      "TOT [ch]", 1024,0.,1024.) );
  list.push_back( HistMakerInfo(kCTime,  "MeanTime",       "time [ns]", 5000,-1000,1000) );
  list.push_back( HistMakerInfo(kHitPat, "HitPat",  "Segment", nsegments+1 ,-0.5,nsegments+0.5) );
  list.push_back( HistMakerInfo(kMulti,  "Multi",   "Multiplicity", nsegments+1 ,-0.5,nsegments+0.5) );

  size_t size=list.size();
  for(size_t i=0; i<size-3;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, list[i].type, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      }
    }
    top_dir->Add(sub_dir);
  }  
  // meantime
  for(size_t i=size-3; i<size-2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=size-2; i<size;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    for(int j = 0; j<3; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      top_dir->Add(createTH1(target_id + j, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
  }
  {
    const char* nameSubDir = "HitPat2D";
    int target_id = getUniqueID(kDET, 0, kHitPat2D, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH2(target_id, title, // 1 origin
			   nsegments,-0.5,nsegments-0.5,
			   nsegments,-0.5,nsegments-0.5,
			   "up segment","down segment"));
  }
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

// -------------------------------------------------------------------------
// createSCH
// -------------------------------------------------------------------------
TList* HistMaker::createSCH(DetectorType kDET, std::string strDet, const int nsegments, const int nud, int nbins, double xmin, double xmax, int nbins2, double xmin2, double xmax2, bool flag_ps){

  name_created_detectors_.push_back(strDet);
  if(flag_ps){
    // name list which are displayed in Ps tab
    name_ps_files_.push_back(strDet);
  }

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  //list.push_back( HistMakerInfo(kADC,    "ADC",    "ADC [ch]", nbins,xmin,xmax) );
  //list.push_back( HistMakerInfo(kADCwTDC,"ADCwTDC",  "ADC [ch]", nbins, xmin,xmax) );
  list.push_back( HistMakerInfo(kTDC,    "TDC",       "TDC [ch]", 1024, 0, 1024) );
  list.push_back( HistMakerInfo(kTDC2D,  "Trailing",  "Trainling [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTOT,    "TOT",      "TOT [ch]", 100, 0., 100.) );
  list.push_back( HistMakerInfo(kTOTwTDC,"TOTwTDC",  "TOT [ch]", 100, 0., 100.) );
  list.push_back( HistMakerInfo(kCTime,  "MeanTime",       "time [ns]", 5000,-1000,1000) );
  list.push_back( HistMakerInfo(kHitPat, "HitPat",  "Segment", nsegments+1 ,-0.5,nsegments+0.5) );
  list.push_back( HistMakerInfo(kMulti,  "Multi",   "Multiplicity", nsegments+1 ,-0.5,nsegments+0.5) );
  
  
  size_t size=list.size();
  for(size_t i=0; i<size-3;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, list[i].type, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      }
    }
    top_dir->Add(sub_dir);
  }
  {
    const char* nameSubDir = "TDCvsTOT";
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, kADC2D, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);

	sub_dir->Add(createTH2(target_id + j+1, title,
			       100,xmin2,xmax2,100,0,20000.,
			       "TDC [ch]", "TOT [ch]") );

      }
    }
  }
  
  for(size_t i=size-3; i<size-2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=size-2; i<size;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    for(int j = 0; j<3; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      top_dir->Add(createTH1(target_id + j, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
  }
  {
    const char* nameSubDir = "HitPat2D";
    int target_id = getUniqueID(kDET, 0, kHitPat2D, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH2(target_id, title, // 1 origin
			   nsegments,-0.5,nsegments-0.5,
			   nsegments,-0.5,nsegments-0.5,
			   "top segment","bottom segment"));
  }
  
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

// -------------------------------------------------------------------------
// createBHT
// -------------------------------------------------------------------------
TList* HistMaker::createBHT(DetectorType kDET, std::string strDet, const int nsegments, const int nud, int nbins2, double xmin2, double xmax2, bool flag_ps)
{
  name_created_detectors_.push_back(strDet);
  if(flag_ps){
    // name list which are displayed in Ps tab
    name_ps_files_.push_back(strDet);
  }

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kTDC,    "TDC",       "TDC [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTDC2D,  "Trailing",  "Trainling [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTOT,    "TOT",       "TOT [ch]", 1024,0.,1024.*32) );
  list.push_back( HistMakerInfo(kCTime,  "MeanTime",  "time [ns]", 5000,-1000,1000) );
  list.push_back( HistMakerInfo(kHitPat, "HitPat",  "Segment", nsegments+1 ,-0.5,nsegments+0.5) );
  list.push_back( HistMakerInfo(kMulti,  "Multi",   "Multiplicity", nsegments+1 ,-0.5,nsegments+0.5) );  

  size_t size=list.size();
  for(size_t i=0; i<size-3;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, list[i].type, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      }
    }
    top_dir->Add(sub_dir);
  }  
  {
    const char* nameSubDir = "TDCvsTOT";
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);    
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, kADC2D, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH2(target_id + j+1, title,
			       100,xmin2,xmax2,100,0,20000.,
			       "TDC [ch]", "TOT [ch]") );
      }
    }
  }
  for(size_t i=size-3; i<size-2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=size-2; i<size;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    for(int j = 0; j<3; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      top_dir->Add(createTH1(target_id + j, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
  }
  {
    const char* nameSubDir = "HitPat2D";
    int target_id = getUniqueID(kDET, 0, kHitPat2D, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH2(target_id, title, // 1 origin
			   nsegments,-0.5,nsegments-0.5,
			   nsegments,-0.5,nsegments-0.5,
			   "top segment","bottom segment"));
  }
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

// -------------------------------------------------------------------------
// createHTOF
// -------------------------------------------------------------------------
TList* HistMaker::createHTOF(DetectorType kDET, std::string strDet, const int nsegments, const int nud, int nbins, double xmin, double xmax,int nbins2, double xmin2, double xmax2, bool flag_ps)
{
  name_created_detectors_.push_back(strDet);
  if(flag_ps){
    // name list which are displayed in Ps tab
    name_ps_files_.push_back(strDet);
  }

  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kADC,    "ADC",    "ADC [ch]", nbins,xmin,xmax) );
  list.push_back( HistMakerInfo(kADCwTDC,"ADCwTDC",  "ADC [ch]", nbins, xmin,xmax) );
  list.push_back( HistMakerInfo(kTDC,    "TDC",       "TDC [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTDC2D,  "Trailing",  "Trainling [ch]", nbins2, xmin2, xmax2) );
  list.push_back( HistMakerInfo(kTOT,    "TOT",      "TOT [ch]", 1024,0.,1024.) );
  list.push_back( HistMakerInfo(kCTime,  "MeanTime",       "time [ns]", 5000,-1000,1000) );
  list.push_back( HistMakerInfo(kHitPat, "HitPat",  "Segment", nsegments+1 ,-0.5,nsegments+0.5) );
  list.push_back( HistMakerInfo(kMulti,  "Multi",   "Multiplicity", nsegments+1 ,-0.5,nsegments+0.5) );
  
  
  size_t size=list.size();
  for(size_t i=0; i<size-3;i++){
  // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, list[i].type, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      }
    }
    top_dir->Add(sub_dir);
  }  
  // meantime
  for(size_t i=size-3; i<size-2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=size-2; i<size;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    for(int j = 0; j<4; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      top_dir->Add(createTH1(target_id + j, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
  }
  {
    const char* nameSubDir = "Threshold";
    int target_id = getUniqueID(kDET, 0, kThreshold, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH1(target_id, title, // 1 origin
			   nbins,xmin,xmax,
			   "ADC [ch]",""));
  }
  {
    const char* nameSubDir = "HitPat2D";
    int target_id = getUniqueID(kDET, 0, kHitPat2D, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH2(target_id, title, // 1 origin
			   nsegments,-0.5,nsegments-0.5,
			   nsegments,-0.5,nsegments-0.5,
			   "up segment","down segment"));
  }

  {
    const char* nameSubDir = "TDCvsTOT";
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);    
    for(int iud=0; iud<nud; ++iud){
      int target_id = getUniqueID(kDET, iud, kADC2D, 0);
      for(int j = 0; j<nsegments; ++j){
	const char* title = NULL;
	title = Form("%s_%s_%d_%d", nameDetector, nameSubDir, j, iud);
	sub_dir->Add(createTH2(target_id + j+1, title,
			       100,xmin2,xmax2,100,0,20000.,
			       "TDC [ch]", "TOT [ch]") );
      }
    }
  }
  
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}



// createMHTDC
// -------------------------------------------------------------------------
TList* HistMaker::createMHTDC(DetectorType kDET, std::string strDet, const int nsegments, double xmax, bool flag_ps)
{
  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  int range=2000;
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kTDC,    CONV_STRING(kTDC),    "TDC [ch]", range,0.,xmax) );
  list.push_back( HistMakerInfo(kTDC2D,  CONV_STRING(kTDC2D),  "Trainling [ch]", range,0.,xmax) );
  list.push_back( HistMakerInfo(kADC,    CONV_STRING(kADC),    "TOT [ch]", 1024,0.,1024.) );
  list.push_back( HistMakerInfo(kHitPat, CONV_STRING(kHitPat), "Segment", nsegments+1 ,-0.5,nsegments+0.5) );
  list.push_back( HistMakerInfo(kMulti,  CONV_STRING(kMulti),  "Multiplicity", nsegments+1 ,-0.5,nsegments+0.5) );
  
  size_t size=list.size();
  for(size_t i=0; i<size-2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=size-2; i<size;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    for(int j = 0; j<2; ++j){
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      top_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
    }
  }
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}

// -------------------------------------------------------------------------
// createBLDC
// -------------------------------------------------------------------------
TList* HistMaker::createBLDC(DetectorType kDET, std::string strDet, int nlayers, int nwires, bool WIRE_RAW, bool ANA,bool flag_ps)
{
  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  
  
  double xmax=2000;
  double nbins=xmax;
  double ntot=500;
  // if(kDET==kVFT){
  //   nbins=1024;
  //   xmax=1024;
  //   ntot=200;
  // }
  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kTDC,      CONV_STRING(kTDC),    "TDC [ch]", nbins,0.,xmax) );
  list.push_back( HistMakerInfo(kTDC2D,    CONV_STRING(kTDC2D),  "Trailing [ch]", nbins,0.,xmax) );
  list.push_back( HistMakerInfo(kTOT,      CONV_STRING(kTOT),    "TOT [ch]", ntot,0.,ntot) );
  list.push_back( HistMakerInfo(kHitPat,   CONV_STRING(kHitPat), "Wire", nwires+1 ,-0.5,nwires+0.5) );
  list.push_back( HistMakerInfo(kMulti,    CONV_STRING(kMulti),  "Multiplicity", nwires+1 ,-0.5,nwires+0.5) );
  list.push_back( HistMakerInfo(kMulti2D,  CONV_STRING(kMulti2D), "Multiplicity2", nwires+1 ,-0.5,nwires+0.5) );
  list.push_back( HistMakerInfo(kResid,    CONV_STRING(kResid),  "Residual", 200,-10,10) );
  
  size_t size=list.size();
  for(size_t i=0; i<size;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);
    // Make histogram and add it
    for(int j = 0; j<nlayers; ++j){
      int target_id = getUniqueID(kDET, 0, list[i].type, 0);
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			     list[i].nbinsx, list[i].xmin, list[i].xmax,
			     list[i].xtitle, list[i].ytitle));
      if(i==0&&WIRE_RAW){
	TString nameSubDir2 = Form("layer%d",j);
	TList *sub_dir2 = new TList;
	sub_dir2->SetName(nameSubDir2);
	top_dir->Add(sub_dir2);
	target_id = getUniqueID(kDET, j+1, list[i].type, 0);
	for(int k = 0; k<nwires; ++k){
	  title = Form("%s_%s_%s_%d", nameDetector, nameSubDir, nameSubDir2.Data(), k);
	  sub_dir2->Add(createTH1(target_id + k+1, title, // 1 origin
				  list[i].nbinsx, list[i].xmin, list[i].xmax,
				  list[i].xtitle, list[i].ytitle));
	}
      }
    }
  }  
  
  if(kDET!=kCDC){
    for(int j = 0; j<nlayers; ++j){
      const char* nameSubDir = "TDCvsTOT";
      TList *sub_dir = new TList;
      sub_dir->SetName(nameSubDir);
      top_dir->Add(sub_dir);
      int target_id = getUniqueID(kDET, 0, kADC2D, 0);
      const char* title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
      sub_dir->Add(createTH2(target_id + j+1, title,
                             100,0,xmax,100,0,ntot,
                             "TDC [ch]", "TOT [ch]") );
    }
  }

  if(!ANA) return top_dir;
  { 
    const char* nameSubDir = "Eff";
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, kEff, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH1(target_id, title, // 1 origin
			   nlayers+1,-0.5,nlayers+0.5,
			   "layer", "Efficiency") );
  }

  { 
    const char* nameSubDir = "WireCorr";
    // Make histogram and add it
    for(int i=0;i<4;i++){
      int target_id = getUniqueID(kDET, 0, kWireCorr, i);
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir,i);
      top_dir->Add(createTH2(target_id, title, // 1 origin
			     nwires+1,-0.5,nwires+0.5,nwires+1,-0.5,nwires+0.5,
			     Form("Wire in layer%d",i), Form("Wire in layer%d",i+4)) );
    }
  }
  {
    const char* nameSubDir = "Profile";
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, kProf, 0);

    int nbins=200;
    const char* title = Form("%s_%s_XY_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+1, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_X_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+2, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_Y_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+3, title, // 1 origin
			   nbins,-200,200,
			   "Y") );
    title = Form("%s_%s_XY_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+4, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_X_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+5, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_Y_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+6, title, // 1 origin
			   nbins,-200,200,
			   "Y") );
    title = Form("%s_%s_XY_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+7, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_X_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+8, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_Y_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+9, title, // 1 origin
			   nbins,-200,200,
			   "Y") );    
    title = Form("%s_%s_A_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+10, title, // 1 origin
			   nbins,-0.1,0.1,
			   "A") );

    title = Form("%s_%s_XYatFF_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+11, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_XatFF_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+12, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_YatFF_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+13, title, // 1 origin
			   nbins,-200,200,
			   "Y") );
    title = Form("%s_%s_XYatFF_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+14, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_XatFF_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+15, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_YatFF_Kaon", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+16, title, // 1 origin
			   nbins,-200,200,
			   "Y") );
    title = Form("%s_%s_XYatFF_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+17, title, // 1 origin
			   nbins,-200,200,nbins,-200,200,
			   "X", "Y") );
    title = Form("%s_%s_XatFF_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+18, title, // 1 origin
			   nbins,-200,200,
			   "X") );
    title = Form("%s_%s_YatFF_Pion", nameDetector, nameSubDir);
    sub_dir->Add(createTH1(target_id+19, title, // 1 origin
			   nbins,-200,200,
			   "Y") );    
    if(kDET==kBPC1||kDET==kBPC2){
      title = Form("%s_%s_XYatRC_All", nameDetector, nameSubDir);
      sub_dir->Add(createTH2(target_id+41, title, // 1 origin
			     nbins,-200,200,nbins,-200,200,
			     "X", "Y") );
      title = Form("%s_%s_XatRC_All", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+42, title, // 1 origin
			     nbins,-200,200,
			     "X") );
      title = Form("%s_%s_YatRC_All", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+43, title, // 1 origin
			     nbins,-200,200,
			     "Y") );
      title = Form("%s_%s_XYatRC_Kaon", nameDetector, nameSubDir);
      sub_dir->Add(createTH2(target_id+44, title, // 1 origin
			     nbins,-200,200,nbins,-200,200,
			     "X", "Y") );
      title = Form("%s_%s_XatRC_Kaon", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+45, title, // 1 origin
			     nbins,-200,200,
			     "X") );
      title = Form("%s_%s_YatRC_Kaon", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+46, title, // 1 origin
			     nbins,-200,200,
			     "Y") );
      title = Form("%s_%s_XYatRC_Pion", nameDetector, nameSubDir);
      sub_dir->Add(createTH2(target_id+47, title, // 1 origin
			     nbins,-200,200,nbins,-200,200,
			     "X", "Y") );
      title = Form("%s_%s_XatRC_Pion", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+48, title, // 1 origin
			     nbins,-200,200,
			     "X") );
      title = Form("%s_%s_YatRC_Pion", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+49, title, // 1 origin
			     nbins,-200,200,
			     "Y") );    
	    
      title = Form("%s_%s_XY_DEF", nameDetector, nameSubDir);
      sub_dir->Add(createTH2(target_id+21, title, // 1 origin
			     nbins,-200,200,nbins,-200,200,
			     "X", "Y") );
      title = Form("%s_%s_X_DEF", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+22, title, // 1 origin
			     nbins,-200,200,
			     "X") );
      title = Form("%s_%s_Y_DEF", nameDetector, nameSubDir);
      sub_dir->Add(createTH1(target_id+23, title, // 1 origin
			     nbins,-200,200,
			     "Y") );
    }
    title = Form("%s_%s_AB_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+31, title, // 1 origin
			   nbins,-0.2,0.2,nbins,-0.2,0.2,
			   "A", "B") );
    title = Form("%s_%s_XA_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+32, title, // 1 origin
			   nbins,-200,200,nbins,-0.2,0.2,
			   "X", "A") );
    title = Form("%s_%s_YB_All", nameDetector, nameSubDir);
    sub_dir->Add(createTH2(target_id+33, title, // 1 origin
			   nbins,-200,200,nbins,-0.2,0.2,
			   "Y", "B") );
  }  
  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}


// -------------------------------------------------------------------------
// createDAQ
// -------------------------------------------------------------------------
TList* HistMaker::createDAQ(bool flag_ps)
{
  // Determine the detector name
  std::string strDet = CONV_STRING(kDAQ);
  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  std::vector<Int_t> vme_fe_id;
  std::vector<Int_t> hul_fe_id;
  std::vector<Int_t> ea0c_fe_id;
  std::vector<Int_t> cobo_fe_id;
  for( auto&& c : gUnpacker.get_root()->get_child_list() ){
    if( !c.second )
      continue;
    TString n = c.second->get_name();
    auto id = c.second->get_id();
    if(n.Contains("vme"))
      vme_fe_id.push_back(id);
    if(n.Contains("hul"))
      hul_fe_id.push_back(id);
    if(n.Contains("easiroc"))
      ea0c_fe_id.push_back(id);
    if(n.Contains("cobo"))
      cobo_fe_id.push_back(id);
  }

  // DAQ infomation --------------------------------------------------
  {
    // Event builder infomation
    int target_id;
    {
      target_id = getUniqueID(kDAQ, kEB, kHitPat, 0);
      top_dir->Add(createTH1(target_id + 1, "Data size EB", // 1 origin
			     8000, 0, 80000,
			     "Data size [words]", ""));
    }
    // Node information
    {
      target_id = getUniqueID(kDAQ, kVME, kHitPat2D, 0);
      auto h = createTH2(target_id + 1, "Data size VME nodes", // 1 origin
			 vme_fe_id.size(), 0, vme_fe_id.size(),
			 100, 0, 1200,
			 "VME node ID", "Data size [words]");
      for(Int_t i=0, n=vme_fe_id.size(); i<n; ++i){
	h->GetXaxis()->SetBinLabel( i+1, "0x"+TString::Itoa(vme_fe_id[i], 16));
      }
      top_dir->Add(h);
    }

    //Easiroc
    {
      target_id = getUniqueID(kDAQ, kEASIROC, kHitPat2D, 0);
      top_dir->Add(createTH2(target_id + 1, "Data size EASIROC nodes", // 1 origin
			     2, 0, 2,
			     50, 0, 100,
			     "EASIROC node ID", "Data size [words]"));

    }

    //HUL
    {
      target_id = getUniqueID(kDAQ, kHUL, kHitPat2D, 0);
      auto h = createTH2(target_id + 1, "Data size HUL nodes", // 1 origin
			 hul_fe_id.size(),0,hul_fe_id.size(),
			 200, 0, 400,
			 "HUL node ID", "Data size [words]");
      for(Int_t i=0, n=hul_fe_id.size(); i<n; ++i){
	h->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(hul_fe_id[i], 16));
      }
      top_dir->Add(h);
    }

    //Cobo
    {
      target_id = getUniqueID(kDAQ, kCoBo, kHitPat2D, 0);
      auto h = createTH2(target_id + 1, "Data size CoBo nodes",
			 cobo_fe_id.size(),0, cobo_fe_id.size(),
			 100, 0, 200000,
			 "CoBo node ID", "Data size [words]");
      for(Int_t i=0, n=cobo_fe_id.size(); i<n; ++i){
	h->GetXaxis()->SetBinLabel(i+1, "0x"+TString::Itoa(cobo_fe_id[i], 16));
      }
      top_dir->Add(h);
    }
    
  }
  return top_dir;
}

