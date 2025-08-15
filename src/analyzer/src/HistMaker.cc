// -*- C++ -*-

#include <iostream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <utility>

#include "DetectorID.hh"
#include "HistHelper.hh"
#include "HistMaker.hh"

#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TDirectory.h>
#include <TString.h>

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
HistMaker::createBVH(Bool_t flag_ps)
{
  const std::string strDet("BVH");
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  { // TDC
    Int_t target_id = getUniqueID(kBVH, 0, kTDC, 0);
    for(Int_t i=0; i<NumOfSegBVH; ++i){
      TString title = Form("%s_TDC_%d", nameDetector, i);
      top_dir->Add(createTH1(++target_id, title,
			     1000, 0, 1000,
			     "TDC [ch]", ""));
    }
  }
  { // TOT
    Int_t target_id = getUniqueID(kBVH, 0, kADC, 0);
    for(Int_t i=0; i<NumOfSegBVH; ++i){
      TString title = Form("%s_TOT_%d", nameDetector, i);
      top_dir->Add(createTH1(++target_id, title,
			     200, 0, 200,
			     "TOT [ch]", ""));
    }
  }
  { // HitPat
    auto h = createTH1(getUniqueID(kBVH, 0, kHitPat),
		       "BVH_HitPat",
		       NumOfSegBVH, -0.5, NumOfSegBVH+0.5,
		       "", "");
    top_dir->Add(h);
  }
  { // Multiplicity
    auto h = createTH1(getUniqueID(kBVH, 0, kMulti),
		       "BVH_Multi",
		       NumOfSegBVH+1, -0.5, NumOfSegBVH+0.5,
		       "", "");
    top_dir->Add(h);
  }
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createT1(Bool_t flag_ps)
{
  const std::string strDet = "T1";
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  { ///// ADC
    TString strSubDir  = "ADC";
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT1, 0, kADC, 0);
    for(Int_t i = 0; i<NumOfSegT1*2; ++i){
      const char* title = NULL;
      if(i < NumOfSegT1){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT1; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     0x1000, 0, 0x1000,
			     "ADC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// ADC w/TDC
    TString strSubDir  = CONV_STRING(kADCwTDC);
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT1, 0, kADCwTDC, 0);
    for( Int_t i=0; i<NumOfSegT1*2; ++i ){
      const char* title = NULL;
      if( i<NumOfSegT1 ){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT1; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     0x1000, 0, 0x1000,
			     "ADC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// TDC
    TString strSubDir  = CONV_STRING(kTDC);
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT1, 0, kTDC, 0);
    for(Int_t i = 0; i<NumOfSegT1*2; ++i){
      const char* title = NULL;
      if(i < NumOfSegT1){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT1; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     //			     10000, 0, 400000,
     			     2000, 0, 3000000,
			     "TDC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// Hit parttern
    Int_t target_id = getUniqueID(kT1, 0, kHitPat, 0);
    top_dir->Add(createTH1(++target_id, "T1_HitPat", // 1 origin
			   NumOfSegT1, -0.5, NumOfSegT1+0.5,
			   "Segment", ""));
  }
  { ///// Multiplicity
    Int_t target_id = getUniqueID(kT1, 0, kMulti, 0);
    top_dir->Add(createTH1(++target_id, "T1_Multi", // 1 origin
			   NumOfSegT1+1, -0.5, NumOfSegT1+1.5,
			   "Multiplicity", ""));
  }
  return top_dir;
}

//_____________________________________________________________________________
TList*
HistMaker::createT2(Bool_t flag_ps)
{
  const std::string strDet = "T2";
  name_created_detectors_.push_back(strDet);
  if(flag_ps) name_ps_files_.push_back(strDet);
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  { ///// ADC
    TString strSubDir  = "ADC";
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT2, 0, kADC, 0);
    for(Int_t i = 0; i<NumOfSegT2*2; ++i){
      const char* title = NULL;
      if(i < NumOfSegT2){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT2; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     0x1000, 0, 0x1000,
			     "ADC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// ADC w/TDC
    TString strSubDir  = CONV_STRING(kADCwTDC);
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT2, 0, kADCwTDC, 0);
    for( Int_t i=0; i<NumOfSegT2*2; ++i ){
      const char* title = NULL;
      if( i<NumOfSegT2 ){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT2; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     0x1000, 0, 0x1000,
			     "ADC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// TDC
    TString strSubDir  = CONV_STRING(kTDC);
    const char* nameSubDir = strSubDir.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    Int_t target_id = getUniqueID(kT2, 0, kTDC, 0);
    for(Int_t i = 0; i<NumOfSegT2*2; ++i){
      const char* title = NULL;
      if(i < NumOfSegT2){
	Int_t seg = i+1; // 1 origin
	title = Form("%s_%s_%dU", nameDetector, nameSubDir, seg);
      }else{
	Int_t seg = i+1-NumOfSegT2; // 1 origin
	title = Form("%s_%s_%dD", nameDetector, nameSubDir, seg);
      }
      sub_dir->Add(createTH1(target_id + i+1, title, // 1 origin
			     //			     10000, 0, 400000,
     			     2000, 0, 3000000,
			     "TDC [ch]", ""));
    }
    top_dir->Add(sub_dir);
  }
  { ///// Hit parttern
    Int_t target_id = getUniqueID(kT2, 0, kHitPat, 0);
    top_dir->Add(createTH1(++target_id, "T2_HitPat", // 1 origin
			   NumOfSegT2, -0.5, NumOfSegT2+0.5,
			   "Segment", ""));
  }
  { ///// Multiplicity
    Int_t target_id = getUniqueID(kT2, 0, kMulti, 0);
    top_dir->Add(createTH1(++target_id, "T2_Multi", // 1 origin
			   NumOfSegT2+1, -0.5, NumOfSegT2+1.5,
			   "Multiplicity", ""));
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
  { // TDC
    Int_t target_id = getUniqueID(kBcOutTracking, 0, 0, 0);
    TString title = Form("%s_X", nameDetector);
    top_dir->Add(createTH1(++target_id, title,
			   1000, 0, 1000,
			   "TDC [ch]", ""));
  }
  return top_dir;
}

// -------------------------------------------------------------------------
// createT98Hist
// -------------------------------------------------------------------------
TList* HistMaker::createT98Hist(DetectorType kDET, std::string strDet, const int nch, const int nbin, bool flag_ps)
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
  int kGRAMS40=40;
  int target_id ;  
  for(int ich=0; ich<nch;ich++){  
    target_id = getUniqueID(kAna, 0, 12, ich+1);  
    top_dir->Add(createTH1(target_id+1, Form("Waveform_Ch%d",ich), nbin, -0.5, (double)nbin-0.5));
    if(ich!=0){
      target_id = getUniqueID(kAna, 0, 12, nch+ich);  
      top_dir->Add(createTH1(target_id+1, Form("ChargeSum_Ch%d",ich), 1000000, -0.5, 1000000-0.5));
      target_id = getUniqueID(kAna, 0, 12, nch*2+ich);  
      top_dir->Add(createTH1(target_id+1, Form("FADCPeak_Ch%d",ich), 30000, -0.5, 30000-0.5));
      target_id = getUniqueID(kAna, 0, 12, nch*3+ich);  
      top_dir->Add(createTH1(target_id+1, Form("ChargeSum_Ch%d",ich), 1000000, -0.5, 1000000-0.5));
    }
  }
  target_id = getUniqueID(kAna, 0, 12, nch*4);  
  top_dir->Add(createTH1(target_id+1, "TIGArCounter", 2000, -0.5, 2000-0.5));
  target_id = getUniqueID(kAna, 0, 12, nch*4+1);  
  top_dir->Add(createTH1(target_id+1, "TIGArCounter_EventCounter", 100, -0.5, 100-0.5));
  
  // Return the TList pointer which is added into TGFileBrowser
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
// createQDC
// -------------------------------------------------------------------------
TList* HistMaker::createSDD(DetectorType kDET, std::string strDet, int nports, int nunits, bool flag_ps)
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
  int nsdds=nports*nunits;
  int nbins=512; double xmin=-0.5, xmax=4095.5, xmaxt=8191.5;
  list.push_back( HistMakerInfo(kADC,    "ADC",    "ADC [ch]", nbins,xmin,xmax) );
  list.push_back( HistMakerInfo(kADCwTDC,"ADCwTDC",  "ADC [ch]", nbins, xmin,xmax) );
  list.push_back( HistMakerInfo(kTDC,    "TDC",       "TDC [ch]", nbins, xmin, xmaxt) );
  list.push_back( HistMakerInfo(kTDC2D,  "Trailing",  "Trainling [ch]", nbins, xmin, xmaxt) );
  list.push_back( HistMakerInfo(kTOT,    "TOT",      "TOT [ch]", 1024,0.,1024.) );
  list.push_back( HistMakerInfo(kResetL, "Reset",       "TDC [ch]", nbins, xmin, xmaxt) );
  list.push_back( HistMakerInfo(kResetT,  "ResetTrailing",  "Trainling [ch]", nbins, xmin, xmaxt) );
  list.push_back( HistMakerInfo(kHitPat, "HitPat",  "Segment", nsdds+1 ,-0.5,nsdds-0.5) );
  list.push_back( HistMakerInfo(kMulti,  "Multi",   "Multiplicity", nsdds+1 ,-0.5,nsdds+0.5) );
  int nsegments=8;
  const TString unitname[4]={"A","B","C","D"};
  for(size_t i=0; i<2;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    
    for(int iport=0; iport<nports; ++iport){
      for(int iunit=0; iunit<nunits; ++iunit){
	int index=iunit+nunits*iport;
	int target_id = getUniqueID(kDET, index, list[i].type, 0);
	for(int j = 0; j<nsegments; ++j){
	  const char* title = NULL;
	  title = Form("%s_%s_%d_%s_%d", nameDetector, nameSubDir, iport,unitname[iunit].Data(),j);
	  sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
				 list[i].nbinsx, list[i].xmin, list[i].xmax,
				 list[i].xtitle, list[i].ytitle));
	}
      }
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=2; i<7;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int iport=0; iport<nports; ++iport){
      for(int iunit=0; iunit<nunits; ++iunit){
	int index=iunit+nunits*iport;
	int target_id = getUniqueID(kDET, index, list[i].type, 0);
	const char* title = NULL;
	title = Form("%s_%s_%d_%s", nameDetector, nameSubDir, iport,unitname[iunit].Data());
	sub_dir->Add(createTH1(target_id +1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      }
    }
    top_dir->Add(sub_dir);
  }  
  for(size_t i=7; i<9;i++){
    const char* nameSubDir = list[i].type_name.Data();
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    const char* title = NULL;
    title = Form("%s_%s", nameDetector, nameSubDir);
    top_dir->Add(createTH1(target_id , title, // 1 origin
			   list[i].nbinsx, list[i].xmin, list[i].xmax,
			   list[i].xtitle, list[i].ytitle));
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
// createPbF2
// -------------------------------------------------------------------------
TList* HistMaker::createPbF2(DetectorType kDET, std::string strDet, const int nsegments, bool flag_ps)
{
  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  std::vector<HistMakerInfo> list;;
  int amax=512;
  list.push_back( HistMakerInfo(kADC,      "QDC",      "QDC [ch]", amax,0.,double(amax)) );
  list.push_back( HistMakerInfo(kADCwTDC,  "QDCwTDC",  "QDC [ch]", amax,0.,double(amax)) );
  list.push_back( HistMakerInfo(kADC2D,    "QDC(NtrigT0)",  "QDC [ch]", amax,0.,double(amax)) );
  list.push_back( HistMakerInfo(kFADC,    "QDCwT(NtrigT0)",  "QDC [ch]", amax,0.,double(amax)) );
  
  size_t size=list.size();
  for(size_t i=0; i<size;i++){
    // Declaration of the sub-directory
    const char* nameSubDir = list[i].type_name.Data();
    // Make histogram and add it
    int target_id = getUniqueID(kDET, 0, list[i].type, 0);
    for(int j = 0; j<nsegments; ++j){
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

  // double nbins=2000;
  double nbins=4096;
  
  double xmax=2000;
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
                             nbins/10,0,xmax,100,0,ntot,
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

#if 1
// -------------------------------------------------------------------------
// createCDC
// -------------------------------------------------------------------------
TList* HistMaker::createCDC(bool flag_ps)
{
  DetectorType kDET = kCDC2;
  // Determine the detector name
  std::string strDet = CONV_STRING(kCDC2);
  // name list of crearted detector
  // Declaration of the directory
  // Just type conversion from std::string to char*
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  int nlayers = 15;
  int nwires  = 200; // temporaly
  const int NumOfCDCWiresInLayer[15] = {72,72,72,90,90,100,100,120,120,150,150,160,160,180,180};

  std::vector<HistMakerInfo> list;;
  list.push_back( HistMakerInfo(kTDC,      CONV_STRING(kTDC),    "TDC [ch]", 2000,0.,2000.) );
  list.push_back( HistMakerInfo(kTDC2D,    CONV_STRING(kTDC2D),  "Trailing [ch]", 2000.,0.,2000.) );
  list.push_back( HistMakerInfo(kTOT,      CONV_STRING(kTOT),    "TOT [ch]", 500,0.,500.) );
  list.push_back( HistMakerInfo(kHitPat,   CONV_STRING(kHitPat), "Wire", nwires+1 ,-0.5,nwires+0.5) );
  list.push_back( HistMakerInfo(kMulti,    CONV_STRING(kMulti),  "Multiplicity", nwires+1 ,-0.5,nwires+0.5) );
  list.push_back( HistMakerInfo(kMulti2D,    CONV_STRING(kMulti2D),  "Multiplicity2", nwires+1 ,-0.5,nwires+0.5) );

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
      if( i<3 )
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       list[i].nbinsx, list[i].xmin, list[i].xmax,
			       list[i].xtitle, list[i].ytitle));
      else
	sub_dir->Add(createTH1(target_id + j+1, title, // 1 origin
			       NumOfCDCWiresInLayer[j]+1, -0.5, NumOfCDCWiresInLayer[j]+0.5,
			       list[i].xtitle, list[i].ytitle));
    }
  }  

  for(int j = 0; j<nlayers; ++j){
    const char* nameSubDir = "TDCvsTOT";
    TList *sub_dir = new TList;
    sub_dir->SetName(nameSubDir);
    top_dir->Add(sub_dir);
 
    int target_id = getUniqueID(kDET, 0, kADC2D, 0);
    const char* title = Form("%s_%s_%d", nameDetector, nameSubDir, j);
    sub_dir->Add(createTH2(target_id + j+1, title,
			   200,0,2000,50,0,500,
			   "TDC [ch]", "TOT [ch]") );
  }


  { 
    const char* nameSubDir = "WireCorr";
    // Make histogram and add it
    int layer1[9]={0,3,5,7,9,11,13,7,7};
    int layer2[9]={1,4,6,8,10,12,14,1,14};
    for(int i=0;i<9;i++){
      int target_id = getUniqueID(kDET, 0, kWireCorr, i);
      int nwires1=NumOfCDCWiresInLayer[layer1[i]];
      int nwires2=NumOfCDCWiresInLayer[layer2[i]];
      const char* title = NULL;
      title = Form("%s_%s_%d", nameDetector, nameSubDir,i);
      top_dir->Add(createTH2(target_id, title, // 1 origin
			     nwires1+1,-0.5,nwires1+0.5,nwires2+1,-0.5,nwires2+0.5,
			     Form("Wire in layer%d",layer1[i]), 
			     Form("Wire in layer%d",layer2[i])) );
    }
  }

  // Return the TList pointer which is added into TGFileBrowser
  return top_dir;
}
#endif

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

  // DAQ infomation --------------------------------------------------
  {
    // Event builder infomation
    int target_id = getUniqueID(kDAQ, kEB, kHitPat, 0);
    top_dir->Add(createTH1(target_id + 1, "Data size EB", // 1 origin
			   5000, 0, 5000,
			   "Data size [words]", ""));

    // Node information
    target_id = getUniqueID(kDAQ, kVME, kHitPat2D, 0);
    top_dir->Add(createTH2(target_id + 1, "Data size VME nodes", // 1 origin
			   15, 0, 15,
			   100, 0, 1200,
			   "VME node ID", "Data size [words]"));

    target_id = getUniqueID(kDAQ, kCLite, kHitPat2D, 0);
    top_dir->Add(createTH2(target_id + 1, "Data size CLite nodes", // 1 origin
			   15, 0, 15,
			   200, 0, 400,
			   "CLite node ID", "Data size [words]"));

    target_id = getUniqueID(kDAQ, kEASIROC, kHitPat2D, 0);
    top_dir->Add(createTH2(target_id + 1, "Data size EASIROC nodes", // 1 origin
			   20, 0, 20,
			   50, 0, 100,
			   "EASIROC node ID", "Data size [words]"));

    target_id = getUniqueID(kDAQ, kCAMAC, kHitPat2D, 0);
    top_dir->Add(createTH2(target_id + 1, "Data size CAMAC nodes", // 1 origin
			   3, 0, 3,
			   100, 0, 200,
			   "CAMAC node ID", "Data size [words]"));

    target_id = getUniqueID(kDAQ, kMiscNode, kHitPat2D, 0);
    top_dir->Add(createTH2(target_id + 1, "Data size Misc nodes", // 1 origin
			   5, 0, 5,
			   100, 0, 200,
			   "Misc node ID", "Data size [words]"));
  }
  return top_dir;
}


TList* HistMaker::createAnaTOF(int type,bool flag_ps)
{
  // Determine the detector name
  std::string strDet = CONV_STRING(kAna);

  // Declaration of the directory
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  // TOF --------------------------------------------------
  {
    int target_id = getUniqueID(kAna, 0, 1, 0);
    TString name[5]={"BHT","T0","DEF","VETO","BTC"};
    for(int i=0;i<5;i++){
      TString hname=Form("TOF %s-T1",name[i].Data());
      top_dir->Add(createTH1(target_id + 1 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if ACHit",name[i].Data());
      top_dir->Add(createTH1(target_id + 2 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if TOF Kaon",name[i].Data());
      top_dir->Add(createTH1(target_id + 3 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if TOF Pion",name[i].Data());
      top_dir->Add(createTH1(target_id + 4 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if TOF Proton",name[i].Data());
      top_dir->Add(createTH1(target_id + 5 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if TOF Deuteron",name[i].Data());
      top_dir->Add(createTH1(target_id + 6 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if B_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 7 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if Pi_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 8 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if K_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 9 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if P_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 10 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if Pion in B_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 11 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if Kaon in B_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 12 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
      hname=Form("TOF %s-T1 if Proton in B_trg",name[i].Data());
      top_dir->Add(createTH1(target_id + 13 + i*20, hname, // 1 origin
			     2000, -100, 100,"TOF(ns)", ""));   
    }
  }
  // AC ------------------
  {
    int target_id = getUniqueID(kAna, 0, 10, 0);
    TString name[4]={"ACsum_All","ACsum_wTDC","ACsum_Kaon","ACsum_Pion"};
    for(int i=0;i<4;i++){
      top_dir->Add(createTH1(target_id + 1 + i, name[i], // 1 origin
			     100,0,4000,
			     "ACsum",""));   
    }
    int target_id2 = getUniqueID(kAna, 0, 11, 0);
    for(int i=0;i<4;i++){
      top_dir->Add(createTH1(target_id2 + 1 + i, name[i]+"_ifBtrg", // 1 origin
			     100,0,4000,
			     "ACsum",""));   
    }
  }
  return top_dir;
}


TList* HistMaker::createAnaTrack(int type,bool flag_ps)
{
  // Determine the detector name
  std::string strDet = CONV_STRING(kAna);

  // Declaration of the directory
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);

  // BLC2 - T0 corr --------------------------------------------------
  {
    // Event builder infomation
    int target_id1 = getUniqueID(kAna, 0, 2, 0);
    int target_id2 = getUniqueID(kAna, 0, 3, 0);
    int target_id3 = getUniqueID(kAna, 0, 4, 0);
    int target_id4 = getUniqueID(kAna, 0, 5, 0);
    for(int i=0;i<5;i++){
      top_dir->Add(createTH2(target_id1 + 1 + i, Form("BLC2a with T0 seg%d", i+1), // 1 origin
			     100,-100,100,100,-100,100,
			     "X", "Y"));   
      top_dir->Add(createTH2(target_id2 + 1 + i, Form("BLC2b with T0 seg%d", i+1), // 1 origin
			     100,-100,100,100,-100,100,
			     "X", "Y"));   
      top_dir->Add(createTH2(target_id3 + 1 + i, Form("BPC1 with DEF seg%d", i+1), // 1 origin
			     100,-100,100,100,-100,100,
			     "X", "Y"));   
      top_dir->Add(createTH2(target_id4 + 1 + i, Form("BPC2 with DEF seg%d", i+1), // 1 origin
			     100,-100,100,100,-100,100,
			     "X", "Y"));   
    }

    int target_id = getUniqueID(kAna, 0, 40, 0);
    top_dir->Add(createTH2(target_id + 1, "X corr BLC1a-BLC1b", // 1 origin
                           100,-100,100,100,-100,100,
                           "BLC1a X", "BLC1b X"));   
    top_dir->Add(createTH2(target_id + 2, "Y corr BLC1a-BLC1b", // 1 origin
                           100,-100,100,100,-100,100,
                           "BLC1a Y", "BLC1b Y"));   
    top_dir->Add(createTH2(target_id + 3, "A corr BLC1a-BLC1b", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BLC1a A", "BLC1b A"));   
    top_dir->Add(createTH2(target_id + 4, "B corr BLC1a-BLC1b", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BLC1a B", "BLC1b B"));   

    top_dir->Add(createTH2(target_id + 11, "X corr BLC2a-BLC2b", // 1 origin
                           100,-100,100,100,-100,100,
                           "BLC2a X", "BLC2b X"));   
    top_dir->Add(createTH2(target_id + 12, "Y corr BLC2a-BLC2b", // 1 origin
                           100,-100,100,100,-100,100,
                           "BLC2a Y", "BLC2b Y"));   
    top_dir->Add(createTH2(target_id + 13, "A corr BLC2a-BLC2b", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BLC2a A", "BLC2b A"));   
    top_dir->Add(createTH2(target_id + 14, "B corr BLC2a-BLC2b", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BLC2a B", "BLC2b B"));   

    top_dir->Add(createTH2(target_id + 21, "X corr BPC1-BPC2", // 1 origin
                           100,-100,100,100,-100,100,
                           "BPC1 X", "BPC2 X"));   
    top_dir->Add(createTH2(target_id + 22, "Y corr BPC1-BPC2", // 1 origin
                           100,-100,100,100,-100,100,
                           "BPC1 Y", "BPC2 Y"));   
    top_dir->Add(createTH2(target_id + 23, "A corr BPC1-BPC2", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BPC1 A", "BPC2 A"));   
    top_dir->Add(createTH2(target_id + 24, "B corr BPC1-BPC2", // 1 origin
                           100,-0.1,0.1,100,-0.1,0.1,
                           "BPC1 B", "BPC2 B"));   
  }


  return top_dir;
}

TList* HistMaker::createAnaDeuteron(int type,bool flag_ps)
{
  // Determine the detector name
  std::string strDet = CONV_STRING(kAna);

  // Declaration of the directory
  const char* nameDetector = strDet.c_str();
  TList *top_dir = new TList;
  top_dir->SetName(nameDetector);
  // deuteron -----------
  {
    TString trig_str[5]={"_all","_ifBtrg","_ifKtrg","_ifPitrg","_ifPtrg"};
    int target_id = getUniqueID(kAna, 0, 21, 0);
    for(int i=0;i<5;i++){
      top_dir->Add(createTH2(target_id + 1+i, "TOF_T1BHT_T1VETO"+trig_str[i],
			     200, -50, 50, 80, -10, 30,
			     "T1BHT TOF(ns)", "T1VETO TOF(ns)"));   
    }
    target_id = getUniqueID(kAna, 0, 22, 0);
    for(int i=0;i<5;i++){
      top_dir->Add(createTH2(target_id + 1+i, "TOF_T1BHT_T1BTC"+trig_str[i],
			     200, -50, 50, 80, -10, 30,
			     "T1BHT TOF(ns)", "T1BTC TOF(ns)"));   
    }
    TString tof_str[5]={"_All","_ifTOFK","_ifTOFPi","_ifTOFP","_ifTOFD"};
    target_id = getUniqueID(kAna, 0, 23, 0);
    for(int i=0;i<5;i++){
      top_dir->Add(createTH2(target_id + 1+i, "TOFVETO_T1ADCsum"+tof_str[i],
			     80,-10,30,100,0,4000,
			     "TOF(ns)", "T1 ADCsum"));   
    }
    target_id = getUniqueID(kAna, 0, 24, 0);
    for(int i=0;i<5;i++){
      top_dir->Add(createTH2(target_id + 1+i, "TOFBTC_T1ADCsum"+tof_str[i],
			     80,-10,30,100,0,4000,
			     "TOF(ns)", "T1 ADCsum"));   
    }    
  }
  return top_dir;
}
