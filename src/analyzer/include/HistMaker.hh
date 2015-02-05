#ifndef HISTMAKER_H
#define HISTMAKER_H 1

#include<vector>
#include<map>
#include<bitset>
#include<string>

#include<TROOT.h>

enum DetectorType{
  kDetectorZero,
  // Detector unique ID in the beam line
  kBH1, kBFT, kBC3, kBC4, kBMW, kBH2, kBAC_SAC,
  // Detector unique ID in the SKS system
  kSDC2, kHDC, kSP0, kSDC3, kSDC4, kTOF, kTOFMT, kSFV_SAC3, kLC,
  // Detector unique ID in Hyperball
  kGe, kPWO,
  // Others
  kTriggerFlag, kCorrelation, kMisc,
  sizeDetectorType,
  factorDetectorType = 10000000
};

enum SubDetectorType{
  kSubDetectorZero,
  // Detector unique sub ID in Counters
  kSP0_L1, kSP0_L2, kSP0_L3, kSP0_L4,
  kSP0_L5, kSP0_L6, kSP0_L7, kSP0_L8,
  // Detector unique sub ID in Hyperball
  kPWO_B, kPWO_E, kPWO_C, kPWO_L,
  sizeSubDetectorType,
  factorSubDetectorType = 100000
};

enum DataType{
  kDataTypeZero,
  // Usual data type
  kADC, kTDC, kHitPat, k2DPlot, kMulti, 
  // Extra data type for Ge detector
  kCRM, kTFA, kPUR, kRST, 
  sizeDataType,
  factorDataType = 1000
};

std::string getStr_FromEnum(const char* c);

class TH1;
class TH2;
class TList;

class HistMaker{
  // Declaration of the private parameters ---------------------------------
  // histogram unique and sequential ID 
  int current_hist_id_; 
  std::map<int, int>         idmap_seq_from_unique_;
  std::map<int, int>         idmap_unique_from_seq_;
  std::map<std::string, int> idmap_seq_from_name_;

  // data set of string which means the created detector name
  std::vector<std::string>      name_created_detectors_;
  std::vector<std::string>      name_ps_files_;

  
public:
  // Public functions ------------------------------------------------------
  ~HistMaker();
  static HistMaker& getInstance();
  static void getListOfDetectors(std::vector<std::string>& vec);
  static void getListOfPsFiles(std::vector<std::string>& vec);

  static int  getUniqueID(int detector_type, int subdetector_type,
			 int data_type, int channel=1 );
  static int  getUniqueID(int sequential_id);
  
  static int  getSequentialID(int detector_type, int subdetector_type,
			     int data_type, int channel=1 );
  static int  getSequentialID(int unique_id);
  static int  getSequentialID(const char* name_hist);

  int  setHistPtr(std::vector<TH1*>& vec);

  TH1* createTH1(int unique_id, const char* title,
		 int nbinx, int xmin, int xmax,
		 const char* xtitle="", const char* ytitle=""
		 );
  
  TH2* createTH2(int unique_id, const char* title,
  		 int nbinx, int xmin, int xmax,
  		 int nbiny, int ymin, int ymax,
  		 const char* xtitle="", const char* ytitle=""
  		 );

  TList* createBH1(bool flag_ps=true);
  TList* createBFT(bool flag_ps=true);
  TList* createBC3(bool flag_ps=true);
  TList* createBC4(bool flag_ps=true);
  TList* createBMW(bool flag_ps=true);
  TList* createBH2(bool flag_ps=true);
  TList* createBAC_SAC(bool flag_ps=true);
  TList* createSDC2(bool flag_ps=true);
  TList* createHDC(bool flag_ps=true);
  TList* createSP0(bool flag_ps=true);
  TList* createSDC3(bool flag_ps=true);
  TList* createSDC4(bool flag_ps=true);
  TList* createTOF(bool flag_ps=true);
  TList* createTOFMT(bool flag_ps=true);
  TList* createSFV_SAC3(bool flag_ps=true);
  TList* createLC(bool flag_ps=true);
  TList* createGe(bool flag_ps=true);
  TList* createPWO(bool flag_ps=true);
  TList* createTriggerFlag(bool flag_ps=true);
  TList* createCorrelation(bool flag_ps=true);

private:
  HistMaker();
  HistMaker(const HistMaker& object);
  HistMaker& operator=(const HistMaker& object);

  ClassDef(HistMaker, 0)
};

// getInstance -----------------------------------------------------------    
inline HistMaker& HistMaker::getInstance()
{
  static HistMaker object;
  return object;
}

// getUniqueID -----------------------------------------------------------    
inline int HistMaker::getUniqueID(int detector_type, int subdetector_type,
				  int data_type, int channel)
{
  int unique_id = 0;
  unique_id += detector_type*factorDetectorType;
  unique_id += subdetector_type*factorSubDetectorType;
  unique_id += data_type*factorDataType;
  unique_id += channel;
  return unique_id;
}

// getUniqueID -----------------------------------------------------------    
inline int HistMaker::getUniqueID(int sequential_id)
{
  HistMaker& g= HistMaker::getInstance();
  return g.idmap_unique_from_seq_[sequential_id];
}

// getSequentialID -------------------------------------------------------    
inline int HistMaker::getSequentialID(int detector_type, int subdetector_type,
				      int data_type, int channel)
{
  HistMaker& g= HistMaker::getInstance();
  int unique_id = g.getUniqueID(detector_type, subdetector_type, data_type, channel);
  return g.idmap_seq_from_unique_[unique_id];
}

// getSequentialID -------------------------------------------------------    
inline int HistMaker::getSequentialID(int unique_id)
{
  HistMaker& g= HistMaker::getInstance();
  return g.idmap_seq_from_unique_[unique_id];
}

// getSequentialID -------------------------------------------------------    
inline int HistMaker::getSequentialID(const char* name)
{
  HistMaker& g= HistMaker::getInstance();
  std::string nameStr = name;
  return g.idmap_seq_from_name_[nameStr];
}


#endif
