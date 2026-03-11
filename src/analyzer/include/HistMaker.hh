#ifndef HISTMAKER_HH
#define HISTMAKER_HH

#include<iostream>
#include<vector>
#include<map>
#include<bitset>
#include<string>

#include<TROOT.h>

enum DetectorType {
  kDetectorZero,
  // Detector unique ID in the beam line
#if 1
  kCDC2,
#endif
  kCDC, kCDH, //3
  kBHD, kBHT=kBHD, kT0, kT0new, kAC, kDEF, kE0, //9
  kBLC1a, kBLC1b, kBLC2a, kBLC2b, kBPCmini, kBPC1=kBPCmini, kBPC, kBPC2=kBPC,  //15
  kPbG, kPbF2, kVeto, kFinger,kBTC,kLeak,kRC, kCNCtest,kT98PMT,kT98MPPC,//25
  kSDD, kSDDGate, kSDDReset, //28
  kCVC, kNC, kVFT, //31
  kBAC, kKVC1, kKVC2, kSAC, kBH2, kCOBO, //34 for E72
  kHTOF, kSAC3, kKVC,
  kBVH, kBVH2, kT1, kT2, kSFV,
  kTPC, kPede,
  kQDC, kQDC1, kQDC2, //37
  kTriggerFlag, kDAQ, kCorrelation, kMisc, kEventDisplay,
  kTimeStamp, kDCEff,kAna, 
  kHRTDC1,kHRTDC2,kHRTDC3, kHRTDC4,kHRTDC5,kHRTDC6, 
  kMHTDC,
  kBcOutTracking,
  sizeDetectorType,
  factorDetectorType = 10000000
};
 
enum SubDetectorType {
  kSubDetectorZero,
  // DAQ unique sub ID
  kEB, kTKO, kVME, kCLite, kEASIROC, kCAMAC,
  kMiscNode,
  sizeSubDetectorType,
  factorSubDetectorType = 50000
}; /// max 200

enum DataType{
  kDataTypeZero,
  // Usual data type
  kADC,   kTDC,   kHitPat,   kMulti, //4
  kADC2D, kTDC2D, kHitPat2D, kMulti2D, //8
  kHitPoly, 
  kADCwTDC, kFADC, kTOT, //11
  kDeltaE, kCTime, kDeltaE2D, kCTime2D, //15
  kResetL, kResetT, //17
  kEnergy1, kEnergy2,  kEnergy3, kEnergy4,  kEnergy5, kEnergy6, //23
  kChisqr, kEff,kProf, kWireCorr, kResid, //28
  kThreshold,
  sizeDataType,
  factorDataType = 1000
};//max 50

enum eUorD { kU, kD, kUorD };

struct HistMakerInfo
{
  DataType type;
  TString type_name;
  TString xtitle;
  int nbinsx;
  double xmin;
  double xmax;
  TString ytitle;
  int nbinsy;
  double ymin;
  double ymax;
  //  HistMakerInfo( void ) {}
  HistMakerInfo( DataType t, std::string n,
		 const char x[20], int nx, double xl, double xu)
    : type(t), type_name(n), 
      xtitle(x), nbinsx(nx), xmin(xl), xmax(xu),
      ytitle(""), nbinsy(999), ymin(-999), ymax(999)
  {
  }
};

std::string getStr_FromEnum(const char* c);

class TH1;
class TH2;
class TList;

class HistMaker {
  // Declaration of the private parameters ---------------------------------
  // histogram unique and sequential ID
  int current_hist_id_;
  std::map<int, int>         idmap_seq_from_unique_;
  std::map<int, int>         idmap_unique_from_seq_;
  std::map<std::string, int> idmap_seq_from_name_;

  // data set of string which means the created detector name
  std::vector<std::string>      name_created_detectors_;
  std::vector<std::string>      name_ps_files_;

  // return value of insert
  typedef std::pair<std::map<int, int>::iterator, bool> TypeRetInsert;

public:
  // Public functions ------------------------------------------------------
  virtual ~HistMaker( void );
  static HistMaker& getInstance( void );
  static void getListOfDetectors(std::vector<std::string>& vec);
  static void getListOfPsFiles(std::vector<std::string>& vec);

  static int  getNofHist();

  static int  getUniqueID(int detector_type, int subdetector_type,
			 int data_type, int channel=1 );
  static int  getUniqueID(int sequential_id);

  static int  getSequentialID(int detector_type, int subdetector_type,
			     int data_type, int channel=1 );
  static int  getSequentialID(int unique_id);
  static int  getSequentialID(const char* name_hist);

  int  setHistPtr(std::vector<TH1*>& vec);

  TH1* createTH1(int unique_id, const char* title,
		 int nbinx, double xmin, double xmax,
		 const char* xtitle="", const char* ytitle=""
		 );
  
  TH2* createTH2(int unique_id, const char* title,
  		 int nbinx, double xmin, double xmax,
  		 int nbiny, double ymin, double ymax,
  		 const char* xtitle="", const char* ytitle=""
  		 );
  TH2*   createTH2Poly( Int_t unique_id, const TString& title,
                        Double_t xmin, Double_t xmax,
                        Double_t ymin, Double_t ymax );
  
  TList* createTriggerFlag(bool flag_ps=true);
  TList* createBVH(bool flag_ps=true);
  TList* createT1(bool flag_ps=true);
  TList* createT2(bool flag_ps=true);
  TList* createTPC( Bool_t flag_ps=true);
  TList* createBTOF( Bool_t flag_ps=true);
  TList* createCorrelation(bool flag_ps=true);
  TList* createEventDisplay(bool flag_ps=true);

  TList* createBcOutTracking(bool flag_ps=true);  

  TList* createTimeStamp(bool flag_ps=true);
  TList* createPbF2(DetectorType kDET, std::string strDet, int nsegments, bool flag_ps=true); 
  TList* createHodo(DetectorType kDET, std::string strDet, int nsegments, int nud, int nbins, double xmin, double xmax, int nbins2, double xmin2, double xmax2, bool flag_ps=true); 
  TList* createBHT(DetectorType kDET, std::string strDet, int nsegments, int nud, int nbins2, double xmin2, double xmax2, bool flag_ps=true); 
  TList* createHTOF(DetectorType kDET, std::string strDet, int nsegments, int nud, int nbins, double xmin, double xmax, int nbins2, double xmin2, double xmax2, bool flag_ps=true); 
  TList* createQDC(DetectorType kDET, std::string strDet, int nsegments, int nbins, double xmin, double xmax, bool flag_ps=true); 
  TList* createSDD(DetectorType kDET, std::string strDet, int nports, int nunits, bool flag_ps=true); 
  TList* createMHTDC(DetectorType kDET, std::string strDet, int nsegments, double xmax=1024., bool flag_ps=true);
  TList* createBLDC(DetectorType kDET, std::string strDet, int nlayers, int nwires, 
                    bool WIRE_RAW = false, bool ANA=false,bool flag_ps=true);
  
#if 1
  TList* createCDC(bool flag_ps=true);
#endif
  TList* createDAQ(bool flag_ps=true);

  TList* createAnaTOF(int type=1, bool flag_ps=true);
  TList* createAnaDeuteron(int type=1, bool flag_ps=true);
  TList* createAnaTrack(int type=1, bool flag_ps=true);
  
  TList* createT98Hist(DetectorType kDET, std::string strDet, int nch,  int nbin, bool flag_ps=true);


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

// getNofHist ------------------------------------------------------------
inline int HistMaker::getNofHist()
{
  HistMaker& g= HistMaker::getInstance();
  return g.idmap_unique_from_seq_.size();
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
