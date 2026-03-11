/*
  ConfMan.hh
*/

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include<string>
#include<vector>
#include<map>
#include<bitset>

class ConfMan {
  // Private parameter declarations ----------------------------------------
  // Arguments of main function
  enum ArgumentList {
    kProcess, kConfPath, kStreamPath, kServPort,
    sizeArgumentList
  };
  std::string path_file_[sizeArgumentList];

  // File path list in configuration file
  std::map<std::string, std::string> name_file_;

  // Flags
  enum Flag {
    kIsGood,
    sizeFlag
  };
  std::bitset<sizeFlag> flag_;

public:
  // Public functions ------------------------------------------------------
  ~ConfMan();
  void initialize(const std::vector<std::string>& argv);
  static ConfMan&  getInstance();
  bool isGood();

  void initializeHodoParamMan();
  void initializeHodoPHCMan();

  void initializeDCGeomMan();
  void initializeDCTdcCalibMan();
  void initializeDCDriftParamMan();

  void initializeDetSizeMan();

  void initializeEMCParamMan();
  void initializeMatrixParamMan();
  void initializeMsTParamMan();

  void initializeUserParamMan();
  void initializeTpcPadHelper();

  void initializeGeAdcCalibMan();

  void initializeBLDCWireMapMan();
  void initializeXTMapMan();
  
private:
  ConfMan();
  ConfMan(const ConfMan& object);
  ConfMan& operator =(const ConfMan& object);
};

// GetInstance
inline ConfMan& ConfMan::getInstance()
{
  static ConfMan object;
  return object;
}

// isGood
inline bool ConfMan::isGood()
{
  return flag_[kIsGood];
}

#endif
