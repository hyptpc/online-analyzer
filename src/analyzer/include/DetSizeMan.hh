// -*- C++ -*-

#ifndef DET_SIZE_MAN_HH
#define DET_SIZE_MAN_HH

#include <string>
#include <map>
#include <vector>

#include <TString.h>

#include "ThreeVector.hh"

//_____________________________________________________________________________
class DetSizeMan
{
public:
  static TString    ClassName();
  static DetSizeMan& GetInstance();
  ~DetSizeMan();

private:
  DetSizeMan();
  DetSizeMan(const DetSizeMan& );
  DetSizeMan& operator =(const DetSizeMan&);

private:
  typedef std::vector<double>          ParamArray;
  typedef std::map<TString, ParamArray> ParamMap;
  typedef ParamMap::const_iterator       PIterator;
  bool   m_is_ready;
  TString m_file_name;
  ParamMap m_param_map;

public:
  bool        Initialize();
  bool        Initialize(const TString& filename);
  bool        IsReady() const { return m_is_ready; }
  void        SetFileName(const TString &file_name);
  double      Get(const TString& key, int i=0) const;
  ThreeVector GetSize(const TString& key) const;
  void          Print() const;
};

//_____________________________________________________________________________
inline TString
DetSizeMan::ClassName()
{
  static TString s_name("DetSizeMan");
  return s_name;
}

//_____________________________________________________________________________
inline DetSizeMan&
DetSizeMan::GetInstance()
{
  static DetSizeMan s_instance;
  return s_instance;
}

#endif
