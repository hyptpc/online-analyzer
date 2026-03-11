// -*- C++ -*-

#include "DetSizeMan.hh"
#include "ConfMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

#include "FuncName.hh"

//_____________________________________________________________________________
DetSizeMan::DetSizeMan()
  : m_is_ready(false),
    m_file_name(),
    m_param_map()
{
}

//_____________________________________________________________________________
DetSizeMan::~DetSizeMan()
{
}

//______________________________________________________________________________
void
ConfMan::initializeDetSizeMan()
{
  if(name_file_["DSIZE:"] != ""){
    DetSizeMan& gSize = DetSizeMan::GetInstance();
    gSize.SetFileName(name_file_["DSIZE:"]);
    flag_[kIsGood] = gSize.Initialize();
  }else{
    std::cout << "#E ConfMan::"
	      << " File path does not exist in " << name_file_["DSIZE:"] 
	      << std::endl;
    flag_.reset(kIsGood);
  }
}

//_____________________________________________________________________________
bool
DetSizeMan::Initialize()
{
  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    std::cerr << "#E " << FUNC_NAME << " "
           << "No such parameter file : " << m_file_name << std::endl;
    return false;
  }

  std::string line;
  while(ifs.good() && std::getline(ifs, line)){
    if(line[0]=='#') continue;
    std::istringstream input_line(line);

    TString first_param;
    input_line >> first_param;

    TString key = first_param;
    ParamArray param_array;
    double   param;
    while(input_line >> param){
      param_array.push_back(param);
    }
    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
bool
DetSizeMan::Initialize(const TString& filename)
{
  m_file_name = filename;
  return Initialize();
}

//_____________________________________________________________________________
void
DetSizeMan::SetFileName(const TString& file_name)
{
  m_file_name = file_name;
}

//_____________________________________________________________________________
double
DetSizeMan::Get(const TString& key, int i) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if(itr==m_param_map.end() ||
     i+1 > (int)itr->second.size()){
    Print();
    TString msg(FUNC_NAME+" No such key : "+key);
    msg += "(i=" + std::to_string(i) + ")";
    throw std::invalid_argument(msg);
  }

  return itr->second.at(i);
}

//_____________________________________________________________________________
ThreeVector
DetSizeMan::GetSize(const TString& key) const
{
  return ThreeVector(Get(key, 0),
                       Get(key, 1),
                       Get(key, 2));
}

//_____________________________________________________________________________
void
DetSizeMan::Print() const
{
  std::cout << FUNC_NAME << std::endl;

  const int w = 20;
  PIterator itr, end=m_param_map.end();
  for(itr=m_param_map.begin(); itr!=end; ++itr){
    std::cout << " key = " << std::setw(w) << std::left
           << itr->first << itr->second.size() << " : ";
    for(int i=0, n=itr->second.size(); i<n; ++i){
      std::cout << std::setw(5) << std::right
             << itr->second.at(i) << " ";
    }
    std::cout << std::endl;
  }
}
