/*
  ConfMan.cc
*/

#include "ConfMan.hh"
#include "UnpackerManager.hh"
#include "Unpacker.hh"
#include "filesystem_util.hh"

#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
//#include "DCDriftParamMan.hh"
#include "UserParamMan.hh"
#include "TpcPadHelper.hh"
//#include "MatrixParamMan.hh"
//#include "MsTParamMan.hh"
//#include "GeAdcCalibMan.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

static const std::string MyName = "ConfMan::";

// Constructor -------------------------------------------------------------
ConfMan::ConfMan()
{
  flag_.reset();
}

// Destructor --------------------------------------------------------------
ConfMan::~ConfMan()
{

}

// initialize --------------------------------------------------------------
void ConfMan::initialize(const std::vector<std::string>& argv)
{
  using namespace hddaq;
  using namespace hddaq::unpacker;
  static const std::string MyFunc = "initialize ";

  std::cout << "#D " << MyName << MyFunc
	    << "argument list" << std::endl;
  std::copy(argv.begin(), argv.end(),
	    std::ostream_iterator<std::string>(std::cout, "\n"));
  std::cout << std::endl;

  int nArg = argv.size();
  if (sizeArgumentList > nArg + 1 )
    {
      std::cerr << "#E " << MyName << MyFunc << std::endl
		<< " Usage: "
		<< basename( argv[kProcess] )
		<< " [config file] [input stream]"
		<< std::endl;
      std::exit(1);
    }


  const std::string& confFile(argv[kConfPath]);
  const std::string& dataSrc(argv[kStreamPath]);
  //  const std::string& servPort(argv[kServPort]);

  std::cout << " config file = " << confFile << std::endl;
  std::string dir = dirname(confFile);
  dir += "/";
  std::cout << " dir = " << dir << std::endl;
  std::ifstream conf(confFile.c_str());
  while (conf.good())
    {
      std::string l;
      std::getline(conf, l);
      std::istringstream iss(l);
      std::istream_iterator<std::string> issBegin(iss);
      std::istream_iterator<std::string> issEnd;
      std::vector<std::string> param(issBegin, issEnd);
      if (param.empty())
	continue;
      if (param[0].empty())
	continue;
      if (param.size()==2)
	{
	  const std::string& name = param[0];
	  std::string value = param[1];
	  std::cout << " key = " << name
		    << " value = " << value << std::endl;
          if (name.find("#")!=std::string::npos)
            continue;
          if (name.find("CDCXT")!=std::string::npos)
            continue;
          if (value.find("_list")!=std::string::npos)
            continue;
	  if (value.find("/")!=0)
	    value = realpath(dir + value);
	  
	  name_file_[name] = value;
	}
    }

  for(const auto& n : name_file_){
    std::cout << n.first << "\t" << n.second << std::endl;
  }

  // std::cout << name_file_["UNPACK:"] << std::endl;
  // std::cout << name_file_["DIGIT:"] << std::endl;
  // std::cout << name_file_["CMAP:"] << std::endl;

  // initialize unpacker system
  UnpackerManager& g_unpacker = GUnpacker::get_instance();
  g_unpacker.set_config_file(name_file_["UNPACK:"],
			     name_file_["DIGIT:"],
			     name_file_["CMAP:"]);
  g_unpacker.set_istream(dataSrc);

  // Initialize of ConfMan and Unpacker were done
  flag_.set(kIsGood);
}
