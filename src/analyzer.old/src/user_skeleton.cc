// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "user_analyzer.hh"
#include "UnpackerManager.hh"
#include "ConfMan.hh"

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  //  GUnpacker::get_instance().set_decode_mode(false);

  return 0;
}

//____________________________________________________________________________
int
process_end()
{
  std::cout << "#D skeleton::process_end End" << std::endl;
  return 0;
}

//____________________________________________________________________________
int
process_event()
{
  sleep(1);
  std::cout << "#D skeleton::process_event" << std::endl;

  return 0;
}

}
