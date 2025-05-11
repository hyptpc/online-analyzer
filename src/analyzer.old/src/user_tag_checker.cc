// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <vector>

#include "user_analyzer.hh"
#include "UnpackerManager.hh"
#include "DAQNode.hh"
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
  GUnpacker::get_instance().set_decode_mode(false);

  return 0;
}

//____________________________________________________________________________
int
process_end()
{
  return 0;
}

//____________________________________________________________________________
int
process_event()
{
  static int event_count = 0;
  if(event_count%1000 == 0){
    std::cout << "#D Tag checker is running" << std::endl;
  }

  ++event_count;

  return 0;
}

}
