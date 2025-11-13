// -*- C++ -*-

// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include "defines.hh"
#include "user_analyzer.hh"
#include "UnpackerManager.hh"
#include "Unpacker.hh"
#include "DAQNode.hh"
#include "ConfMan.hh"

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    ConfMan& gConfMan = ConfMan::getInstance();
    UnpackerManager& gUnpacker = GUnpacker::get_instance();
    //std::vector<std::string> target = { "hul01", "hul03" };
  // std::vector<std::string> target = { "vme01" };
   std::vector<std::string> target = { "hul-hrtdc-125" };
  }

//______________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  gConfMan.initialize(argv);

  // for(auto&& c: gUnpacker.get_root()->get_child_list()){
  //   c.second->set_dump_mode(defines::k_hex);
  // }

  for( int i=0, n=target.size(); i<n; ++i ){
    int node_id = gUnpacker.get_fe_id( target.at(i) );
    Unpacker *node = gUnpacker.get_root()->get_child(node_id);
    if( !node ) continue;
    node->set_dump_mode(defines::k_hex);
  }

  return 0;
}

//______________________________________________________________________________
int
process_end( void )
{
  return 0;
}

//______________________________________________________________________________
int
process_event( void )
{
  // std::cout<<gUnpacker.get_run_number()<<std::endl;

  // int node_id = gUnpacker.get_fe_id( target.at(0) );
  
//  gUnpacker.dump_data_fe(node_id, 0xff010000);
  /*
    int node_id = gUnpacker.get_fe_id("k18eb"); // Event builder
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_magic) << std::dec << std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_data_size) << std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_event_number) << std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_run_number) << std::hex <<std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_node_id) << std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_node_type) << std::dec <<std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_number_of_blocks) << std::dec <<std::endl;
    std::cout << g_unpacker.get_node_header(node_id, DAQNode::k_unix_time) << std::dec <<std::endl;

    std::cout << "EEE" << std::endl;
    std::cout << g_unpacker.get_root()->get_run_number() << std::endl;
  */
  return 0;
}

}
