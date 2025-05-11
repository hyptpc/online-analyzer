// Author: Tomonori Takahashi

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <vector>

#include "user_analyzer.hh"
#include "Unpacker.hh"
#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include <Rtypes.h>

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;
  namespace
  {
    UnpackerManager& gUnpacker = GUnpacker::get_instance();
    enum eType { kTag, kLock, nType };
    enum eData { kEvent, kSpill, kSerial, kTime, nData };
  }  
//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  ConfMan& gConfMan = ConfMan::getInstance();
  gConfMan.initialize(argv);
  //  GUnpacker::get_instance().set_decode_mode(false);

  //  new THttpServer("http:8080");
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
  static Int_t  event_count = 0;
  ++event_count;
  if(event_count%100!=0) return 0;
  std::cout << "\033[2J" << std::endl;
  int event_number = gUnpacker.get_event_number();
  int run_number = gUnpacker.get_run_number();
  std::cout << std::left  << std::setw(16) << "RUN"
	    << std::right << std::setw(16) <<  run_number <<" : "
	    << std::left  << std::setw(16) << "Event Number"
	    << std::right << std::setw(16) <<  event_number  << " : "
	    << std::endl << std::endl;

  int event_tag[6];
  for(int i=0;i<6;i++) event_tag[i]=-1;
  int index=-1;
  {
    const char* name="VME-RM";
    const int k_device = gUnpacker.get_device_id(name);    
    for(int plane=0;plane<1;plane++){
      //      if(plane==1) continue;
      int spill=-1;
      int evnum=-1;
      for( int ch=0; ch<nType; ++ch ){
	for( int data=0; data<nData; ++data ){
	  int nhit = gUnpacker.get_entries( k_device, plane, 0, ch, data );  
	  if( nhit<=0 ) continue;
	  unsigned int val = gUnpacker.get( k_device, plane, 0, ch, data );  
	  //std::cout<<ch<<" "<<data<<"  "<<val<<std::endl;
	  if( ch==kLock && val!=0 ) val = 1;
	  if( ch==kTag && data==kSpill ) spill=val;
	  if( ch==kTag && data==kEvent ) evnum=val-1;
	}
      }
      std::cout<<std::setw(10)<<name
	       <<std::setw(10)<<"plane"
	       <<std::setw(5)<<plane<<" : "
	       <<std::setw(10)<<"spill"
	       <<std::setw(8)<<spill<<" : "
	       <<std::setw(15)<<"Event_number"
	       <<std::setw(8)<<evnum
	       <<std::endl;
      event_tag[++index]=evnum;
    }
  }
  {
    const char* name="HUL-RM";
    const int k_device = gUnpacker.get_device_id(name);  
    for(int seg=0;seg<2;seg++){
      int spill=-1;
      int evnum=-1;
      int nhit = gUnpacker.get_entries(k_device, 0, seg, 0, 4);
      if(nhit==1){
	spill = gUnpacker.get(k_device, 0, seg, 0, 4);
      }
      nhit = gUnpacker.get_entries(k_device, 0, seg, 0, 3);
      if(nhit==1){
	evnum = gUnpacker.get(k_device, 0, seg, 0, 3);
      }
      std::cout<<std::setw(10)<<name
	       <<std::setw(10)<<"seg"
	       <<std::setw(5)<<seg<<" : "
	       <<std::setw(10)<<"spill"
	       <<std::setw(8)<<spill<<" : "
	       <<std::setw(15)<<"Event_number"
	       <<std::setw(8)<<evnum
	       <<std::endl;
      event_tag[++index]=evnum;
    }
  }
  {
    const char* name="Scaler";
    const int k_device = gUnpacker.get_device_id(name);  
    for(int plane=0;plane<1;plane++){
      int spill=-1;
      int evnum=-1;
      int nhit = gUnpacker.get_entries(k_device, plane, 0, 0, 2);
      if(nhit==1){
	spill = gUnpacker.get(k_device, plane, 0, 0, 2);
      }
      nhit = gUnpacker.get_entries(k_device, plane, 0, 0, 1);
      if(nhit==1){
	evnum = gUnpacker.get(k_device, plane, 0, 0, 1);
      }
      std::cout<<std::setw(10)<<name
	       <<std::setw(10)<<"plane"
	       <<std::setw(5)<<plane<<" : "
	       <<std::setw(10)<<"spill"
	       <<std::setw(8)<<spill<<" : "
	       <<std::setw(15)<<"Event_number"
	       <<std::setw(8)<<evnum
	       <<std::endl;
      event_tag[++index]=evnum;
    }
  }
  bool SLIP=false;
  for(int i=0;i<6;i++){
    if(event_tag[i]==-1) break;
    if(event_tag[i]!=event_tag[0]){
      SLIP=true;
      break;
    } 
  }
  std::ofstream ofs("/home/oper/.daq/rm.dat");
  if(SLIP){
    std::cout<<std::setw(70)<<"!!!!! Event SLIP !!!!!" <<std::endl;
    std::cout<<std::setw(70)<<"!!!!! Event SLIP !!!!!" <<std::endl;
    ofs<<"!!!!! Event SLIP !!!!!"<<std::endl;
  }else{
    std::cout<<std::setw(70)<<"OK" <<std::endl;
    ofs<<"OK"<<std::endl;
  }
  ofs.close();
  return 0;
}

}
