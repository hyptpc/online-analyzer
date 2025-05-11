// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include <THttpServer.h>

#include "user_analyzer.hh"
#include "UnpackerManager.hh"
#include "ConfMan.hh"
#include <Rtypes.h>

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;
//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  ConfMan::getInstance().initialize(argv);
  //  GUnpacker::get_instance().set_decode_mode(false);

  new THttpServer("http:8080");
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
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();
  sleep(1);
  std::cout << "\033[2J" << std::endl;
  std::cout << "#D skeleton::process_event" << std::endl;
  gUnpacker.dump_data_fe(1638);
  gUnpacker.dump_data_fe(1641);
  // int list[2]={513,515};
  // for(int i=0;i<2;i++){
  //   gUnpacker.dump_data_fe(list[i]);
  // }
  // for(int i=0;i<17;i++){
  //   gUnpacker.dump_data_fe(1637+i);
  // }
  // gUnpacker.dump_data_fe(1657);
  // gUnpacker.dump_data_fe(1661);

  // static const int k_device = gUnpacker.get_device_id("SDD2");
  // gUnpacker.dump_data_device(k_device);

  // static const int k_device2 = gUnpacker.get_device_id("V1290");
  // gUnpacker.dump_data_device(k_device2);

  // static const int k_device3 = gUnpacker.get_device_id("Scaler2");
  // gUnpacker.dump_data_device(k_device3);

  // static const int k_device = gUnpacker.get_device_id("Counter2");
  // static const int k_l    = 0;
  // static const int k_t    = 1;
  // static const int NumOfSegments=32;
  // for(int l = 0; l<NumOfSegments; ++l){
  //   int nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_l);
  //   std::cout<<"-----------"<<std::endl;
  //   std::cout<<"leading nhit: "<<nhit<<std::endl;
  //   if( nhit ){
  //     // This wire fired at least one times.
  //     for( int m=0; m<nhit; ++m ){
  // 	int tdc = gUnpacker.get(k_device, 0, l, 0, k_l, m);
  //     }	
  //   }
  //   nhit = gUnpacker.get_entries(k_device, 0, l, 0, k_t);
  //   std::cout<<"trailing nhit: "<<nhit<<std::endl;
  //   if( nhit ){
  //     for( int m=0; m<nhit; ++m ){
  // 	int tdc = gUnpacker.get(k_device, 0, l, 0, k_t , m);
  //     } 	
  //   }
  //   nhit = gUnpacker.get_entries(k_device, 0, l, 0, 2);
  //   std::cout<<"data3 nhit: "<<nhit<<std::endl;
  //   if( nhit ){
  //     for( int m=0; m<nhit; ++m ){
  // 	int tdc = gUnpacker.get(k_device, 0, l, 0, 2 , m);
  //     } 	
  //   }
  // }

  // int module_id=0;
  // for(int channel=0;channel<31;channel++){
  //   for(int data_type=0;data_type<3;data_type++){
  //     Int_t nhit = gUnpacker.get_entries( device_id3, module_id, channel,0, data_type );
  //     if( nhit<=0 ) continue;
  //     if( nhit>2 ){
  // 	for(int ihit=0;ihit<nhit; ihit++){
  // 	  int val = gUnpacker.get( device_id3, module_id, channel, 0,  data_type, ihit );
  // 	  std::cout<<channel<<"  "<<ihit<<" / "<<nhit<<"  "<<data_type<<"  "<<val<<std::endl;
  // 	}
  //     }
  //   }
  // }
  return 0;
}

}
