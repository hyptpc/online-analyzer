// -*- C++ -*-

// Author: Tomonori Takahashi

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <TSystem.h>

#include "DAQNode.hh"
#include "UnpackerManager.hh"
#include "Unpacker.hh"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "PrintHelper.hh"
#include "ScalerAnalyzer.hh"
#include "user_analyzer.hh"

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    UnpackerManager& gUnpacker = GUnpacker::get_instance();
    ScalerAnalyzer&  gScaler   = ScalerAnalyzer::GetInstance();
    ///// Number of spill for Scaler Sheet
    Scaler nspill_scaler_sheet  = 1;
    const std::string filename=Form("%s/k18br_analyzer/e73/param/CMAP/scaler_e73_202404.param",std::getenv("HOME"));
  }

//____________________________________________________________________________
Int_t
process_begin( const std::vector<std::string>& argv )
{
  ConfMan::getInstance().initialize(argv);

  // flags
  gScaler.SetNRows(17);
  gScaler.SetFlag( ScalerAnalyzer::kSeparateComma );
  gScaler.SetFlag( ScalerAnalyzer::kSemiOnline );
  for( Int_t i=0, n=argv.size(); i<n; ++i ){
    TString v = argv[i];
    if( v.Contains("--print") ){
      gScaler.SetFlag( ScalerAnalyzer::kScalerSheet );
    }
    if( v.Contains("--spill=") ){
      nspill_scaler_sheet = v.ReplaceAll("--spill=","").Atoi();
    }
    if( v.Contains(":") ){
      gScaler.SetFlag( ScalerAnalyzer::kSemiOnline, false );
    }
    if( v.Contains("--spill-by-spill") ){
      gScaler.SetFlag( ScalerAnalyzer::kSpillBySpill );
    }
  }

  std::cout<<filename<<std::endl;
  std::ifstream ifs(filename);
  if(!ifs){
    std::cout<<"#E "<<filename<<" cannot be opend"<<std::endl;
    exit(0);
  }
  std::cout<<filename<<" opened"<<std::endl;
  std::string name;
  int plane,channel;
  int column,row;
  while(ifs>>name>>plane>>channel>>column>>row){
#if 1
    std::cout<<std::setw(20)<<name
             <<std::setw(5)<<plane
             <<std::setw(5)<<channel
             <<std::setw(5)<<column
             <<std::setw(5)<<row
             <<std::endl;
#endif
    if(column<0||row<0) continue;
    gScaler.Set( column, row, ScalerInfo(name,plane,channel) );
  }
  std::cout<<""<<__func__<<" finished"<<std::endl;
  ifs.close();
  return 0;
  //////////////////// Set Channels
  // ScalerAnalylzer::Set( Int_t column,
  //                       Int_t raw,
  //                       ScalerInfo( name, module, channel ) );
  // scaler information is defined from here.
  // please do not use a white space character.
  Int_t i = 0;
  {
    Int_t c = ScalerAnalyzer::kLeft;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "                  FT",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                 TM1",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                SYIM",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                 BHD",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                  T0",       0, i++ ) );//4
    gScaler.Set( c, r++, ScalerInfo( "                  AC",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               T0new",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                 DEF",       0, i++ ) );
    // gScaler.Set( c, r++, ScalerInfo( "               Veto0",       0, i++ ) );//8
    // gScaler.Set( c, r++, ScalerInfo( "               Veto1",       0, i++ ) );
    // gScaler.Set( c, r++, ScalerInfo( "                PbF2",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                Veto",       0, i++ ) );//8
    gScaler.Set( c, r++, ScalerInfo( "                PbF2",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                 BTC",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                CDH1",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "                CDH2",       0, i++ ) );//12
    gScaler.Set( c, r++, ScalerInfo( "                CDH3",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           real time",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           dead time",       0, i++ ) );
  }
  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "       Beam1(BHDxT0)",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "    Beam2(BHDxT0new)",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "    Beam3(T0newxDEF)",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon1",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon2",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon3",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               Pion1",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "               Pion2",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           KaonxCDH1",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           KaonxCDH2",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           KaonxCDH3",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "     KaonxCDH1xGamma",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "            PionxCDH",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "             Request",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "              Accept",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "         Clock 10kHz",       0, i++ ) );
  }
  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "         Spill start",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "           Spill end",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "              Beam/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "              Pion/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "             Kaon1/f",       0, i++ ) ); //4
    gScaler.Set( c, r++, ScalerInfo( "             Kaon2/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "         KaonxCDH1/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "         KaonxCDH2/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "         KaonxCDH3/f",       0, i++ ) );//8
    gScaler.Set( c, r++, ScalerInfo( "     KaonxCDH1xGamma",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "          KaonxGamma",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "            PionxCDH",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "         PionxPbF2/f",       0, i++ ) );//12
    gScaler.Set( c, r++, ScalerInfo( "     ElectronxPbF2/f",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "          CDH cosmic",       0, i++ ) );
    gScaler.Set( c, r++, ScalerInfo( "        no beam(10s)",       0, i++ ) );
  }

  gScaler.PrintFlags();

  return 0;
}

//____________________________________________________________________________
Int_t
process_end( void )
{
  if( gScaler.GetFlag( ScalerAnalyzer::kScalerSheet ) )
    return 0;

  gScaler.Print();

  return 0;

}

//____________________________________________________________________________
Int_t
process_event( void )
{
  static Int_t  event_count = 0;
  static Bool_t en_disp     = false;

  //  std::cout<<gUnpacker.get_root()->get_run_number()<<std::endl;
  if( gScaler.GetFlag( ScalerAnalyzer::kScalerSheet ) && event_count==0 )
    std::cout << "waiting spill end " << std::flush;

  ++event_count;

  gScaler.Decode();

  if( gScaler.GetFlag( ScalerAnalyzer::kSemiOnline ) ){
    if( event_count%300 == 0 ) en_disp = true;
  } else {
    if( event_count%10 == 0 ) en_disp = true;
  }

  if( gScaler.IsSpillEnd() )
    en_disp = true;

  if( gScaler.GetFlag( ScalerAnalyzer::kScalerSheet ) &&
      !gScaler.IsSpillEnd() ){
    if( event_count%100==0 )
      std::cout << "." << std::flush;
    return 0;
  }

  if( en_disp ){
    gScaler.Print();
    en_disp = false;
  }

  if( gScaler.IsSpillEnd() &&
      gScaler.GetFlag( ScalerAnalyzer::kScalerSheet ) ){
    std::cout << "found spill end "
    	      << gScaler.Get("Spill") << "/" << nspill_scaler_sheet
	      << std::endl;

    if( gScaler.Get("Spill") == nspill_scaler_sheet ){
      gScaler.PrintScalerSheet();
      return -1;
    }

    if( gScaler.Get("Spill") > nspill_scaler_sheet ){
      std::cout << "something is wrong!" << std::endl;
      return -1;
    }

    std::cout << "waiting spill end " << std::flush;
  }

  return 0;
}

}
