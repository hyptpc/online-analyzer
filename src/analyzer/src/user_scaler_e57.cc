// -*- C++ -*-

// Author: Tomonori Takahashi

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

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
  }

//____________________________________________________________________________
Int_t
process_begin( const std::vector<std::string>& argv )
{
  ConfMan::getInstance().initialize(argv);

  // flags
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

  //////////////////// Set Channels
  // ScalerAnalylzer::Set( Int_t column,
  //                       Int_t raw,
  //                       ScalerInfo( name, module, channel ) );
  // scaler information is defined from here.
  // please do not use a white space character.
  {
    Int_t c = ScalerAnalyzer::kLeft;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "                  FT",       0,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "                 TM1",       0,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "                 TM2",       0,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "                SYIM",       0,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "                 BHD",       0,  4 ) );
    gScaler.Set( c, r++, ScalerInfo( "                  T0",       0,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "                  AC",       0,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "               T0new",       0,  7 ) );
    gScaler.Set( c, r++, ScalerInfo( "                  E0",       0,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "                 DEF",       0,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "                CDH1",       0, 10 ) );
    gScaler.Set( c, r++, ScalerInfo( "                CDH2",       0, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "                CDH3",       0, 12 ) );
    gScaler.Set( c, r++, ScalerInfo( "                 SDD",       0, 13 ) );
    gScaler.Set( c, r++, ScalerInfo( "           real time",       0, 14 ) );
    gScaler.Set( c, r++, ScalerInfo( "           dead time",       0, 15 ) );
    gScaler.Set( c, r++, ScalerInfo( "       Beam1(BHDxT0)",       0, 16 ) );
    gScaler.Set( c, r++, ScalerInfo( "    Beam2(BHDxT0new)",       0, 17 ) );
    gScaler.Set( c, r++, ScalerInfo( "    Beam3(T0newxDEF)",       0, 18 ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon1",       0, 19 ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon2",       0, 20 ) );
    gScaler.Set( c, r++, ScalerInfo( "               Kaon3",       0, 21 ) );
    gScaler.Set( c, r++, ScalerInfo( "               Pion1",       0, 22 ) );
    gScaler.Set( c, r++, ScalerInfo( "               Pion2",       0, 23 ) );
    gScaler.Set( c, r++, ScalerInfo( "           KaonxCDH1",       0, 24 ) );
    gScaler.Set( c, r++, ScalerInfo( "           KaonxCDH2",       0, 25 ) );
    gScaler.Set( c, r++, ScalerInfo( "           PionxCDH1",       0, 26 ) );
    gScaler.Set( c, r++, ScalerInfo( "           PionxCDH2",       0, 27 ) );
    gScaler.Set( c, r++, ScalerInfo( "             Request",       0, 28 ) );
    gScaler.Set( c, r++, ScalerInfo( "              Accept",       0, 29 ) );
    gScaler.Set( c, r++, ScalerInfo( "                None",       0, 30 ) );
    gScaler.Set( c, r++, ScalerInfo( "          Clock 1kHz",       0, 31 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "SDD gate  1",       1,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  2",       1,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  3",       1,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  4",       1,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  5",       1,  4 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  6",       1,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  7",       1,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate  8",       1,  7 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD gate  9",       1,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 10",       1,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 11",       1, 10 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 12",       1, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 13",       1, 12 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 14",       1, 13 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 15",       1, 14 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 16",       1, 15 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD gate 17",       1, 16 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 18",       1, 17 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 19",       1, 18 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 20",       1, 19 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 21",       1, 20 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 22",       1, 21 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 23",       1, 22 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 24",       1, 23 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 25",       1, 24 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 26",       1, 25 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD gate 27",       1, 26 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 28",       1, 27 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 29",       1, 28 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 30",       1, 29 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 31",       1, 30 ) );
    gScaler.Set( c, r++, ScalerInfo( "    gate 32",       1, 31 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "SDD reset  1",      1, 32 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  2",      1, 33 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  3",      1, 34 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  4",      1, 35 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  5",      1, 36 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  6",      1, 37 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  7",      1, 38 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset  8",      1, 39 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD reset  9",      1, 40 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 10",      1, 41 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 11",      1, 42 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 12",      1, 43 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 13",      1, 44 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 14",      1, 45 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 15",      1, 46 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 16",      1, 47 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD reset 17",      1, 48 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 18",      1, 49 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 19",      1, 50 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 20",      1, 51 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 21",      1, 52 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 22",      1, 53 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 23",      1, 54 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 24",      1, 55 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 25",      1, 56 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 26",      1, 57 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD reset 27",      1, 58 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 28",      1, 59 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 29",      1, 60 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 30",      1, 61 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 31",      1, 62 ) );
    gScaler.Set( c, r++, ScalerInfo( "    reset 32",      1, 63 ) );
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
