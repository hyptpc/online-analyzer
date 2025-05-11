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
    gScaler.Set( c, r++, ScalerInfo( "PbF2  1",       0,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  2",       0,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  3",       0,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  4",       0,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  5",       0,  4 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  6",       0,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  7",       0,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  8",       0,  7 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2  9",       0,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 10",       0,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 11",       0, 10 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 12",       0, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 13",       0, 12 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 14",       0, 13 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 15",       0, 14 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 16",       0, 15 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 17",       0, 16 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 18",       0, 17 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 19",       0, 18 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 20",       0, 19 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "PbF2 21",       0, 20 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 22",       0, 21 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 23",       0, 22 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 24",       0, 23 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 25",       0, 24 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 26",       0, 25 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 27",       0, 26 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 28",       0, 27 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 29",       0, 28 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 30",       0, 29 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 31",       0, 30 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 32",       0, 31 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 33",       0, 32 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 34",       0, 33 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 35",       0, 34 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 36",       0, 35 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 37",       0, 36 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 38",       0, 37 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 39",       0, 38 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 40",       0, 39 ) );
    gScaler.Set( c, r++, ScalerInfo( "PbF2 or",       0, 56 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0; 
    gScaler.Set( c, r++, ScalerInfo( "request",       0, 40 ) );
    gScaler.Set( c, r++, ScalerInfo( "accept",       0, 41 ) );
    gScaler.Set( c, r++, ScalerInfo( "spill start",       0, 42 ) );
    gScaler.Set( c, r++, ScalerInfo( "spill end",       0, 43 ) );
    gScaler.Set( c, r++, ScalerInfo( "beam trig",       0, 44 ) );
    gScaler.Set( c, r++, ScalerInfo( "cosmic trig",       0, 45 ) );
    gScaler.Set( c, r++, ScalerInfo( "clock trig",       0, 46 ) );
    gScaler.Set( c, r++, ScalerInfo( "1kHz",       0, 47 ) );
    gScaler.Set( c, r++, ScalerInfo( "trig_x",       0, 48 ) );
    gScaler.Set( c, r++, ScalerInfo( "tirg_y",       0, 49 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto1-U",       0, 50 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto1-D",       0, 51 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto2-U",       0, 52 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto2-D",       0, 53 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto1-coin",       0, 54 ) );
    gScaler.Set( c, r++, ScalerInfo( "veto2-coin",       0, 55 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma1-U",       0, 58 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma1-D",       0, 59 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma2-U",       0, 60 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma2-D",       0, 61 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma1-coin",       0, 57 ) );
    gScaler.Set( c, r++, ScalerInfo( "sakuma2-coin",       0, 63 ) );
    //    gScaler.Set( c, r++, ScalerInfo( "",       0, 62 ) );
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
