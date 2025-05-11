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
    //    gScaler.Set( c, r++, ScalerInfo( "FT",     0, 16 ) );
    gScaler.Set( c, r++, ScalerInfo( "FT",       1, 35 ) );
    gScaler.Set( c, r++, ScalerInfo( "SYIM",     0, 17 ) );
    gScaler.Set( c, r++, ScalerInfo( "TM",       0, 18 ) );
    gScaler.Set( c, r++, ScalerInfo( "BHD",      0, 20 ) );
    gScaler.Set( c, r++, ScalerInfo( "T0",       0, 21 ) );
    gScaler.Set( c, r++, ScalerInfo( "AC",  0, 22 ) );
    gScaler.Set( c, r++, ScalerInfo( "T0new",     0, 23 ) );
    gScaler.Set( c, r++, ScalerInfo( "E0",     0, 24 ) );
    gScaler.Set( c, r++, ScalerInfo( "DEF",    0, 25 ) );
    gScaler.Set( c, r++, ScalerInfo( "Start",  0, 26 ) );
    gScaler.Set( c, r++, ScalerInfo( "Stop",   0, 27 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDD",    0, 28 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDDReset",    0, 39 ) );
    gScaler.Set( c, r++, ScalerInfo( "SDDVeto",0, 29 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clock1kHz",   0, 19 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "Beam1",  0, 30 ) );
    gScaler.Set( c, r++, ScalerInfo( "Beam2" , 0, 31 ) );
    gScaler.Set( c, r++, ScalerInfo( "Beam3" , 0, 32 ) );
    gScaler.Set( c, r++, ScalerInfo( "Pion1" , 0, 33 ) );
    gScaler.Set( c, r++, ScalerInfo( "Pion2" , 0, 34 ) );
    gScaler.Set( c, r++, ScalerInfo( "Kaon1",  0, 35 ) );
    gScaler.Set( c, r++, ScalerInfo( "Kaon2",  0, 37 ) );
    gScaler.Set( c, r++, ScalerInfo( "Kaon3",  0, 36 ) );
    // gScaler.Set( c, r++, ScalerInfo( "PixStart" ,   0, 42 ) );
    // gScaler.Set( c, r++, ScalerInfo( "PixStartStop",0, 43 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Beam1" ,      0,  2 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Beam2",       0,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "Beam/f",      0,  4 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Pion",        0,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "Pion/f",      0,  6 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Kaon1",       0,  7 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Kaon2" ,      0,  8 ) );
    // gScaler.Set( c, r++, ScalerInfo( "Kaon3" ,      0,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "Kaon2/f" ,    0, 10) );
    gScaler.Set( c, r++, ScalerInfo( "Kaon3/f" ,    0, 11 ) );
    gScaler.Set( c, r++, ScalerInfo( "StartxStop",  0, 38 ) );
    gScaler.Set( c, r++, ScalerInfo( "KxStart" ,    0, 42 ) );
    gScaler.Set( c, r++, ScalerInfo( "KxStartStop" ,0, 43 ) );
    // gScaler.Set( c, r++, ScalerInfo( "KxStart" ,    0, 12 ) );
    // gScaler.Set( c, r++, ScalerInfo( "KxStartxStop",0, 13 ) );
    // gScaler.Set( c, r++, ScalerInfo( "StartxStop" , 0, 14 ) );
    //    gScaler.Set( c, r++, ScalerInfo( "AC" ,         0, 15 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    gScaler.Set( c, r++, ScalerInfo( "Request" ,    0, 44 ) );
    gScaler.Set( c, r++, ScalerInfo( "Accept" ,     0, 45 ) );
    //    gScaler.Set( c, r++, ScalerInfo( "real time" ,  0, 46 ) );
    //    gScaler.Set( c, r++, ScalerInfo( "live time" ,  0, 47 ) );
    gScaler.Set( c, r++, ScalerInfo( "SpillStart" , 0,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "SpillEnd" ,   0,  1 ) ); 
    gScaler.Set( c, r++, ScalerInfo( "Gate1",    1,  0 ) );
    gScaler.Set( c, r++, ScalerInfo( "Gate2",    1,  1 ) );
    gScaler.Set( c, r++, ScalerInfo( "Gate3",   1,  2 ) );
    gScaler.Set( c, r++, ScalerInfo( "Gate4",  1,  3 ) );
    gScaler.Set( c, r++, ScalerInfo( "Reset1",    1,  4 ) );
    gScaler.Set( c, r++, ScalerInfo( "Reset2",    1,  5 ) );
    gScaler.Set( c, r++, ScalerInfo( "Reset3", 1,  6 ) );
    gScaler.Set( c, r++, ScalerInfo( "Reset4",    1,  7 ) );    
    gScaler.Set( c, r++, ScalerInfo( "Clock_480ns",    1,  8 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clock_7p2us" ,   1,  9 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clock_480ns",    1,  32 ) );
    gScaler.Set( c, r++, ScalerInfo( "Clock_7p2us" ,   1,  33 ) );
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
