// -*- C++ -*-

#include "PsMaker.hh"
#include "GuiPs.hh"
#include "HistMaker.hh"
#include "DetectorID.hh"
#include "HistHelper.hh"
#include "UnpackerManager.hh"
#include "Unpacker.hh"

#include <iostream>
#include <algorithm>

#include <TROOT.h>
#include <TList.h>
#include <TIterator.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPad.h>
#include <TPostScript.h>
#include <TTimeStamp.h>
#include <TH1.h>
#include <TH2.h>
#include <TText.h>
#include <TString.h>
#include <TStyle.h>

#define CONV_STRING(x) getStr_FromEnum(#x)

namespace
{
  const std::string MyName = "PsMaker::";
  using namespace hddaq::gui;
}

// Constructor --------------------------------------------------------------
PsMaker::PsMaker():
  cps_(NULL), ps_(NULL)
{
  // Initialize option list
  name_option_.resize(sizeOptionList);

  name_option_[kExpDataSheet]  = "ExpDataSheet";
  name_option_[kFixXaxis]      = "FixXaxis";
  name_option_[kLogyADC]       = "LogyADC";
  name_option_[kLogyTDC]       = "LogyTDC";
  //  name_option_[kLogyHitPat]    = "LogyHitPat";
  //  name_option_[kLogyMulti]     = "LogyMulti";
}

// Destructor --------------------------------------------------------------
PsMaker::~PsMaker()
{

}

// -------------------------------------------------------------------------
// getListOfOption
// -------------------------------------------------------------------------
void PsMaker::getListOfOption(std::vector<std::string>& vec)
{
  PsMaker& g = PsMaker::getInstance();
  std::copy(g.name_option_.begin(),
	    g.name_option_.end(),
	    back_inserter(vec)
	    );
}

// -------------------------------------------------------------------------
// makePs
// -------------------------------------------------------------------------
void PsMaker::makePs()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1110);
  gStyle->SetTitleW(.4);
  gStyle->SetTitleH(.1);
  gStyle->SetStatW(.32);
  gStyle->SetStatH(.25);

  // make ps file instance
  const std::string& filename = GuiPs::getFilename();
  ps_  = new TPostScript(filename.c_str(), 112);
  cps_ = new TCanvas("cps","");

  std::cout << std::endl << "PSFile = " << filename
	    << std::endl << std::endl;

  // Make title page with run number
  drawRunNumber();

  // start draw histograms of each detector
  std::vector<std::string> name_detectors;
  HistMaker::getListOfPsFiles(name_detectors);
  const int NofDetector = name_detectors.size();
  for(int i = 0; i<NofDetector; ++i){
    if(GuiPs::isDevOn(i) || GuiPs::isOptOn(kExpDataSheet)){
      create(name_detectors[i]);
    }
  }

  // close
  if(ps_){
    ps_->Close();
    ps_ = NULL;
  }

  if(cps_){
    delete cps_;
    cps_ = NULL;
  }

  gROOT->SetStyle("Classic");
}

// -------------------------------------------------------------------------
// create
// -------------------------------------------------------------------------
void PsMaker::create(std::string& name)
{
  std::vector<int> id_list;
  std::vector<int> par_list(sizeParameterList);
  bool flag_xaxis = false; // case by case in each detector
  bool flag_log   = false; // case by case in each detector

  // BH1 ----------------------------------------------------------------
  if(name == CONV_STRING(kBH1)){
    int base_id = 0;
    // ADC
    par_list[kXdiv] = 3; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    // ADC U
    base_id = HistMaker::getUniqueID(kBH1, 0, kADC);
    for(int i = 0; i<NumOfSegBH1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // ADC D
    base_id = HistMaker::getUniqueID(kBH1, 0, kADC, NumOfSegBH1+1);
    for(int i = 0; i<NumOfSegBH1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 3; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC U
    base_id = HistMaker::getUniqueID(kBH1, 0, kTDC);
    for(int i = 0; i<NumOfSegBH1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC D
    base_id = HistMaker::getUniqueID(kBH1, 0, kTDC, NumOfSegBH1+1);
    for(int i = 0; i<NumOfSegBH1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kBH1, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kBH1, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // BFT ----------------------------------------------------------------
  if(name == CONV_STRING(kBFT)){
    // TDC
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    // par_list[kXrange_min] = 0; par_list[kXrange_max] = 1000; // for TDC
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC,      1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC,      2));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC,     11));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC,     12));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC2D,    1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC2D,    2));
    GHist::get( HistMaker::getUniqueID(kBFT, 0, kTDC2D, 1) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    GHist::get( HistMaker::getUniqueID(kBFT, 0, kTDC2D, 2) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    drawOneCanvas(id_list, par_list, flag_xaxis, false, "colz");

    // TOT
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    // par_list[kXrange_min] = -50; par_list[kXrange_max] = 150; // for TOT
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC,      1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC,      2));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC,     11));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC,     12));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC2D,    1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC2D,    2));
    drawOneCanvas(id_list, par_list, flag_xaxis, false, "colz");

    // HitPat & multi
    par_list[kXdiv] = 3; par_list[kYdiv] = 2;
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kHitPat,   1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kHitPat,   2));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kMulti,    1));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kHitPat,  11));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kHitPat,  12));
    id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kMulti,   11));
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    // id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC2D,    1));
    // id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kTDC2D,    2));
    // id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC2D,    1));
    // id_list.push_back(HistMaker::getUniqueID(kBFT, 0, kADC2D,    2));
    // drawOneCanvas(id_list, par_list, flag_xaxis, false, "box");
  }

  // BC3 ----------------------------------------------------------------
  if(name == CONV_STRING(kBC3)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // BC3 TDC
    par_list[kXrange_min] = 256; par_list[kXrange_max] = 1000;
    base_id = HistMaker::getUniqueID(kBC3, 0, kTDC);
    for(int i = 0; i<NumOfLayersBC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // BC3 HitPat
    base_id = HistMaker::getUniqueID(kBC3, 0, kHitPat);
    for(int i = 0; i<NumOfLayersBC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // BC3 Multi
    base_id = HistMaker::getUniqueID(kBC3, 0, kMulti);
    for(int i = 0; i<NumOfLayersBC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

  // BC4 ----------------------------------------------------------------
  if(name == CONV_STRING(kBC4)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // BC4 TDC
    par_list[kXrange_min] = 256; par_list[kXrange_max] = 1000;
    base_id = HistMaker::getUniqueID(kBC4, 0, kTDC);
    for(int i = 0; i<NumOfLayersBC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // BC4 HitPat
    base_id = HistMaker::getUniqueID(kBC4, 0, kHitPat);
    for(int i = 0; i<NumOfLayersBC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // BC4 Multi
    base_id = HistMaker::getUniqueID(kBC4, 0, kMulti);
    for(int i = 0; i<NumOfLayersBC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

//  // BMW ----------------------------------------------------------------
//  if(name == CONV_STRING(kBMW)){
//    // For all
//    par_list[kXdiv] = 2; par_list[kYdiv] = 1;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kLogyTDC) | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC && TDC
//    id_list.push_back(HistMaker::getUniqueID(kBMW, 0, kADC));
//    id_list.push_back(HistMaker::getUniqueID(kBMW, 0, kTDC));
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }

  // BH2 ----------------------------------------------------------------
  if(name == CONV_STRING(kBH2)){
    int base_id = 0;
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;

    // ADC
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    // ADC Up
    base_id = HistMaker::getUniqueID(kBH2, 0, kADC);
    for(int i = 0; i<NumOfSegBH2; ++i){id_list.push_back(base_id + i);}
    // ADC Down
    base_id = HistMaker::getUniqueID(kBH2, 0, kADC, NumOfSegBH2 + 1);
    for(int i = 0; i<NumOfSegBH2; ++i){id_list.push_back(base_id + i);}

    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    // TDC
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC Up
    base_id = HistMaker::getUniqueID(kBH2, 0, kTDC);
    for(int i = 0; i<NumOfSegBH2; ++i){id_list.push_back(base_id + i);}
    // TDC Down
    base_id = HistMaker::getUniqueID(kBH2, 0, kTDC, NumOfSegBH2 + 1);
    for(int i = 0; i<NumOfSegBH2; ++i){id_list.push_back(base_id + i);}

    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    id_list.push_back(HistMaker::getUniqueID(kBH2, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kBH2, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // CFT ----------------------------------------------------------------
  if(name == CONV_STRING(kCFT)){
    int base_id = 0;
    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kTDC, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kTDC2D ,1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // TOT
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kADC, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TOT_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kADC2D ,1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CTOT
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kADC, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // CTOT_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kADC2D ,11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // HighGain
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kHighGain, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // Pede 
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kPede, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // LowGain
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kCFT, 0, kLowGain, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HighGain_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHighGain, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // Pede_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kPede, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // LowGain_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kLowGain, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CHighGain
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHighGain, 21);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CLowGain
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kLowGain, 21);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CHighGain_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHighGain, 31);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CLowGain_2D
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kLowGain, 31);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CHighGain_x_TOT
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHighGain, 41);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CLowGain_x_TOT
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kLowGain, 41);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // HitPattern 
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHitPat, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CHitPattern 
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kHitPat, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // Multiplicity 
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kMulti, 1);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // CMultiplicity 
    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kCFT, 0, kMulti, 11);
    for(int i = 0; i<NumOfLayersCFT; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // HitPattern Multiplicity 
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    id_list.push_back(HistMaker::getUniqueID(kCFT, 0, kMulti,  9));
    id_list.push_back(HistMaker::getUniqueID(kCFT, 0, kMulti, 10));
    id_list.push_back(HistMaker::getUniqueID(kCFT, 0, kMulti, 19));
    id_list.push_back(HistMaker::getUniqueID(kCFT, 0, kMulti, 20));
    drawOneCanvas(id_list, par_list, false, false);

  }

  // BGO ----------------------------------------------------------------
  if(name == CONV_STRING(kBGO)){
    int base_id = 0;
    // FADC  1-12       
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    base_id = HistMaker::getUniqueID(kBGO, 0, kFADC);
    for(int i = 0; i<12; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
    // FADC  13-24         
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    for(int i = 0; i<12; ++i){id_list.push_back(base_id + i + 12);}
    drawOneCanvas(id_list, par_list, false, false);

  }
  // PiID ----------------------------------------------------------------
  if(name == CONV_STRING(kPiID)){
    int base_id = 0;
    // ADC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    // ADC 1-16
    base_id = HistMaker::getUniqueID(kPiID, 0, kADC);
    for(int i = 0; i<16; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // ADC 17-32 
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    base_id = HistMaker::getUniqueID(kPiID, 0, kADC, 16);
    for(int i = 0; i<16; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC 1-16 
    base_id = HistMaker::getUniqueID(kPiID, 0, kTDC);
    for(int i = 0; i<16; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC 17-32 
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    base_id = HistMaker::getUniqueID(kPiID, 0, kTDC, 16);
    for(int i = 0; i<16; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kPiID, 0, kHitPat);
    id_list.push_back(base_id);
    id_list.push_back(base_id + 1);
    base_id = HistMaker::getUniqueID(kPiID, 0, kMulti);
    id_list.push_back(base_id);
    id_list.push_back(base_id + 1);
    drawOneCanvas(id_list, par_list, false, false);
  }
  // SFT ----------------------------------------------------------------
  if(name == CONV_STRING(kSFT)){
    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    // par_list[kXrange_min] = 0; par_list[kXrange_max] = 1000; // for TDC
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,      3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,      4));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,      1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,      2));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,     13));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,     14));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,     11));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC,     12));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    4));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    2));
    GHist::get( HistMaker::getUniqueID(kSFT, 0, kTDC2D, 3) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    GHist::get( HistMaker::getUniqueID(kSFT, 0, kTDC2D, 4) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    GHist::get( HistMaker::getUniqueID(kSFT, 0, kTDC2D, 1) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    GHist::get( HistMaker::getUniqueID(kSFT, 0, kTDC2D, 2) )
      ->GetYaxis()->SetRangeUser( 400., 1000. );
    drawOneCanvas(id_list, par_list, flag_xaxis, false, "colz");

    // TOT
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    // par_list[kXrange_min] = -50; par_list[kXrange_max] = 150; // for TOT
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,      3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,      4));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,      1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,      2));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,     13));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,     14));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,     11));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC,     12));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    4));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    2));
    drawOneCanvas(id_list, par_list, flag_xaxis, false, "colz");

    // HitPat & multi
    par_list[kXdiv] = 5; par_list[kYdiv] = 3;
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,   3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,   4));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,   1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,   2));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,    2));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,    3));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,    1));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,  13));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,  14));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,  11));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kHitPat,  12));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,   12));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,   13));
    id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kMulti,   11));
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    // id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    1));
    // id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kTDC2D,    2));
    // id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    1));
    // id_list.push_back(HistMaker::getUniqueID(kSFT, 0, kADC2D,    2));
    // drawOneCanvas(id_list, par_list, flag_xaxis, false, "box");
  }

//  // BAC ----------------------------------------------------------------
//  if(name == CONV_STRING(kBAC)){
//    int base_id = 0;
//
//    par_list[kXdiv] = 3; par_list[kYdiv] = 2;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC
//    base_id = HistMaker::getUniqueID(kBAC, 0, kADC);
//    for(int i = 0; i<NumOfSegBAC; ++i){id_list.push_back(base_id + i);}
//    // HitPat
//    base_id = HistMaker::getUniqueID(kBAC, 0, kHitPat);
//    id_list.push_back(base_id);
//    // TDC
//    base_id = HistMaker::getUniqueID(kBAC, 0, kTDC);
//    for(int i = 0; i<NumOfSegBAC; ++i){id_list.push_back(base_id + i);}
//    // Multi
//    base_id = HistMaker::getUniqueID(kBAC, 0, kMulti);
//    id_list.push_back(base_id);
//
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // BH2_E07 ----------------------------------------------------------------
//  if(name == CONV_STRING(kBH2_E07)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    base_id = HistMaker::getUniqueID(kBH2_E07, 0, kADC);
//    for(int i = 0; i<NumOfSegBH2_E07*2; ++i){id_list.push_back(base_id + i);}
//    base_id = HistMaker::getUniqueID(kBH2_E07, 0, kTDC);
//    for(int i = 0; i<NumOfSegBH2_E07*2; ++i){id_list.push_back(base_id + i);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // HitPat/Multi
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kBH2_E07, 0, kHitPat);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kBH2_E07, 0, kMulti);
//    id_list.push_back(base_id);
//    drawOneCanvas(id_list, par_list, false, false);
//  }
//
//  // BAC_E07 ----------------------------------------------------------------
//  if(name == CONV_STRING(kBAC_E07)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 3; par_list[kYdiv] = 4;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    base_id = HistMaker::getUniqueID(kBAC_E07, 0, kADC);
//    id_list.push_back(base_id + 0);
//    id_list.push_back(base_id + 2);
//    id_list.push_back(base_id + 3);
//    id_list.push_back(base_id + 1);
//    id_list.push_back(base_id + 4);
//    id_list.push_back(base_id + 5);
//    base_id = HistMaker::getUniqueID(kBAC_E07, 0, kTDC);
//    id_list.push_back(base_id + 0);
//    id_list.push_back(base_id + 2);
//    id_list.push_back(base_id + 3);
//    id_list.push_back(base_id + 1);
//    id_list.push_back(base_id + 4);
//    id_list.push_back(base_id + 5);
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // HitPat/Multi
//    // par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    // base_id = HistMaker::getUniqueID(kBAC_E07, 0, kHitPat);
//    // id_list.push_back(base_id);
//    // base_id = HistMaker::getUniqueID(kBAC_E07, 0, kMulti);
//    // id_list.push_back(base_id);
//    // drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // FBH ----------------------------------------------------------------
//  if(name == CONV_STRING(kFBH)){
//    int base_id = 0;
//    // TDC
//    par_list[kXdiv] = 8; par_list[kYdiv] = 4;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 1000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    base_id    = HistMaker::getUniqueID(kFBH, 0, kTDC);
//    for(int i=0; i<NumOfSegFBH; ++i){id_list.push_back(base_id +i);}
//    for(int i=0; i<NumOfSegFBH; ++i){id_list.push_back(base_id +i +NumOfSegFBH);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, false);
//
//    // TOT
//    par_list[kXdiv] = 8; par_list[kYdiv] = 4;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kADC);
//    for(int i=0; i<NumOfSegFBH; ++i){id_list.push_back(base_id +i);}
//    for(int i=0; i<NumOfSegFBH; ++i){id_list.push_back(base_id +i +NumOfSegFBH);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // TDC/TOT 2D, HitPat and Multi
//    par_list[kXdiv] = 3; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kFBH, 0, kTDC, NumOfSegFBH*2+1);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kADC, NumOfSegFBH*2+1);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kHitPat);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kTDC2D);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kADC2D);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFBH, 0, kMulti);
//    id_list.push_back(base_id);
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//  }
//
//  // SSD0 ----------------------------------------------------------------
//  if(name == CONV_STRING(kSSD0)){
//    int base_id = 0;
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = false;
//
//    // ADC
//    base_id = HistMaker::getUniqueID(kSSD0, 0, kADC2D);
//    for(int i=0; i<NumOfLayersSSD0; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    // TDC
//    base_id = HistMaker::getUniqueID(kSSD0, 0, kTDC2D);
//    for(int i=0; i<NumOfLayersSSD0; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // HitPat/Multi
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD0, 0, kHitPat);
//    for(int i=0; i<NumOfLayersSSD0; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD0, 0, kMulti);
//    for(int i=0; i<NumOfLayersSSD0; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // SSD1 ----------------------------------------------------------------
//  if(name == CONV_STRING(kSSD1)){
//    int base_id = 0;
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = false;
//
//    // ADC
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kADC2D);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//    // TDC
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kTDC2D);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // HitPat/Multi
//    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kHitPat);
//    for(int i=0; i<NumOfLayersSSD1; ++i){
//      id_list.push_back(base_id+2*i);
//    }
//    for(int i=0; i<NumOfLayersSSD1; ++i){
//      id_list.push_back(base_id+2*i+1);
//    }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kMulti);
//    for(int i=0; i<NumOfLayersSSD1; ++i){
//      id_list.push_back(base_id+2*i);
//    }
//    for(int i=0; i<NumOfLayersSSD1; ++i){
//      id_list.push_back(base_id+2*i+1);
//    }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // DeltaE
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kDeltaE);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kDeltaE2D);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // CTime
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kCTime);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kCTime2D);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // Chisqr
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD1, 0, kChisqr);
//    for(int i=0; i<NumOfLayersSSD1; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // SSD2 ----------------------------------------------------------------
//  if(name == CONV_STRING(kSSD2)){
//    int base_id = 0;
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = false;
//
//    // ADC
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kADC2D);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//    // TDC
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kTDC2D);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // HitPat/Multi
//    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kHitPat);
//    for(int i=0; i<NumOfLayersSSD2; ++i){
//      id_list.push_back(base_id+2*i);
//    }
//    for(int i=0; i<NumOfLayersSSD2; ++i){
//      id_list.push_back(base_id+2*i+1);
//    }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kMulti);
//    for(int i=0; i<NumOfLayersSSD2; ++i){
//      id_list.push_back(base_id+2*i);
//    }
//    for(int i=0; i<NumOfLayersSSD2; ++i){
//      id_list.push_back(base_id+2*i+1);
//    }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // DeltaE
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kDeltaE);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kDeltaE2D);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // CTime
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kCTime);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kCTime2D);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
//
//    // Chisqr
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kSSD2, 0, kChisqr);
//    for(int i=0; i<NumOfLayersSSD2; ++i){ id_list.push_back(base_id++); }
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // PVAC ----------------------------------------------------------------
//  if(name == CONV_STRING(kPVAC)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    base_id = HistMaker::getUniqueID(kPVAC, 0, kADC);
//    for(int i = 0; i<NumOfSegPVAC; ++i){id_list.push_back(base_id + i);}
//    base_id = HistMaker::getUniqueID(kPVAC, 0, kTDC);
//    for(int i = 0; i<NumOfSegPVAC; ++i){id_list.push_back(base_id + i);}
//
//    // HitPat/Multi
//    base_id = HistMaker::getUniqueID(kPVAC, 0, kHitPat);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kPVAC, 0, kMulti);
//    id_list.push_back(base_id);
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // FAC ----------------------------------------------------------------
//  if(name == CONV_STRING(kFAC)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    base_id = HistMaker::getUniqueID(kFAC, 0, kADC);
//    for(int i = 0; i<NumOfSegFAC; ++i){id_list.push_back(base_id + i);}
//    base_id = HistMaker::getUniqueID(kFAC, 0, kTDC);
//    for(int i = 0; i<NumOfSegFAC; ++i){id_list.push_back(base_id + i);}
//
//    // HitPat/Multi
//    base_id = HistMaker::getUniqueID(kFAC, 0, kHitPat);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kFAC, 0, kMulti);
//    id_list.push_back(base_id);
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }
//
//  // SAC1 ----------------------------------------------------------------
//  if(name == CONV_STRING(kSAC1)){
//    // ADC
//    par_list[kXdiv] = 2; par_list[kYdiv] = 1;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    id_list.push_back( HistMaker::getUniqueID(kSAC1, 0, kADC) );
//    id_list.push_back( HistMaker::getUniqueID(kSAC1, 0, kTDC) );
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // HitPat/Multi
//    // par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    // id_list.push_back( HistMaker::getUniqueID(kSAC1, 0, kHitPat) );
//    // id_list.push_back( HistMaker::getUniqueID(kSAC1, 0, kMulti) );
//    // drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }

  // SCH ----------------------------------------------------------------
  if(name == CONV_STRING(kSCH)){
    int base_id = 0;
    // TDC
    par_list[kXdiv] = 8; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 1000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    base_id    = HistMaker::getUniqueID(kSCH, 0, kTDC);
    for(int i=0; i<NumOfSegSCH/2; ++i){id_list.push_back(base_id +i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
    for(int i=0; i<NumOfSegSCH/2; ++i){id_list.push_back(base_id +i +NumOfSegSCH/2);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TOT
    par_list[kXdiv] = 8; par_list[kYdiv] = 4;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    base_id = HistMaker::getUniqueID(kSCH, 0, kADC);
    for(int i=0; i<NumOfSegSCH/2; ++i){id_list.push_back(base_id +i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
    for(int i=0; i<NumOfSegSCH/2; ++i){id_list.push_back(base_id +i +NumOfSegSCH/2);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC/TOT 2D, HitPat and Multi
    par_list[kXdiv] = 3; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kSCH, 0, kTDC, NumOfSegSCH+1);
    id_list.push_back(base_id);
    base_id = HistMaker::getUniqueID(kSCH, 0, kADC, NumOfSegSCH+1);
    id_list.push_back(base_id);
    base_id = HistMaker::getUniqueID(kSCH, 0, kHitPat);
    id_list.push_back(base_id);
    base_id = HistMaker::getUniqueID(kSCH, 0, kTDC2D);
    id_list.push_back(base_id);
    base_id = HistMaker::getUniqueID(kSCH, 0, kADC2D);
    id_list.push_back(base_id);
    base_id = HistMaker::getUniqueID(kSCH, 0, kMulti);
    id_list.push_back(base_id);
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "colz");
  }

//  // KIC ----------------------------------------------------------------
//  if(name == CONV_STRING(kKIC)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 4; par_list[kYdiv] = 2;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC/TDC
//    base_id = HistMaker::getUniqueID(kKIC, 0, kADC);
//    for(int i = 0; i<NumOfSegKIC; ++i){id_list.push_back(base_id + i);}
//    base_id = HistMaker::getUniqueID(kKIC, 0, kTDC);
//    for(int i = 0; i<NumOfSegKIC; ++i){id_list.push_back(base_id + i);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // HitPat/Multi
//    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
//    base_id = HistMaker::getUniqueID(kKIC, 0, kHitPat);
//    id_list.push_back(base_id);
//    base_id = HistMaker::getUniqueID(kKIC, 0, kMulti);
//    id_list.push_back(base_id);
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//  }

  // SDC1 ----------------------------------------------------------------
  if(name == CONV_STRING(kSDC1)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // SDC1 TDC
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 1000;
    base_id = HistMaker::getUniqueID(kSDC1, 0, kTDC);
    for(int i = 0; i<NumOfLayersSDC1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // SDC1 HitPat
    base_id = HistMaker::getUniqueID(kSDC1, 0, kHitPat);
    for(int i = 0; i<NumOfLayersSDC1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // SDC1 Multi
    base_id = HistMaker::getUniqueID(kSDC1, 0, kMulti);
    for(int i = 0; i<NumOfLayersSDC1; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

  // SAC_E40 ------------------------------------------------------------
  if(name == CONV_STRING(kSAC)){
    int base_id = 0;

    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    // ADC
    base_id = HistMaker::getUniqueID(kSAC, 0, kADC);
    for(int i = 0; i<NumOfRoomsSAC; ++i){id_list.push_back(base_id + i);}
    // TDC
    base_id = HistMaker::getUniqueID(kSAC, 0, kTDC);
    for(int i = 0; i<NumOfRoomsSAC; ++i){id_list.push_back(base_id + i);}
    // HitPat
    base_id = HistMaker::getUniqueID(kSAC, 0, kHitPat);
    id_list.push_back(base_id);
    // Multi
    base_id = HistMaker::getUniqueID(kSAC, 0, kMulti);
    id_list.push_back(base_id);

    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
  }

  // SDC2 ----------------------------------------------------------------
  if(name == CONV_STRING(kSDC2)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // SDC2 TDC
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 1500;
    base_id = HistMaker::getUniqueID(kSDC2, 0, kTDC);
    for(int i = 0; i<NumOfLayersSDC2; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // SDC2 HitPat
    base_id = HistMaker::getUniqueID(kSDC2, 0, kHitPat);
    for(int i = 0; i<NumOfLayersSDC2; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // SDC2 Multi
    base_id = HistMaker::getUniqueID(kSDC2, 0, kMulti);
    for(int i = 0; i<NumOfLayersSDC2; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

//  // HDC ----------------------------------------------------------------
//  if(name == CONV_STRING(kHDC)){
//    // For all
//    int base_id = 0;
//    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//
//    // HDC TDC
//    par_list[kXrange_min] = 256; par_list[kXrange_max] = 1000;
//    base_id = HistMaker::getUniqueID(kHDC, 0, kTDC);
//    for(int i = 0; i<NumOfLayersHDC; ++i){id_list.push_back(base_id + i);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, false);
//
//    // HDC HitPat
//    base_id = HistMaker::getUniqueID(kHDC, 0, kHitPat);
//    for(int i = 0; i<NumOfLayersHDC; ++i){id_list.push_back(base_id + i);}
//    drawOneCanvas(id_list, par_list, false, false);
//
//    // HDC Multi
//    base_id = HistMaker::getUniqueID(kHDC, 0, kMulti);
//    for(int i = 0; i<NumOfLayersHDC; ++i){id_list.push_back(base_id + i);}
//    drawOneCanvas(id_list, par_list, false, false);
//  }
//
//  // SP0 -----------------------------------------------------------------
//  if(name == CONV_STRING(kSP0)){
//    int base_id = 0;
//    // ADC
//    par_list[kXdiv] = 5; par_list[kYdiv] = 4;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    for(int draw = 0; draw<4; ++draw){
//      for(int l = 0; l<2; ++l){
//	// ADC U
//	base_id = HistMaker::getUniqueID(kSP0, kSP0_L1+l+2*draw, kADC);
//	for(int i = 0; i<NumOfSegSP0; ++i){id_list.push_back(base_id + i);}
//
//	// ADC D
//	base_id = HistMaker::getUniqueID(kSP0, kSP0_L1+l+2*draw, kADC, NumOfSegSP0+1);
//	for(int i = 0; i<NumOfSegSP0; ++i){id_list.push_back(base_id + i);}
//      }
//      drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    }
//
//    // TDC
//    par_list[kXdiv] = 5; par_list[kYdiv] = 4;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    for(int draw = 0; draw<4; ++draw){
//      for(int l = 0; l<2; ++l){
//	// TDC U
//	base_id = HistMaker::getUniqueID(kSP0, kSP0_L1+l+2*draw, kTDC);
//	for(int i = 0; i<NumOfSegSP0; ++i){id_list.push_back(base_id + i);}
//
//	// TDC D
//	base_id = HistMaker::getUniqueID(kSP0, kSP0_L1+l+2*draw, kTDC, NumOfSegSP0+1);
//	for(int i = 0; i<NumOfSegSP0; ++i){id_list.push_back(base_id + i);}
//      }
//      drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//    }
//  }

  // SDC3 ----------------------------------------------------------------
  if(name == CONV_STRING(kSDC3)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // SDC3 TDC
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 1500;
    base_id = HistMaker::getUniqueID(kSDC3, 0, kTDC);
    for(int i = 0; i<NumOfLayersSDC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // SDC3 HitPat
    base_id = HistMaker::getUniqueID(kSDC3, 0, kHitPat);
    for(int i = 0; i<NumOfLayersSDC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // SDC3 Multi
    base_id = HistMaker::getUniqueID(kSDC3, 0, kMulti);
    for(int i = 0; i<NumOfLayersSDC3; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

  // SDC4 ----------------------------------------------------------------
  if(name == CONV_STRING(kSDC4)){
    // For all
    int base_id = 0;
    par_list[kXdiv] = 2; par_list[kYdiv] = 3;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);

    // SDC4 TDC
    par_list[kXrange_min] = 700; par_list[kXrange_max] = 1600;
    base_id = HistMaker::getUniqueID(kSDC4, 0, kTDC);
    for(int i = 0; i<NumOfLayersSDC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    // SDC4 HitPat
    base_id = HistMaker::getUniqueID(kSDC4, 0, kHitPat);
    for(int i = 0; i<NumOfLayersSDC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);

    // SDC4 Multi
    base_id = HistMaker::getUniqueID(kSDC4, 0, kMulti);
    for(int i = 0; i<NumOfLayersSDC4; ++i){id_list.push_back(base_id + i);}
    drawOneCanvas(id_list, par_list, false, false);
  }

  // TOF ----------------------------------------------------------------
  if(name == CONV_STRING(kTOF)){
    int base_id = 0;
    int index   = 0;
    // ADC
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    // ADC U
    base_id = HistMaker::getUniqueID(kTOF, 0, kADC);
    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // ADC D
    index = 0;
    base_id = HistMaker::getUniqueID(kTOF, 0, kADC, NumOfSegTOF+1);
    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC U
    index = 0;
    base_id = HistMaker::getUniqueID(kTOF, 0, kTDC);
    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC D
    index = 0;
    base_id = HistMaker::getUniqueID(kTOF, 0, kTDC, NumOfSegTOF+1);
    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kTOF, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kTOF, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // TOF_HT ----------------------------------------------------------------
  if(name == CONV_STRING(kTOF_HT)){
    int base_id = 0;
    int index   = 0;
//    // ADC
//    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC U
//    base_id = HistMaker::getUniqueID(kTOF_HT, 0, kADC);
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // ADC D
//    index = 0;
//    base_id = HistMaker::getUniqueID(kTOF_HT, 0, kADC, NumOfSegTOF_HT+1);
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC U
    index = 0;
    base_id = HistMaker::getUniqueID(kTOF_HT, 0, kTDC);
    for(int i = 0; i<NumOfSegTOF_HT; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

//    // TDC D
//    index = 0;
//    base_id = HistMaker::getUniqueID(kTOF_HT, 0, kTDC, NumOfSegTOF_HT+1);
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegTOF_HT/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kTOF_HT, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kTOF_HT, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // TOFMT --------------------------------------------------------------
  if(name == CONV_STRING(kTOFMT)){
    int base_id = 0;
    int index   = 0;

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC MT
    index = 0;
    base_id = HistMaker::getUniqueID(kTOFMT, 0, kTDC);
    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegTOF/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kTOFMT, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kTOFMT, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // LAC ----------------------------------------------------------------
  if(name == CONV_STRING(kLAC)){
    int base_id = 0;
    int index   = 0;
    // ADC
    par_list[kXdiv] = 5; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    index = 0;
    base_id = HistMaker::getUniqueID(kLAC, 0, kADC);
    for(int i = 0; i<NumOfSegLAC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegLAC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 5; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    index = 0;
    base_id = HistMaker::getUniqueID(kLAC, 0, kTDC);
    for(int i = 0; i<NumOfSegLAC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegLAC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kLAC, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kLAC, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // LC ----------------------------------------------------------------
  if(name == CONV_STRING(kLC)){
    int base_id = 0;
    int index   = 0;
//    // ADC
//    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
//    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2000;
//    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
//    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);
//
//    // ADC U
//    base_id = HistMaker::getUniqueID(kLC, 0, kADC);
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    // ADC D
//    index = 0;
//    base_id = HistMaker::getUniqueID(kLC, 0, kADC, NumOfSegLC+1);
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 4000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC U
    index = 0;
    base_id = HistMaker::getUniqueID(kLC, 0, kTDC);
    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

//    // TDC D
//    index = 0;
//    base_id = HistMaker::getUniqueID(kLC, 0, kTDC, NumOfSegLC+1);
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
//
//    for(int i = 0; i<NumOfSegLC/2; ++i){id_list.push_back(base_id + index++);}
//    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kLC, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kLC, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // MsT ----------------------------------------------------------------
  if(name == CONV_STRING(kMsT)){
    int base_id = 0;
    int index   = 0;

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 2048;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    // TDC U
    index = 0;
    base_id = HistMaker::getUniqueID(kMsT, 0, kTDC);
    for(int i = 0; i<NumOfSegMsT/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
    for(int i = 0; i<NumOfSegMsT/2; ++i){id_list.push_back(base_id + index++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and Flag
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    base_id = HistMaker::getUniqueID(kMsT, 0, kHitPat, 0);
    id_list.push_back(base_id);
    id_list.push_back(base_id+1);
    id_list.push_back(HistMaker::getUniqueID(kMsT, 0, kTDC2D));
    drawOneCanvas(id_list, par_list, false, false, "colz");
    par_list[kXdiv] = 1; par_list[kYdiv] = 1;
    id_list.push_back(HistMaker::getUniqueID(kMsT, 0, kHitPat2D, 0));
    drawOneCanvas(id_list, par_list, false, false, "text");
  }

  // PWO ----------------------------------------------------------------
  if(name == CONV_STRING(kPWO)){
    const int box_id = 3;
    int base_id = 0;
    // ADC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x800;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyADC)  | GuiPs::isOptOn(kExpDataSheet);

    base_id = HistMaker::getUniqueID(kPWO, 0, kADC);
    for(int i=0; i<NumOfUnitPWO[box_id]; ++i){ id_list.push_back(base_id++);}
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // TDC
    par_list[kXdiv] = 4; par_list[kYdiv] = 4;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;
    flag_xaxis = GuiPs::isOptOn(kFixXaxis) | GuiPs::isOptOn(kExpDataSheet);
    flag_log   = GuiPs::isOptOn(kLogyTDC)  | GuiPs::isOptOn(kExpDataSheet);

    base_id = HistMaker::getUniqueID(kPWO, 0, kTDC);
    for(int i=0; i<NumOfUnitPWO[box_id]; ++i){ id_list.push_back(base_id++); }
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);

    // HitPat and multiplicity
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;

    id_list.push_back(HistMaker::getUniqueID(kPWO, 0, kHitPat));
    id_list.push_back(HistMaker::getUniqueID(kPWO, 0, kMulti));
    drawOneCanvas(id_list, par_list, false, false);
  }

  // Correlation ------------------------------------------------------------
  if(name == CONV_STRING(kCorrelation)){
    par_list[kXdiv] = 2; par_list[kYdiv] = 2;
    id_list.push_back(HistMaker::getUniqueID(kCorrelation, 0, 0, 1));
    id_list.push_back(HistMaker::getUniqueID(kCorrelation, 0, 0, 2));
    id_list.push_back(HistMaker::getUniqueID(kCorrelation, 0, 0, 3));
    id_list.push_back(HistMaker::getUniqueID(kCorrelation, 0, 0, 4));
    drawOneCanvas(id_list, par_list, false, false, "box");
  }

  // TriggerFlag------------------------------------------------------------
  if(name == CONV_STRING(kTriggerFlag)){
    int base_id = 0;
    par_list[kXdiv] = 6; par_list[kYdiv] = 5;
    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;

    base_id = HistMaker::getUniqueID( kTriggerFlag, 0, kTDC );
    for( int i=0; i<30; ++i){ id_list.push_back(base_id++); }
    drawOneCanvas(id_list, par_list, flag_xaxis, false);

    par_list[kXdiv] = 1; par_list[kYdiv] = 1;
    //    par_list[kXrange_min] = 0; par_list[kXrange_max] = 0x1000;

    base_id = HistMaker::getUniqueID(kTriggerFlag, 0, kHitPat);
    id_list.push_back(base_id);
    // base_id = HistMaker::getUniqueID(kTriggerFlag, 1, kHitPat);
    // id_list.push_back(base_id); //for(E07)
    drawOneCanvas(id_list, par_list, false, false);
  }

  //// EMC ----------------------------------------------------------------
  //if(name == CONV_STRING(kEMC)){
  //  int base_id = 0;
  //  // XYpos
  //  par_list[kXdiv] = 1; par_list[kYdiv] = 1;
  //  base_id = HistMaker::getUniqueID(kEMC, 0, kXYpos);
  //  for(int i = 0; i<NumOfSegEMC; ++i){ id_list.push_back( base_id++ ); }
  //  drawOneCanvas(id_list, par_list, flag_xaxis, flag_log, "box" );
  //}

  // TimeStamp ----------------------------------------------------------------
  if(name == CONV_STRING(kTimeStamp)){
    int base_id = 0;
    par_list[kXdiv] = 4; par_list[kYdiv] = 3;
    base_id = HistMaker::getUniqueID(kTimeStamp, 0, kTDC);
    for(int i=0; i<NumOfPlaneVmeRm; ++i ){ id_list.push_back( base_id++ ); }
    drawOneCanvas(id_list, par_list, flag_xaxis, flag_log);
  }

  // DC Eff ----------------------------------------------------------------
  if(name == CONV_STRING(kDCEff)){
    drawDCEff();
  }
}

// -------------------------------------------------------------------------
// drawOneChanvas
// -------------------------------------------------------------------------
void PsMaker::drawOneCanvas(std::vector<int>& id_list,
			    std::vector<int>& par_list,
			    bool flag_xaxis,
			    bool flag_log,
			    const char* optDraw
			    )
{
  static const std::string MyFunc = "drawOneCanvas ";
  // error check
  int NofPad  = par_list[kXdiv]*par_list[kYdiv];
  int NofHist = id_list.size();
  if(NofHist > NofPad){
    std::cerr << "#E: " << MyName << MyFunc
	      << "Num of pad is not enough " << NofHist << "/" << NofPad << std::endl
	      << " 1st unique ID: " << id_list[0] << std::endl;
    return;
  }

  // make new ps page
  ps_->NewPage();
  cps_->Divide(par_list[kXdiv], par_list[kYdiv]);
  for(int i = 0; i<(int)id_list.size(); ++i){
    if(id_list[i]==-1) continue;
    cps_->cd(i + 1)->SetGrid();
    if(flag_log){
      // log scale flag
      cps_->GetPad(i + 1)->SetLogy(1);
    }

    TH1* h = GHist::get(id_list[i]);
    if(!h){
      std::cerr << "#E: " << MyName << MyFunc
		<< "Pointer is NULL\n"
		<< " Unique ID: " << id_list[i] << std::endl;
      return;
    }

    if(flag_xaxis){
      // xaxis fix flag
      h->GetXaxis()->SetRangeUser(par_list[kXrange_min],
				  par_list[kXrange_max]);
    }

    // Rebin
    TString hclass = h->ClassName();
    TString hname  = Form("hclone_%d", id_list[i]);
    if( hclass.Contains("TH2") &&
	h->GetNbinsX() * h->GetNbinsY() > 200000 ){
      if( gROOT->FindObject(hname) )
	gROOT->FindObject(hname)->Delete();
      TH2 *hclone = (TH2*)h->Clone(hname);
      if( h->GetNbinsX() > 200 )
	hclone->RebinX(h->GetNbinsX()/200);
      if( h->GetNbinsY() > 200 )
	hclone->RebinY(h->GetNbinsY()/200);
      hclone->SetLineColor(1);
      hclone->Draw(optDraw);
    } else {
      h->SetLineColor(1);
      h->Draw(optDraw);
    }
  }

  cps_->Update();
  cps_->cd();
  clearOneCanvas(par_list[kXdiv]*par_list[kYdiv]);

  id_list.clear();
}

// -------------------------------------------------------------------------
// clearOneCanvas
// -------------------------------------------------------------------------
void PsMaker::clearOneCanvas(int npad)
{
  for(int i = 0; i<npad; ++i){
    cps_->GetPad(i + 1)->Close();
  }
  return;
}

// -------------------------------------------------------------------------
// drawRunNumber
// -------------------------------------------------------------------------
void PsMaker::drawRunNumber()
{
  hddaq::unpacker::UnpackerManager& g_unpacker =
    hddaq::unpacker::GUnpacker::get_instance();
  int runno = g_unpacker.get_root()->get_run_number();

  ps_->NewPage();
  cps_->Divide(1,1);
  cps_->cd(1);
  cps_->GetPad(1)->Range(0,0,100,100);

  TText text;
  text.SetTextSize(0.2);
  text.SetTextAlign(22);
  text.DrawText(50.,65., Form("Run# %d", runno));
  TTimeStamp stamp;
  stamp.Add( -stamp.GetZoneOffset() );
  text.SetTextSize(0.1);
  text.DrawText(50.,40., stamp.AsString("s") );

  cps_->Update();
  cps_->cd();
  cps_->GetPad(1)->Close();
}

// -------------------------------------------------------------------------
// drawDCEff
// -------------------------------------------------------------------------
void PsMaker::drawDCEff( void )
{
  static const std::string MyFunc = "drawDCEff ";

  TLatex text;
  text.SetNDC();
  text.SetTextSize( 0.08 );

  // BcOut
  {
    ps_->NewPage();
    int x = 4; int y = 3;
    cps_->Divide( x, y );
    int bc3_id = HistMaker::getUniqueID( kBC3, 0, kMulti );
    for( int i=0; i<NumOfLayersBC3; ++i ){
      cps_->cd( i+1 );
      TH1* h = GHist::get( bc3_id+i+NumOfLayersBC3);
      if( !h ){
	std::cerr << "#E: " << MyName << MyFunc
		  << "Pointer is NULL\n"
		  << " Unique ID: " << bc3_id+i+NumOfLayersBC3 << std::endl;
	return;
      }
      h->SetLineColor(1);
      h->Draw();
      double nof0   = h->GetBinContent(1);
      double nofall = h->GetEntries();
      double eff    = 1. - nof0/nofall;
      double xpos   = 0.4;
      double ypos   = 0.5;
      text.SetTextSize( 0.08 );
      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
    }
    int bc4_id = HistMaker::getUniqueID( kBC4, 0, kMulti );
    for( int i=0; i<NumOfLayersBC4; ++i ){
      cps_->cd( i+1+NumOfLayersBC3 );
      TH1* h = GHist::get( bc4_id+i+NumOfLayersBC4 );
      if( !h ){
	std::cerr << "#E: " << MyName << MyFunc
		  << "Pointer is NULL\n"
		  << " Unique ID: " << bc4_id+i+NumOfLayersBC4 << std::endl;
	return;
      }
      h->SetLineColor(1);
      h->Draw();
      double nof0   = h->GetBinContent(1);
      double nofall = h->GetEntries();
      double eff    = 1. - nof0/nofall;
      double xpos   = 0.4;
      double ypos   = 0.5;
      text.SetTextSize( 0.08 );
      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
    }
    cps_->Update();
    cps_->cd();
    clearOneCanvas( x*y );
  }

  // SdcIn
  {
    ps_->NewPage();
    int x = 4; int y = 2;
    cps_->Divide( x, y );
    int sdc1_id = HistMaker::getUniqueID( kSDC1, 0, kMulti );
    for( int i=0; i<NumOfLayersSDC1; ++i ){
      cps_->cd( i+1 );
      TH1* h = GHist::get( sdc1_id+i+NumOfLayersSDC1 );
      if( !h ){
	std::cerr << "#E: " << MyName << MyFunc
		  << "Pointer is NULL\n"
		  << " Unique ID: " << sdc1_id+i+NumOfLayersSDC1 << std::endl;
	return;
      }
      h->SetLineColor(1);
      h->Draw();
      double nof0   = h->GetBinContent(1);
      double nofall = h->GetEntries();
      double eff    = 1. - nof0/nofall;
      double xpos   = 0.35;
      double ypos   = 0.5;
      text.SetTextSize( 0.08 );
      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
    }
    cps_->Update();
    cps_->cd();
    clearOneCanvas( x*y );
  }

  // SdcOut
  {
    ps_->NewPage();
    int x = 4; int y = 2;
    cps_->Divide( x, y );
    int sdc2_id = HistMaker::getUniqueID( kSDC2, 0, kMulti );
    for( int i=0; i<NumOfLayersSDC2; ++i ){
      cps_->cd( i+1 );
      TH1* h = GHist::get( sdc2_id+i+NumOfLayersSDC2 );
      if( !h ){
	std::cerr << "#E: " << MyName << MyFunc
		  << "Pointer is NULL\n"
		  << " Unique ID: " << sdc2_id+i+NumOfLayersSDC2 << std::endl;
	return;
      }
      h->SetLineColor(1);
      h->Draw();
      double nof0   = h->GetBinContent(1);
      double nofall = h->GetEntries();
      double eff    = 1. - nof0/nofall;
      double xpos   = 0.35;
      double ypos   = 0.5;
      text.SetTextSize( 0.08 );
      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
    }
    int sdc3_id = HistMaker::getUniqueID( kSDC3, 0, kMulti );
    for( int i=0; i<NumOfLayersSDC3; ++i ){
      cps_->cd( i+1+NumOfLayersSDC2 );
      TH1* h = GHist::get( sdc3_id+i+NumOfLayersSDC3 );
      if( !h ){
	std::cerr << "#E: " << MyName << MyFunc
		  << "Pointer is NULL\n"
		  << " Unique ID: " << sdc3_id+i+NumOfLayersSDC3 << std::endl;
	return;
      }
      h->SetLineColor(1);
      h->Draw();
      double nof0   = h->GetBinContent(1);
      double nofall = h->GetEntries();
      double eff    = 1. - nof0/nofall;
      double xpos   = 0.35;
      double ypos   = 0.5;
      text.SetTextSize( 0.08 );
      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
    }
    cps_->Update();
    cps_->cd();
    clearOneCanvas( x*y );
  }

  // SSD
//  {
//    ps_->NewPage();
//    int x = 4; int y = 2;
//    cps_->Divide( x, y );
//    int ssd1_id = HistMaker::getUniqueID( kSSD1, 0, kMulti );
//    for( int i=0; i<NumOfLayersSSD1; ++i ){
//      cps_->cd( i+1 );
//      TH1* h = GHist::get( ssd1_id+i*2+1 );
//      if( !h ){
//	std::cerr << "#E: " << MyName << MyFunc
//		  << "Pointer is NULL\n"
//		  << " Unique ID: " << ssd1_id+i*2+1 << std::endl;
//	continue;
//      }
//      h->SetLineColor(1);
//      h->GetXaxis()->SetRangeUser(0,30);
//      h->Draw();
//      double nof0   = h->GetBinContent(1);
//      double nofall = h->GetEntries();
//      double eff    = 1. - nof0/nofall;
//      double xpos   = 0.35;
//      double ypos   = 0.5;
//      text.SetTextSize( 0.08 );
//      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
//    }
//    int ssd2_id = HistMaker::getUniqueID( kSSD2, 0, kMulti );
//    for( int i=0; i<NumOfLayersSSD2; ++i ){
//      cps_->cd( i+1+NumOfLayersSSD1 );
//      TH1* h = GHist::get( ssd2_id+i*2+1 );
//      if( !h ){
//	std::cerr << "#E: " << MyName << MyFunc
//		  << "Pointer is NULL\n"
//		  << " Unique ID: " << ssd2_id+i*2+1 << std::endl;
//	continue;
//      }
//      h->SetLineColor(1);
//      h->GetXaxis()->SetRangeUser(0,30);
//      h->Draw();
//      double nof0   = h->GetBinContent(1);
//      double nofall = h->GetEntries();
//      double eff    = 1. - nof0/nofall;
//      double xpos   = 0.35;
//      double ypos   = 0.5;
//      text.SetTextSize( 0.08 );
//      text.DrawLatex( xpos, ypos, Form("plane eff. %.2f", eff ) );
//    }
//    cps_->Update();
//    cps_->cd();
//    clearOneCanvas( x*y );
//  }

}
