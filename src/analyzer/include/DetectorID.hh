// -*- C++ -*-
 
#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <TString.h>

// Counters ___________________________________________________________
const int DetIdCDH    =  0;
const int DetIdBHD    =  1;
const int DetIdBHT    =  1;
const int DetIdT0     =  2;
const int DetIdAC     =  3;
const int DetIdT0new  =  4;
const int DetIdT1     =  4;
const int DetIdE0     =  5;
const int DetIdDEF    =  6;
const int DetIdStart  =  11;
const int DetIdStop   =  12;
const int DetIdCVC    =  14;
const int DetIdSDD    =  20;
const int DetIdSDDVeto=  25;
const int DetIdPbG    =  29;
const int DetIdPbF2   =  30;
const int DetIdVeto1  =  31;
const int DetIdVeto   =  31;
const int DetIdFinger =  32;
const int DetIdNC     =  32;
const int DetIdBTC    =  33;
const int DetIdGRAMS40    =  40;
const int DetIdVFT    =  53;
const int NumOfSegBHD    =  11;
const int NumOfSegT0     =  5;
const int NumOfSegT0new  =  1;
const int DetIdScaler    =  91;
const int DetIdRC        =  98;
const int DetIdCNCtest        =  80;
const int DetIdTrigFlag  =  99;

// Chambers
const int DetIdCDC = 100;
const int DetIdBLC1a  = 101;
const int DetIdBLC1b  = 102;
const int DetIdBLC2a  = 103;
const int DetIdBLC2b  = 104;
const int DetIdBPC = 105;
const int DetIdBPC2 = 105;
const int DetIdBPC1 = 106;
const int DetIdSDC = 106;
const int DetIdFDC = 107;
const int DetIdBLC1= 111;
const int DetIdBLC2= 112;

// 
const int DetIdVmeRm=91;
const int NumOfPlaneVmeRm=3;

namespace trigger
{
enum ETriggerFlag
  {
    kSpillOnEnd,
    kSpillOffEnd,
    kL1SpillOn,
    kL1SpillOff,
    kBHT, // kMtx2D1,
    kMtx2D2,
    kMtx3D,
    kBeamA,
    kBeamB,
    kBeamC,
    kBeamD,
    kBeamE,
    kBeamF,
    kTrigA,
    kTrigB,
    kTrigC,
    kTrigD,
    kTrigE,
    kTrigF,
    NTriggerFlag,
  };

const std::vector<TString> STriggerFlag =
  {
    "SpillOnEnd",
    "SpillOffEnd",
    "L1SpillOn",
    "L1SpillOff",
    "BHT", // "Mtx2D1",
    "Mtx2D2",
    "Mtx3D",
    "BeamA",
    "BeamB",
    "BeamC",
    "BeamD",
    "BeamE",
    "BeamF",
    "TrigA",
    "TrigB",
    "TrigC",
    "TrigD",
    "TrigE",
    "TrigF",
  };
}

const Int_t NumOfSegTFlag = trigger::NTriggerFlag;

#endif
