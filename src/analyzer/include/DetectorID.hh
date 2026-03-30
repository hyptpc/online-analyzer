// -*- C++ -*-
 
#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <TString.h>
#include <map>

const Int_t NumOfSegBHT = 63;
const Int_t NumOfSegBH2 = 15;
const Int_t NumOfSegBVH = 4;
const Int_t NumOfSegT1  = 1;
const Int_t NumOfSegT2  = 1;
const Int_t NumOfSegSCH = 64;

const Int_t NumOfSegHTOF = 34;
const Int_t NumOfSegKVC = 8;
const Int_t NumOfSegBAC = 1;

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
const int DetIdBcOut = 120;

// HypTPC
const Int_t DetIdTPC       = 70;
const Int_t NumOfLayersTPC = 32;
const Int_t NumOfTimeBucket = 170;
const Int_t NumOfAsadTPC    = 31;
const Int_t NumOfPadTPC     = 5768;

// 
const int DetIdVmeRm=91;
const int NumOfPlaneVmeRm=3;

//BcInOut
const Int_t PlMinBcIn        =   0;
const Int_t PlMaxBcIn        =  15;
const Int_t PlMinBcOut       =  16;
const Int_t PlMaxBcOut       =  31;
const Int_t PlOffsBc         = 100;
const Int_t NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;

inline const std::map<TString, std::vector<TString>> DCNameList =
{
  {"BcIn", { "BLC1a", "BLC1b" }},
  {"BcOut", { "BLC2a", "BLC2b" }},
};


namespace trigger
{
enum ETriggerFlag
  {
    // kSpillOnEnd,
    // kSpillOffEnd,
    // kL1SpillOn,
    // kL1SpillOff,
    kBHT, // kMtx2D1,
    kBH2OR,
    kBAC,
    kCVC,
    kSAC3,
    kSFV,
    kBH2,
    kHTOFFwd,
    kHTOFPiMtight,
    kKVC,
    kT1,
    kHTOFNIMMp2,
    kHTOFPiMwide,
    kTOF24,
    kCHANNEL14,
    kCHANNEL15,
    kTrigA_PS,
    kTrigB_PS,
    kTrigC_PS,
    kTrigD_PS,
    kTrigE_PS,
    kTrigF_PS,
    kTrigOR_A_PS,
    kTrigOR_B_PS,
    kClock_PS,
    NTriggerFlag,
  };

const std::vector<TString> STriggerFlag =
  {
    // "SpillOnEnd",
    // "SpillOffEnd",
    // "L1SpillOn",
    // "L1SpillOff",
    "BHT", // "Mtx2D1",
    "BH2OR",
    "BAC",
    "CVC",
    "SAC3",
    "SFV",
    "BH2",
    "HTOF-Fwd",
    "HTOF-PiMinus-Tight",
    "KVC",
    "T1",
    "HTOF-NIM-Mp",
    "HTOF-PiMinus-Wide",
    "TOF24",
    "CHANNEL14",
    "CHANNEL15",
    "TrigA-PS",
    "TrigB-PS",
    "TrigC-PS",
    "TrigD-PS",
    "TrigE-PS",
    "TrigF-PS",
    "TrigOR_A_PS",
    "TrigOR_B_PS",
    "Clock-PS",
  };
}

const Int_t NumOfSegTFlag = trigger::NTriggerFlag;

#endif
