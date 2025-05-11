// -*- C++ -*-

#ifndef MACRO_BUILDER_HH
#define MACRO_BUILDER_HH

#include "HistMaker.hh"

class TCanvas;
class TObject;
class TString;

namespace analyzer
{

//______________________________________________________________________________
namespace macro
{
TObject* Get( TString name );
}

//______________________________________________________________________________
// for HttpServer
namespace http
{
TCanvas*  BAC();
TCanvas*  KVC1();
TCanvas*  KVC2();
TCanvas*  SAC();
TCanvas*  BH2ADC(Int_t ch);
TCanvas*  BH2ADCU();
TCanvas*  BH2ADCD();
TCanvas*  BH2TDC(Int_t ch);
TCanvas*  BH2TDCU();
TCanvas*  BH2TDCD();
TCanvas*  BH2TDCMT();
TCanvas*  HTOFADC();
TCanvas*  HTOFTDCMT();
TCanvas*  HTOFTDC2();
TCanvas*  Multiplicity();
TCanvas*  HitPat(); 
TCanvas*  BLDCTDC( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCTOT( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCTDCvsTOT( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCXYProf( DetectorType det, std::string strDet, int offs);
TCanvas*  BLDCProf( DetectorType det, std::string strDet);
TCanvas*  BeamAxis();
TCanvas*  BLDCHitPat( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCResid( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCMulti( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BLDCEff( DetectorType det, std::string strDet);
TCanvas*  BLDCCorr(int ichm, std::string strDet);

TCanvas*  BLDCWIRE( DetectorType det, std::string strDet, int layer, int nwires, int nx, int ny );
TCanvas*  BHTTOT( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  BHTTDCvsTOT( DetectorType det, std::string strDet, int nlayers, int init=0, int nx=4, int ny=2 );
TCanvas*  QDC( DetectorType det, std::string strDet, int subDet, int begin, int nlayers, int nx, int ny, double xmin=0., double xmax=500. );
TCanvas*  PbG( DetectorType det, std::string strDet, int subDet, int begin, int nlayers, int nx, int ny, double xmin=0., double xmax=500. );
TCanvas*  SDDADC( DetectorType det, std::string strDet, int nlayers, double xmin=0., double xmax=500., int rebin=8 );
TCanvas*  SDD( DetectorType det, std::string strDet, int subDet, int dataType, int nlayers, int nx, int ny, double xmin=0., double xmax=500. );
TCanvas*  MHTDCTDC( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny );
TCanvas*  SDDMHTDC( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny );
TCanvas*  MHTDCMeanTime( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny );
TCanvas*  MHTDCHitPatMulti( DetectorType det, std::string strDet );
TCanvas*  Ana(int type);
TCanvas*  Check(int type, std::string strTrig);
TCanvas*  TOF(int type, std::string strTrig);
TCanvas*  TOF_Btrg(int type, std::string strTrig);
TCanvas*  TOF2(int type, std::string strTrig);
TCanvas*  TOF2D(int type, std::string strTrig);
TCanvas*  AC(int type, std::string strTrig);

void      UpdateBLDCEfficiency();
void      UpdateCounterEfficiency();
}

};

#endif
