// -*- C++ -*-

#include <cmath>

#include <TMath.h>

#include "DCGeomRecord.hh"
#include "FuncName.hh"

ClassImp(DCGeomRecord);

//______________________________________________________________________________
DCGeomRecord::DCGeomRecord( Int_t id, const TString& name,
			    Double_t x, Double_t y, Double_t z, Double_t ta,
			    Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
			    Double_t w0, Double_t dd, Double_t ofs )
  : TObject(),
    id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
    rotAngle1_(ra1), rotAngle2_(ra2),
    length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
{
  calcVectors();
}

//______________________________________________________________________________
DCGeomRecord::DCGeomRecord( Int_t id, const TString& name,
			    const ThreeVector pos, Double_t ta,
			    Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
			    Double_t w0, Double_t dd, Double_t ofs )
  : TObject(),
    id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
    rotAngle1_(ra1), rotAngle2_(ra2),
    length_(length), resol_(resol), w0_(w0), dd_(dd), ofs_(ofs)
{
  calcVectors();
}

//______________________________________________________________________________
DCGeomRecord::~DCGeomRecord( void )
{
}

//______________________________________________________________________________
void
DCGeomRecord::calcVectors( void )
{
  Double_t ct1=cos(rotAngle1_*TMath::DegToRad()), st1=sin(rotAngle1_*TMath::DegToRad());
  Double_t ct2=cos(rotAngle2_*TMath::DegToRad()), st2=sin(rotAngle2_*TMath::DegToRad());
  Double_t ct0=cos(tiltAngle_*TMath::DegToRad()), st0=sin(tiltAngle_*TMath::DegToRad());

  dxds_ =  ct0*ct2-st0*ct1*st2;
  dxdt_ = -st0*ct2-ct0*ct1*st2;
  dxdu_ =  st1*st2;

  dyds_ =  ct0*st2+st0*ct1*ct2;
  dydt_ = -st0*st2+ct0*ct1*ct2;
  dydu_ = -st1*ct2;

  dzds_ =  st0*st1;
  dzdt_ =  ct0*st1;
  dzdu_ =  ct1;


  dsdx_ =  ct0*ct2-st0*ct1*st2;
  dsdy_ =  ct0*st2+st0*ct1*ct2;
  dsdz_ =  st0*st1;

  dtdx_ = -st0*ct2-ct0*ct1*st2;
  dtdy_ = -st0*st2+ct0*ct1*ct2;
  dtdz_ =  ct0*st1;

  dudx_ =  st1*st2;
  dudy_ = -st1*ct2;
  dudz_ =  ct1;
}

//______________________________________________________________________________
Int_t
DCGeomRecord::WireNumber( Double_t pos ) const
{
  Double_t dw = ((pos-ofs_)/dd_)+w0_;
  Int_t    iw = (Int_t)(dw);
  if( (dw-(Double_t)(iw))>0.5 )
    return iw+1;
  else
    return iw;
}
