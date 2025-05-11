// -*- C++ -*-

#ifndef DC_GEOM_RECORD_HH
#define DC_GEOM_RECORD_HH

#include <iostream>
#include <string>
#include <functional>

#include "ThreeVector.hh"
#include "TString.h"
//______________________________________________________________________________
class DCGeomRecord : public TObject
{
public:
  DCGeomRecord( Int_t id, const TString &name,
                Double_t x, Double_t y, Double_t z, Double_t ta,
                Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
		Double_t w0, Double_t dd, Double_t ofs );
  DCGeomRecord( Int_t id, const TString &name,
                const ThreeVector pos, Double_t ta,
                Double_t ra1, Double_t ra2, Double_t length, Double_t resol,
		Double_t w0, Double_t dd, Double_t ofs );
  ~DCGeomRecord( void );

  DCGeomRecord( const DCGeomRecord& );
  DCGeomRecord& operator=( const DCGeomRecord& );

private:
  Int_t       id_;
  TString     name_;
  ThreeVector pos_;
  Double_t    tiltAngle_;
  Double_t    rotAngle1_;
  Double_t    rotAngle2_;
  Double_t    length_;
  Double_t    resol_;
  Double_t    w0_;
  Double_t    dd_;
  Double_t    ofs_;

  Double_t dxds_, dxdt_, dxdu_;
  Double_t dyds_, dydt_, dydu_;
  Double_t dzds_, dzdt_, dzdu_;

  Double_t dsdx_, dsdy_, dsdz_;
  Double_t dtdx_, dtdy_, dtdz_;
  Double_t dudx_, dudy_, dudz_;

public:
  const ThreeVector & Position( void ) const { return pos_; }
  ThreeVector NormalVector( void ) const
  { return ThreeVector( dxdu_, dydu_, dzdu_ ); }
  ThreeVector UnitVector( void ) const
  { return ThreeVector( dxds_, dyds_, dzds_ ); }

  Double_t dsdx( void ) const { return dsdx_; }
  Double_t dsdy( void ) const { return dsdy_; }
  Double_t dsdz( void ) const { return dsdz_; }
  Double_t dtdx( void ) const { return dtdx_; }
  Double_t dtdy( void ) const { return dtdy_; }
  Double_t dtdz( void ) const { return dtdz_; }
  Double_t dudx( void ) const { return dudx_; }
  Double_t dudy( void ) const { return dudy_; }
  Double_t dudz( void ) const { return dudz_; }

  Double_t dxds( void ) const { return dxds_; }
  Double_t dxdt( void ) const { return dxdt_; }
  Double_t dxdu( void ) const { return dxdu_; }
  Double_t dyds( void ) const { return dyds_; }
  Double_t dydt( void ) const { return dydt_; }
  Double_t dydu( void ) const { return dydu_; }
  Double_t dzds( void ) const { return dzds_; }
  Double_t dzdt( void ) const { return dzdt_; }
  Double_t dzdu( void ) const { return dzdu_; }

  Double_t WirePos( Double_t wire ) const { return dd_*(wire - w0_)+ofs_; }
  Int_t    WireNumber( Double_t pos ) const;

private:
  void calcVectors( void );

  friend class DCGeomMan;
  friend struct DCGeomRecordComp;

  ClassDef(DCGeomRecord,0);
};

//______________________________________________________________________________
struct DCGeomRecordComp
  : public std::binary_function <DCGeomRecord *, DCGeomRecord *, Bool_t>
{
  Bool_t operator()( const DCGeomRecord * const p1,
		   const DCGeomRecord * const p2 ) const
  { return p1->id_ < p2->id_; }
};

#endif
