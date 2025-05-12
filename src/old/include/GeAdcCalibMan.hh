// -*- C++ -*-

#ifndef GE_ADC_CALIB_MAN_HH
#define GE_ADC_CALIB_MAN_HH

#include <map>
#include <string>

#include <TObject.h>
#include <TString.h>

struct GeAdcCalibMap;
class GeHit;

//______________________________________________________________________________
class GeAdcCalibMan : public TObject
{
public:
  ~GeAdcCalibMan( void );
  static GeAdcCalibMan& GetInstance( void );

private:
  GeAdcCalibMan( void );
  GeAdcCalibMan( const GeAdcCalibMan & );
  GeAdcCalibMan & operator = ( const GeAdcCalibMan & );

private:
  TString ParamFileName_;
  mutable std::map <unsigned int, GeAdcCalibMap *> Cont_;

public:
  bool Initialize( const TString& filename );
  double CalibAdc( int SegId, int adc, double &adccal );

private:
  GeAdcCalibMap * getMap( int SegmentId ) const;
  void clearElements( void );

  ClassDef(GeAdcCalibMan,0);
};

//______________________________________________________________________________
inline GeAdcCalibMan&
GeAdcCalibMan::GetInstance( void )
{
  static GeAdcCalibMan g_instance;
  return g_instance;
}

#endif
