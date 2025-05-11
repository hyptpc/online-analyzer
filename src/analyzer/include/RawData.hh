/**
 *  file: RawData.hh
 *  date: 2017.04.10
 *
 */

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include "DetectorID.hh"
#include <vector>

class HodoRawHit;
class DCRawHit;

typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;

//______________________________________________________________________________
class RawData
{
public:
  RawData( void );
  ~RawData( void );

private:
  RawData( const RawData& );
  RawData& operator=( const RawData& );

private:
  bool              m_is_decoded;
  HodoRHitContainer m_BHDRawHC;
  HodoRHitContainer m_T0RawHC;
  HodoRHitContainer m_T1RawHC;
  HodoRHitContainer m_E0RawHC;
  HodoRHitContainer m_DEFRawHC;
  HodoRHitContainer m_CDHRawHC;
  HodoRHitContainer m_PbF2RawHC;
  HodoRHitContainer m_PbGRawHC;
  HodoRHitContainer m_VetoRawHC;
  HodoRHitContainer m_BTCRawHC;
  HodoRHitContainer m_RCRawHC;  
  HodoRHitContainer m_CVCRawHC;  
  HodoRHitContainer m_NCRawHC;  
  std::vector<DCRHitContainer> m_BLC1aRawHC;
  std::vector<DCRHitContainer> m_BLC1bRawHC;
  std::vector<DCRHitContainer> m_BLC2aRawHC;
  std::vector<DCRHitContainer> m_BLC2bRawHC;
  std::vector<DCRHitContainer> m_BPC1RawHC;
  std::vector<DCRHitContainer> m_BPC2RawHC;
  
  HodoRHitContainer m_VmeCalibRawHC;


public:
  void                     ClearAll( void );
  bool                     DecodeHits( void );
  bool                     DecodeCalibHits( void );
  const HodoRHitContainer& GetHodoRawHC( const int &detid ) const;
  const DCRHitContainer&   GetDCRawHC( const int &detid, int layer ) const;
  const HodoRHitContainer& GetVmeCalibRawHC( void ) const;

};

//______________________________________________________________________________
//______________________________________________________________________________
inline const HodoRHitContainer&
RawData::GetHodoRawHC( const int &detid ) const
{
  switch(detid){
  case DetIdBHD:
    return m_BHDRawHC;
  case DetIdT0:
    return m_T0RawHC;
  case DetIdT1:
    return m_T1RawHC;
  case DetIdE0:
    return m_E0RawHC;
  case DetIdDEF:
    return m_DEFRawHC;
  case DetIdCDH:
    return m_CDHRawHC;
  case DetIdPbF2:
    return m_PbF2RawHC;
  case DetIdPbG:
    return m_PbGRawHC;
  case DetIdVeto:
    return m_VetoRawHC;
  case DetIdRC:
    return m_RCRawHC;
  case DetIdBTC:
    return m_BTCRawHC;
  case DetIdCVC:
    return m_CVCRawHC;
  case DetIdNC:
    return m_NCRawHC;
  }
}

//______________________________________________________________________________
inline const DCRHitContainer&
RawData::GetDCRawHC( const int &detid, int layer ) const
{
  if( layer<0 || layer>8 ) layer = 0;
  switch(detid){
  case DetIdBPC1:
    return m_BPC1RawHC[layer];
  case DetIdBPC2:
    return m_BPC2RawHC[layer];
  case DetIdBLC1a:
    return m_BLC1aRawHC[layer];
  case DetIdBLC1b:
    return m_BLC1bRawHC[layer];
  case DetIdBLC2a:
    return m_BLC2aRawHC[layer];
  case DetIdBLC2b:
    return m_BLC2bRawHC[layer];
  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
  }
}

inline const HodoRHitContainer&
RawData::GetVmeCalibRawHC( void ) const
{
  return m_VmeCalibRawHC;
}

#endif
