/**
 *  file: HodoAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "Hodo2Hit.hh"
#include "BHTHit.hh"

class RawData;
class Hodo1Hit;
class Hodo2Hit;
class BHTHit;

typedef std::vector <Hodo1Hit*> Hodo1HitContainer;
typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <BHTHit*> BHTHitContainer;

//______________________________________________________________________________
class HodoAnalyzer
{
public:
  HodoAnalyzer( void );
  ~HodoAnalyzer( void );

  static HodoAnalyzer& GetInstance( void );

private:
  HodoAnalyzer( const HodoAnalyzer& );
  HodoAnalyzer& operator =( const HodoAnalyzer& );

private:
#if T98
  BHTHitContainer         m_BHDCont;
#elif E73_2024
  BHTHitContainer         m_BHDCont;
#elif E72
  BHTHitContainer         m_BHDCont;
#else
  Hodo2HitContainer     m_BHDCont;
#endif
  Hodo2HitContainer     m_T0Cont;
  Hodo2HitContainer     m_T0newCont;
  Hodo2HitContainer     m_E0Cont;
  Hodo2HitContainer     m_DEFCont;
#if T98
  Hodo2HitContainer       m_VetoCont;
  Hodo2HitContainer       m_BTCCont;
  Hodo2HitContainer       m_RCCont;
#elif E73_2024
  Hodo2HitContainer       m_VetoCont;
  Hodo2HitContainer       m_BTCCont;
#endif

public:
  bool DecodeRawHits( RawData* rawData );
  bool DecodeHodoHits(const int &detid, Hodo2HitContainer &m_Cont, RawData *rawData );
  bool DecodeBHTHits(const int &detid, BHTHitContainer &m_Cont, RawData *rawData );

  inline int  GetNHits( int detID )  const;
  int  GetNHitsBHD( void )  const { return m_BHDCont.size();  };
  int  GetNHitsT0( void )  const { return m_T0Cont.size();  };
  int  GetNHitsT0new( void )  const { return m_T0newCont.size();  };

  inline Hodo2Hit * GetHit( int detID, std::size_t i )  const;
  inline Hodo2Hit * GetHitBHD( std::size_t i )  const;
  inline Hodo2Hit * GetHitT0( std::size_t i )  const;
  inline Hodo2Hit * GetHitT0new( std::size_t i )  const;
  
  bool ReCalcBHDHits( bool applyRecursively=false );
  bool ReCalcT0Hits( bool applyRecursively=false );
  bool ReCalcT0newHits( bool applyRecursively=false );
  bool ReCalcAll( void );

  void TimeCutBHD(double tmin, double tmax);
  void TimeCutT0(double tmin, double tmax);
  void TimeCutT0new(double tmin, double tmax);

private:
  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, double tmin, double tmax);

};

inline int
HodoAnalyzer::GetNHits( int detID )  const
{
  switch(detID){
  case DetIdBHD:
    return m_BHDCont.size();
  case DetIdT0:
    return m_T0Cont.size();
  case DetIdT0new:
    return m_T0newCont.size();
  case DetIdE0:
    return m_E0Cont.size();
  case DetIdDEF:
    return m_DEFCont.size();
#if T98
  case DetIdVeto:
    return m_VetoCont.size();
  case DetIdBTC:
    return m_BTCCont.size();
  case DetIdRC:
    return m_RCCont.size();
#endif
  default: 
    return 0;
  }
}
//______________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHit( int detID, std::size_t i ) const
  {
  switch(detID){
  case DetIdBHD:
  if( i<m_BHDCont.size() )
    return m_BHDCont[i];
  else 
      return 0;
  case DetIdT0:
    if( i<m_T0Cont.size() )
      return m_T0Cont[i];
    else 
      return 0;
  case DetIdT0new:
    if( i<m_T0newCont.size() )
      return m_T0newCont[i];
    else 
      return 0;
  case DetIdE0:
    if( i<m_E0Cont.size() )
      return m_E0Cont[i];
    else 
      return 0;
  case DetIdDEF:
    if( i<m_DEFCont.size() )
      return m_DEFCont[i];
    else 
      return 0;
#if T98
  case DetIdRC:
    if( i<m_RCCont.size() )
      return m_RCCont[i];
    else 
      return 0;
  case DetIdVeto:
    if( i<m_VetoCont.size() )
      return m_VetoCont[i];
    else 
      return 0;
  case DetIdBTC:
    if( i<m_BTCCont.size() )
      return m_BTCCont[i];
    else 
      return 0;
#endif
  default: 
    return 0;
  }
}
//______________________________________________________________________________
#endif
