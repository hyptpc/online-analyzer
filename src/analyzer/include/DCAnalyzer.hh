/**
 *  file: DCAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

class DCHit;
class DCLocalTrack;
class RawData;
class HodoCluster;

class Hodo1Hit;
class Hodo2Hit;
class HodoAnalyzer;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<DCLocalTrack*> DCLocalTrackContainer;

//______________________________________________________________________________
class DCAnalyzer
{
public:
  DCAnalyzer( void );
  ~DCAnalyzer( void );

private:
  DCAnalyzer( const DCAnalyzer& );
  DCAnalyzer& operator =( const DCAnalyzer& );

private:
  enum e_type { k_BLC1a, k_BLC1b, k_BLC2a, k_BLC2b, k_BPC1, k_BPC2, k_BcIn, k_BcOut,
		n_type };
  std::vector<bool>     m_is_decoded;
  std::vector<int>      m_much_combi;
  std::vector<DCHitContainer>       m_BLC1aHC;
  std::vector<DCHitContainer>       m_BLC1bHC;
  std::vector<DCHitContainer>       m_BLC2aHC;
  std::vector<DCHitContainer>       m_BLC2bHC;
  std::vector<DCHitContainer>       m_BPC1HC;
  std::vector<DCHitContainer>       m_BPC2HC;
  std::vector<DCHitContainer>       m_BcOutHC;

  DCLocalTrackContainer m_BcOutTC;

public:
  bool DecodeRawHits( RawData* rawData );
  // bool DecodeFiberHits( FiberCluster* FiberCl, int layer );
  bool DecodeRawHits( RawData* rawData, e_type k_type,const int &detid );
  bool DecodeDCHits( RawData* rawData , const int &detid);
  bool DecodeBcOutHits( RawData* rawData );
  

  inline const DCHitContainer& GetDCHC( const int &detid, int layer ) const;

  bool ReCalcDCHits( std::vector<DCHitContainer>& cont,
		     bool applyRecursively=false );
  bool ReCalcDCHits( bool applyRecursively=false );

  bool ReCalcTrack( DCLocalTrackContainer& cont, bool applyRecursively=false );

  bool ReCalcAll( void );

  bool TrackSearchBcOut (void);
  int GetNtracksBcOut( void )  const { return m_BcOutTC.size(); }

protected:
  void ClearDCHits( void );
  void ClearDCHits( const int &detid );
  void ChiSqrCut( DCLocalTrackContainer& cont, double chisqr );
  void ClearBcOutHits( void );
  void ClearTracksBcOut( void );
};

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetDCHC( const int &detid,int layer ) const
{
  if( layer>8 ) layer=0;
  switch(detid){
  case DetIdBPC1:
    return m_BPC1HC[layer];
  case DetIdBPC2:
    return m_BPC2HC[layer];
  case DetIdBLC1a:
    return m_BLC1aHC[layer];
  case DetIdBLC1b:
    return m_BLC1bHC[layer];
  case DetIdBLC2a:
    return m_BLC2aHC[layer];
  case DetIdBLC2b:
    return m_BLC2bHC[layer];
  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
    exit(0);
  }
}
#endif
