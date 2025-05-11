/**
 *  file: HodoAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
//#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"
//#include "HodoCluster.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("HodoAnalyzer");
  const double MaxTimeDifBHD   =  2.0;
  const double MaxTimeDifT0    =  2.0;
  const double MaxTimeDifT0new =  2.0;
}

#define Cluster 0

//______________________________________________________________________________
HodoAnalyzer::HodoAnalyzer( void )
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer( void )
{
  del::ClearContainer( m_BHDCont );
  del::ClearContainer( m_T0Cont );
  del::ClearContainer( m_T0newCont );
  del::ClearContainer( m_E0Cont );
  del::ClearContainer( m_DEFCont );
#if T98
  del::ClearContainer( m_VetoCont );
  del::ClearContainer( m_BTCCont );
  del::ClearContainer( m_RCCont );
#elif E73_2024
  del::ClearContainer( m_VetoCont );
  del::ClearContainer( m_BTCCont );
#endif
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
//______________________________________________________________________________
bool
HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
#if T98
  DecodeBHTHits( DetIdBHT, m_BHDCont, rawData );
  DecodeHodoHits( DetIdT0new, m_T0newCont, rawData );
  DecodeHodoHits( DetIdE0  ,  m_E0Cont,  rawData );
  DecodeHodoHits( DetIdDEF ,  m_DEFCont, rawData );
  DecodeHodoHits( DetIdVeto,  m_VetoCont, rawData );
  DecodeHodoHits( DetIdBTC ,  m_BTCCont , rawData );
  DecodeHodoHits( DetIdRC  ,  m_RCCont  , rawData );
#elif E73_2024
  DecodeBHTHits( DetIdBHT, m_BHDCont, rawData );
  DecodeHodoHits( DetIdT0new, m_T0newCont, rawData );
  DecodeHodoHits( DetIdE0  ,  m_E0Cont,  rawData );
  DecodeHodoHits( DetIdDEF ,  m_DEFCont, rawData );
  DecodeHodoHits( DetIdVeto,  m_VetoCont, rawData );
  DecodeHodoHits( DetIdBTC ,  m_BTCCont , rawData );
#elif E72
  DecodeBHTHits( DetIdBHT, m_BHDCont, rawData );
#else
  DecodeHodoHits( DetIdBHD ,  m_BHDCont, rawData );
#endif
  // DecodeHodoHits( DetIdT0  ,  m_T0Cont,  rawData );
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeHodoHits(const int &detid, Hodo2HitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    Hodo2Hit *hp = new Hodo2Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}
//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBHTHits(const int &detid, BHTHitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    BHTHit *hp = new BHTHit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcBHDHits( bool applyRecursively )
{
  int n = m_BHDCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit = m_BHDCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}
//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcT0Hits( bool applyRecursively )
{
  int n = m_T0Cont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit = m_T0Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}
//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcT0newHits( bool applyRecursively )
{
  int n = m_T0newCont.size();
  for( int i=0; i<n; ++i ){
    Hodo2Hit *hit = m_T0newCont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcAll( void )
{
  ReCalcBHDHits();
  ReCalcT0Hits();
  ReCalcT0newHits();
  return true;
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutBHD( double tmin, double tmax )
{
  //  TimeCut( m_BHDClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutT0( double tmin, double tmax )
{
  //  TimeCut( m_T0ClCont, tmin, tmax );
}

//______________________________________________________________________________
void
HodoAnalyzer::TimeCutT0new( double tmin, double tmax )
{
  //  TimeCut( m_T0newClCont, tmin, tmax );
}

//______________________________________________________________________________
//Implementation of Time cut for the cluster container
template <typename TypeCluster>
void
HodoAnalyzer::TimeCut( std::vector<TypeCluster>& cont,
		       double tmin, double tmax )
{
  std::vector<TypeCluster> DeleteCand;
  std::vector<TypeCluster> ValidCand;
  std::size_t size = cont.size();
  for( std::size_t i=0; i<size; ++i ){
    double ctime = cont.at(i)->CMeanTime();
    if(tmin < ctime && ctime < tmax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }
  del::ClearContainer( DeleteCand );

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}
