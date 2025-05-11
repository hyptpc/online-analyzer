/**
 *  file: DCAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "DCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "ConfMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"

// Tracking routine selection __________________________________________________
namespace
{
  // using namespace ;
  const std::string& class_name("DCAnalyzer");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  //______________________________________________________________________________
}

//______________________________________________________________________________
DCAnalyzer::DCAnalyzer( void )
  : m_is_decoded(n_type),
    m_much_combi(n_type),
    m_BLC1aHC(8+1),
    m_BLC1bHC(8+1),
    m_BLC2aHC(8+1),
    m_BLC2bHC(8+1),
    m_BPC1HC(8+1),
    m_BPC2HC(8+1)
{
  for( int i=0; i<n_type; ++i ){
    m_is_decoded[i] = false;
    m_much_combi[i] = 0;
  }
  debug::ObjectCounter::increase(class_name);
}

DCAnalyzer::~DCAnalyzer( void )
{
  ClearDCHits(DetIdBLC1a);
  ClearDCHits(DetIdBLC1b);
  ClearDCHits(DetIdBLC2a);
  ClearDCHits(DetIdBLC2b);
  ClearDCHits(DetIdBPC1);
  ClearDCHits(DetIdBPC2);
  debug::ObjectCounter::decrease(class_name);
}
//______________________________________________________________________________
//______________________________________________________________________________
bool
DCAnalyzer::DecodeRawHits( RawData *rawData, e_type k_det, const int &detid )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_decoded[k_det] ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded" << std::endl;
    return true;
  }
  ClearDCHits(detid);
  for( int layer=1; layer<=8; ++layer ){
    const DCRHitContainer &RHitCont=rawData->GetDCRawHC(detid,layer-1);
    int nh = RHitCont.size();
    // if(detid==DetIdBPC)
    //    std::cout<<func_name<<"  "<<detid<<"  "<<nh<<std::endl;
    for( int i=0; i<nh; ++i ){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *hit   = new DCHit( detid, rhit->PlaneId(), rhit->WireId() );
      int       nhtdc = rhit->GetTdcSize();
      if(!hit) continue;
      for( int j=0; j<nhtdc; ++j ){
	hit->SetTdcVal( rhit->GetTdc(j) );
      }      
      if( hit->CalcDCObservables() ){
	switch(detid){
	case DetIdBPC1:
	    m_BPC1HC[layer].push_back(hit);
	    break;
	case DetIdBPC2:
	    m_BPC2HC[layer].push_back(hit);
	    break;
	case DetIdBLC1a:
	    m_BLC1aHC[layer].push_back(hit);
	    break;
	case DetIdBLC1b:
	    m_BLC1bHC[layer].push_back(hit);
	    break;
	case DetIdBLC2a:
	    m_BLC2aHC[layer].push_back(hit);
	    break;
	case DetIdBLC2b:
	    m_BLC2bHC[layer].push_back(hit);
	    break;
	default:
	  std::cout<<"E# invalid detector id "<< detid<<std::endl;
	  return false;
	}
      }
      else
	delete hit;
    }
  }
  m_is_decoded[k_det] = true;
  return true;
}


//______________________________________________________________________________
bool
DCAnalyzer::DecodeRawHits( RawData *rawData )
{
  ClearDCHits();
#if E73_2024
  DecodeRawHits( rawData, k_BPC1, DetIdBPC1 );
  DecodeRawHits( rawData, k_BPC2, DetIdBPC2 );
#endif
  DecodeRawHits( rawData, k_BLC1a, DetIdBLC1a );
  DecodeRawHits( rawData, k_BLC1b, DetIdBLC1b );
  DecodeRawHits( rawData, k_BLC2a, DetIdBLC2a );
  DecodeRawHits( rawData, k_BLC2b, DetIdBLC2b );
  return true;
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCHits( const int &detid )
{
  switch(detid){
  case DetIdBPC1:
    del::ClearContainerAll( m_BPC1HC );
  case DetIdBPC2:
    del::ClearContainerAll( m_BPC2HC );
  case DetIdBLC1a:
    del::ClearContainerAll( m_BLC1aHC );
  case DetIdBLC1b:
    del::ClearContainerAll( m_BLC1bHC );
  case DetIdBLC2a:
    del::ClearContainerAll( m_BLC2aHC );
  case DetIdBLC2b:
    del::ClearContainerAll( m_BLC2bHC );
  }
}

//______________________________________________________________________________
void
DCAnalyzer::ClearDCHits( void)
{
  del::ClearContainerAll( m_BPC1HC );
  del::ClearContainerAll( m_BPC2HC );
  del::ClearContainerAll( m_BLC1aHC );
  del::ClearContainerAll( m_BLC1bHC );
  del::ClearContainerAll( m_BLC2aHC );
  del::ClearContainerAll( m_BLC2bHC );
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcDCHits( std::vector<DCHitContainer>& cont,
			  bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t l=0; l<n; ++l ){
    const std::size_t m = cont[l].size();
    for( std::size_t i=0; i<m; ++i ){
      DCHit *hit = (cont[l])[i];
      if( !hit ) continue;
      hit->ReCalcDC(applyRecursively);
    }
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcDCHits( bool applyRecursively )
{
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcTrack( DCLocalTrackContainer& cont,
			 bool applyRecursively )
{
  const std::size_t n = cont.size();
  for( std::size_t i=0; i<n; ++i ){
    DCLocalTrack *track = cont[i];
    if( track ) track->ReCalc( applyRecursively );
  }
  return true;
}

//______________________________________________________________________________
bool
DCAnalyzer::ReCalcAll( void )
{
  ReCalcDCHits();
  return true;
}
//______________________________________________________________________________
void
DCAnalyzer::ChiSqrCut( DCLocalTrackContainer& TrackCont,
		       double chisqr )
{
  DCLocalTrackContainer DeleteCand;
  DCLocalTrackContainer ValidCand;
  int NofTrack = TrackCont.size();
  for(int i = NofTrack-1; i>=0; --i){
    DCLocalTrack* tempTrack = TrackCont.at(i);
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  del::ClearContainer( DeleteCand );

  TrackCont.clear();
  TrackCont.resize( ValidCand.size() );
  std::copy( ValidCand.begin(), ValidCand.end(), TrackCont.begin() );
  ValidCand.clear();
}
