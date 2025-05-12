// -*- C++ -*-

#include "BH2Filter.hh"
#include "DCAnalyzer.hh"
#include "Exception.hh"
#include "EventAnalyzer.hh"
#include "FuncName.hh"
#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

ClassImp(analyzer::EventAnalyzer);

namespace analyzer
{

namespace
{
  const UserParamMan& gUser = UserParamMan::GetInstance();
}

//______________________________________________________________________________
EventAnalyzer::EventAnalyzer( void )
  : TObject(),
    m_bh2filter_is_applied( false ),
    m_raw_data( new RawData ),
    m_dc_analyzer( new DCAnalyzer ),
    m_hodo_analyzer( new HodoAnalyzer )
{
}

//______________________________________________________________________________
EventAnalyzer::~EventAnalyzer( void )
{
  delete m_raw_data;
  delete m_dc_analyzer;
  delete m_hodo_analyzer;
}

//______________________________________________________________________________
void
EventAnalyzer::ApplyBH2Filter( void )
{
  BH2Filter& gBH2Filter = BH2Filter::GetInstance();
  gBH2Filter.Apply( *m_hodo_analyzer, *m_dc_analyzer, m_bh2filter_list );
  m_bh2filter_is_applied = true;
}

//______________________________________________________________________________
void
EventAnalyzer::DecodeRawData( void )
{
  try {
    m_raw_data->DecodeHits();
  } catch( const std::exception& e ){
    throw Exception( FUNC_NAME+" "+e.what() );
  }
}

//______________________________________________________________________________
void
EventAnalyzer::DecodeDCAnalyzer( void )
{
  try {
    m_dc_analyzer->DecodeRawHits( m_raw_data );
  } catch( const std::exception& e ){
    throw Exception( FUNC_NAME+" "+e.what() );
  }
}

//______________________________________________________________________________
void
EventAnalyzer::DecodeHodoAnalyzer( void )
{
  try {
    m_hodo_analyzer->DecodeRawHits( m_raw_data );
  } catch( const std::exception& e ){
    throw Exception( FUNC_NAME+" "+e.what() );
  }
}

//______________________________________________________________________________
void
EventAnalyzer::TimeCutBFT( void )
{
  static const Double_t MinTimeBFT = gUser.GetParameter("TimeBFT", 0);
  static const Double_t MaxTimeBFT = gUser.GetParameter("TimeBFT", 1);
  m_hodo_analyzer->TimeCutBFT(MinTimeBFT, MaxTimeBFT);
}

//______________________________________________________________________________
Bool_t
EventAnalyzer::TrackSearchBcOut( void )
{
  if( !m_dc_analyzer )
    return false;
  if( !m_bh2filter_is_applied )
    return m_dc_analyzer->TrackSearchBcOut();
  else
    return m_dc_analyzer->TrackSearchBcOut( m_bh2filter_list );
}

//______________________________________________________________________________
Bool_t
EventAnalyzer::TrackSearchSdcIn( void )
{
  if( !m_dc_analyzer )
    return false;

  return m_dc_analyzer->TrackSearchSdcIn();
}

//______________________________________________________________________________
Bool_t
EventAnalyzer::TrackSearchSdcOut( void )
{
  if( !m_dc_analyzer )
    return false;

  Hodo2HitContainer TOFCont;
  Int_t nhTof = m_hodo_analyzer->GetNHitsTOF();
  for( Int_t i=0; i<nhTof; ++i ){
    Hodo2Hit *hit = m_hodo_analyzer->GetHitTOF(i);
    if( !hit ) continue;
    TOFCont.push_back( hit );
  }

  // return m_dc_analyzer->TrackSearchSdcOut();
  return m_dc_analyzer->TrackSearchSdcOut( TOFCont );
}

//______________________________________________________________________________
Bool_t
EventAnalyzer::TrackSearchKurama( void )
{
  if( !m_dc_analyzer )
    return false;

  return m_dc_analyzer->TrackSearchKurama();
}

}
