// -*- C++ -*-

#ifndef EVENT_ANALYZER_HH
#define EVENT_ANALYZER_HH

#include <TObject.h>

#include "BH2Filter.hh"

class RawData;
class DCAnalyzer;
class DCLocalTrack;
class HodoAnalyzer;
class BH2Hit;

namespace analyzer
{

//______________________________________________________________________________
class EventAnalyzer : public TObject
{
public:
  EventAnalyzer( void );
  ~EventAnalyzer( void );

protected:
  Bool_t                m_bh2filter_is_applied;
  RawData*              m_raw_data;
  DCAnalyzer*           m_dc_analyzer;
  HodoAnalyzer*         m_hodo_analyzer;
  BH2Filter::FilterList m_bh2filter_list;

public:
  void          ApplyBH2Filter( void );
  void          DecodeRawData( void );
  void          DecodeDCAnalyzer( void );
  void          DecodeHodoAnalyzer( void );
  void          TimeCutBFT( void );

  RawData*      GetRawData( void ) const { return m_raw_data; }
  DCAnalyzer*   GetDCAnalyzer( void ) const { return m_dc_analyzer; }
  HodoAnalyzer* GetHodoAnalyzer( void ) const { return m_hodo_analyzer; }

  Bool_t        TrackSearchBcOut( void );
  Bool_t        TrackSearchSdcIn( void );
  Bool_t        TrackSearchSdcOut( void );
  Bool_t        TrackSearchKurama( void );

  ClassDef(EventAnalyzer,0);
};

}

#endif
