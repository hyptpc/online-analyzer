// -*- C++ -*-

#ifndef EVENT_ANALYZER_HH
#define EVENT_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"

class DCAnalyzer;
class RawData;

//_____________________________________________________________________________
class EventAnalyzer
{
public:
  EventAnalyzer();
  ~EventAnalyzer();

private:

public:

  void DCRawHit(const TString& dcname,const RawData& rawData/*,
                beam::EBeamFlag beam_flag=beam::kAll*/);
  void DCHit(const TString& dcname, const DCAnalyzer& dcAna/*,
             beam::EBeamFlag beam_flag=beam::kAll*/);
  void BcInTracking(DCAnalyzer& dcAna/*, beam::EBeamFlag beam_flag=beam::kAll*/);
  void BcOutTracking(DCAnalyzer& dcAna/*, beam::EBeamFlag beam_flag=beam::kAll*/);

private:
};

#endif
