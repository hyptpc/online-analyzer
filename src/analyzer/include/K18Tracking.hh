// -*- C++ -*-

#ifndef K18_TRACKING_HH
#define K18_TRACKING_HH

#include <string>
#include <vector>

class DCAnalyzer;
class TH1;

//_____________________________________________________________________________
// K18 momentum: BcIn (BLC1) and BcOut (BLC2) local tracks fitted together
// through the D5 transfer matrix, ported from UserK18Tracking.cc (offline).
// Histograms are made by HistMaker::createK18Tracking().
namespace k18
{
// call in process_begin after the DCGEO and TM parameters are initialized
// (ConfMan::initializeK18Tracking() passes the conf key "D5FIELD:")
void Initialize(const std::string& d5_file);
// call in process_event after TrackSearchBcIn() and TrackSearchBcOut()
void Fill(const DCAnalyzer& dcAna, std::vector<TH1*>& hptr_array);
// re-read the D5 field file (called at run change and by HttpServer::ResetAll)
void ReloadD5Field();
}

#endif
