// -*- C++ -*-

#ifndef K18_TRACKING_HH
#define K18_TRACKING_HH

#include <string>
#include <vector>

#include <TString.h>

#include "ThreeVector.hh"

class DCAnalyzer;
class TGraph;
class TH1;

//_____________________________________________________________________________
// K18 momentum: BcIn (BLC1) and BcOut (BLC2) local tracks fitted together
// through the D5 transfer matrix, ported from UserK18Tracking.cc (offline).
// Histograms are made by HistMaker::createK18Tracking().
namespace k18
{
// beam track: BcIn-BcOut pair fitted through the D5 matrix
struct BeamTrack
{
  Double_t p;       // [GeV/c] = p0*(1+delta)
  Double_t delta;
  Double_t chisqr;  // reduced
  Double_t vo[4];   // fitted state at VO in the BcOut frame {x, u, y, v}
};

// max beam tracks drawn on the 2D event display
const Int_t kMaxDisplayTracks = 5;

// call in process_begin after the DCGEO and TM parameters are initialized
// (ConfMan::initializeK18Tracking() passes the conf key "D5FIELD:")
void Initialize(const std::string& d5_file);
// call in process_event after TrackSearchBcIn() and TrackSearchBcOut()
// particle = kPi / kKaon (beam TOF), or kAll when neither: the momentum is
// filled to Beam_Momentum (all) and to Beam_Momentum_Pi / _K
void Fill(const DCAnalyzer& dcAna, std::vector<TH1*>& hptr_array,
          Int_t particle);
// fitted beam tracks, best chisqr first (exclusive BcIn-BcOut pairs)
std::vector<BeamTrack> FitBeamTracks(const DCAnalyzer& dcAna);
// RK trajectory through the HS field map from VO to z_end, HypTPC local [mm]
Bool_t TraceBeam(const BeamTrack& track, Double_t z_end,
                 std::vector<ThreeVector>& points);
// z of a DCGEO entry in the HypTPC local frame, NaN if missing
Double_t LocalZ(const TString& name);
// RK beam lines for the 2D event display (z vs x), drawn by http::SHS2D_EVT
std::vector<TGraph*>& EventDisplayGraphs();
// call in process_event after TrackSearchBcIn() and TrackSearchBcOut()
void UpdateEventDisplay(const DCAnalyzer& dcAna);
// re-read the D5 field file (called at run change and by HttpServer::ResetAll)
void ReloadD5Field();
}

#endif
