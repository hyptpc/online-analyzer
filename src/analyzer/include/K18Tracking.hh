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
class TText;

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
// htof_segs = HTOF segments with U&D hits in the HTOF_TDC window: the T2
// counts (Beam_T2_Count_Pi/K) use those behind T2 as the real T2 arrivals
void Fill(const DCAnalyzer& dcAna, std::vector<TH1*>& hptr_array,
          Int_t particle, const std::vector<Int_t>& htof_segs);
// T2 box center x, HypTPC local [mm] (DCGEO T2), 0 without T2: the T2
// histogram x range is centered on it (HistMaker::createK18Tracking)
Double_t T2CenterX();
// T2 efficiency text (eff. = InBox&HTOF / HTOF) for kPi or kKaon, drawn by http::Beam
TText* T2CountText(Int_t particle);
// HTOF segments used for the T2 efficiency (0-based), drawn by http::Beam
TText* T2SegText(Int_t particle);
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
// h_target_x (optional): filled with the RK x at the target plane (DCGEO
// Target) of every drawn track, e.g. TPC_BeamProfile_RK
void UpdateEventDisplay(const DCAnalyzer& dcAna, TH1* h_target_x = nullptr);
// re-read the D5 field file (called at run change and by HttpServer::ResetAll)
void ReloadD5Field();
}

#endif
