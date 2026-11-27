// -*- C++ -*-

#include "K18Tracking.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TGraph.h>
#include <TH1.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TText.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "DetSizeMan.hh"
#include "FieldMan.hh"
#include "HistMaker.hh"
#include "TransferMatrixMan.hh"
#include "ThreeVector.hh"

#include <UnpackerManager.hh>

namespace
{
auto& gTM = TransferMatrixMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gField = FieldMan::GetInstance();

const Int_t    kNParam      = 5;
const Int_t    kMatrixOrder = 2;  // with D5Matrix2nd corrections
// VI = BLC1a/b pair center, VO = BLC2a/b pair center in the local tracking
// frames [mm]. Taken from DCGEO in k18::Initialize() (offline geometry gives
// VI = 0, VO = -1299.815; the online geometry has the BcIn origin at BLC1a).
Double_t kZVI = 0.;
Double_t kZVO = -1299.815;
const Int_t    kMaxNTrack = 10;   // max local tracks per side in the fit
// D5 field for the central momentum, read from the file given by the conf
// key "D5FIELD:" (one value [T], '#' starts a comment). The file is read at
// start-up and again at every run change; it overrides the TM file value:
// p0 [GeV/c] = 0.299792458 * |B| [T] * R [m] * FieldFactor (as D5FieldMan)
std::string g_d5_file;
Double_t    g_d5_field  = 0.;    // [T], last value read
const Double_t kD5Radius    = 1.8984;    // [m]
const Double_t kFieldFactor = 1.08;

//_____________________________________________________________________________
// beam charge from the D5 field sign (-1 when D5 is off), as offline
// UserK18Tracking (UpdateCentralMomentum)
Int_t
BeamCharge()
{
  static const Double_t kMinD5Field = 0.01; // [T], below this D5 is off
  if(TMath::Abs(g_d5_field) < kMinD5Field) return -1;
  return g_d5_field > 0. ? 1 : -1;
}

//_____________________________________________________________________________
// read the D5 field and set the central momentum; on failure the previous
// central momentum is kept
void
UpdateD5Field()
{
  if(!gTM.IsReady()) return;
  if(g_d5_file.empty()){
    std::cout << "#W K18Tracking no D5FIELD: key in the conf file, keep"
              << " central momentum " << gTM.GetCentralMomentum()
              << " GeV/c" << std::endl;
    return;
  }
  std::ifstream ifs(g_d5_file);
  if(!ifs.is_open()){
    std::cout << "#W K18Tracking cannot open " << g_d5_file << ", keep"
              << " central momentum " << gTM.GetCentralMomentum()
              << " GeV/c" << std::endl;
    return;
  }
  std::string line;
  while(std::getline(ifs, line)){
    const auto comment = line.find('#');
    if(comment != std::string::npos) line.erase(comment);
    std::istringstream iss(line);
    Double_t field;
    if(!(iss >> field)) continue;
    g_d5_field = field;
    gTM.SetCentralMomentum(0.299792458*TMath::Abs(field)*kD5Radius*kFieldFactor);
    std::cout << "#D K18Tracking D5 = " << field << " T from " << g_d5_file
              << " -> central momentum " << gTM.GetCentralMomentum()
              << " GeV/c" << std::endl;
    return;
  }
  std::cout << "#W K18Tracking no D5 value in " << g_d5_file << ", keep"
            << " central momentum " << gTM.GetCentralMomentum()
            << " GeV/c" << std::endl;
}

struct FitHit
{
  Double_t z;     // hit z in the local tracking frame [mm]
  Double_t tilt;  // wire tilt angle [rad]
  Double_t s;     // measured local hit position [mm]
  Double_t reso;  // resolution [mm]
};
std::vector<FitHit> g_blc1_hits;
std::vector<FitHit> g_blc2_hits;

struct K18Fit
{
  Int_t    ibcin;
  Int_t    ibcout;
  Double_t chisqr;  // reduced
  Double_t par[kNParam];  // {xin, uin, yin, vin, delta}
  Double_t vo[4];         // matrix output at VO {xout, uout, yout, vout}
};

//_____________________________________________________________________________
Double_t
HitChisqr(const std::vector<FitHit>& hits, Double_t z0,
          Double_t x, Double_t u, Double_t y, Double_t v)
{
  Double_t chisqr = 0.;
  for(const auto& h : hits){
    const Double_t dz   = h.z - z0;
    const Double_t scal = (x+u*dz)*TMath::Cos(h.tilt)
                        + (y+v*dz)*TMath::Sin(h.tilt);
    const Double_t res  = (h.s - scal)/h.reso;
    chisqr += res*res;
  }
  return chisqr;
}

//_____________________________________________________________________________
// chi2 = BLC1 residuals (at VI) + BLC2 residuals (propagated to VO)
void
K18FitFCN(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f,
          Double_t* par, Int_t /*flag*/)
{
  Double_t parblc2[kNParam] = {0., 0., 0., 0., 0.};
  gTM.CalcBLC1toBLC2(par, parblc2, kMatrixOrder);
  f = HitChisqr(g_blc1_hits, kZVI, par[0], par[1], par[2], par[3])
    + HitChisqr(g_blc2_hits, kZVO,
                parblc2[0], parblc2[1], parblc2[2], parblc2[3]);
}

//_____________________________________________________________________________
void
FillFitHits(const DCLocalTrack* track, std::vector<FitHit>& out)
{
  out.clear();
  for(Int_t i=0, n=track->GetNHit(); i<n; ++i){
    const DCLTrackHit* h = track->GetHit(i);
    out.push_back({h->GetZ(), h->GetTiltAngle()*TMath::DegToRad(),
                   h->GetLocalHitPos(), h->GetResolution()});
  }
}

//_____________________________________________________________________________
// fit one BcIn-BcOut combination; true if MIGRAD converged
Bool_t
DoK18Fit(const DCLocalTrack* blc1, const DCLocalTrack* blc2, K18Fit& fit)
{
  FillFitHits(blc1, g_blc1_hits);
  FillFitHits(blc2, g_blc2_hits);

  const Int_t ndf = (Int_t)(g_blc1_hits.size()+g_blc2_hits.size()) - kNParam;
  if(ndf <= 0) return false;

  TMinuit minuit(kNParam);
  minuit.SetPrintLevel(-1);
  minuit.SetFCN(K18FitFCN);
  minuit.DefineParameter(0, "x",  blc1->GetX0() + blc1->GetU0()*kZVI, 0.1, 0., 0.);
  minuit.DefineParameter(1, "u",  blc1->GetU0(), 0.001, 0., 0.);
  minuit.DefineParameter(2, "y",  blc1->GetY0() + blc1->GetV0()*kZVI, 0.1, 0., 0.);
  minuit.DefineParameter(3, "v",  blc1->GetV0(), 0.001, 0., 0.);
  minuit.DefineParameter(4, "dp", 0., 0.001, 0., 0.);

  Double_t arglist[2] = {1000., 1.};
  Int_t ierflg = 0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
  if(ierflg != 0) return false;

  Double_t err;
  for(Int_t i=0; i<kNParam; ++i)
    minuit.GetParameter(i, fit.par[i], err);

  Double_t fmin, fedm, errdef;
  Int_t npari, nparx, istat;
  minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);
  fit.chisqr = fmin/ndf;

  Double_t p2[kNParam] = {0., 0., 0., 0., 0.};
  gTM.CalcBLC1toBLC2(fit.par, p2, kMatrixOrder);
  for(Int_t i=0; i<4; ++i) fit.vo[i] = p2[i];
  return true;
}

//_____________________________________________________________________________
// Beam RK through the HS field map, as offline HSTrack::Propagate
// (RK4 of NIM 160 (1979) 43, RK::PropagateOnce). Only the track state
// (x, y, u, v) is propagated: the transport-matrix derivatives of the
// offline code do not change it and are not needed here.
// Frame: HypTPC local [mm]; the start is the fitted VO state converted
// from the BcOut (FF) frame through DCGEO "HypTPC" (HSTrack::GlobalToHypTPC).
const Double_t kCHLB    = 2.99792458e-4; // [GeV/c / (T mm)]
const Double_t kRKStepZ = 10.;           // [mm], as HSTrack
const Int_t    kRKMaxStep = 10000;
// RK q = kShsMapChargeSign * charge / p (offline HSTrack, SHS map)
const Double_t kShsMapChargeSign = -1.;
// histogram planes in the HypTPC local frame [mm]
// Target plane: DCGEO "Target" converted to the HypTPC local frame in
// k18::Initialize() (DCGEO: global z = -293.1 -> local z = -143.0, the
// offline tpc::Z_TARGET); NaN without a Target entry
Double_t g_z_target = std::nan("");
// T2 plane: DCGEO "T2" converted to the HypTPC local frame in
// k18::Initialize() (online DCGEO: global z = 200 -> local z = 350.1);
// NaN without a T2 entry
Double_t g_z_t2 = std::nan("");
// BAC plane: DCGEO "BAC" in the HypTPC local frame (online DCGEO: global
// z = -648.4 -> local z = -498.3), NaN without a BAC entry
Double_t g_z_bac = std::nan("");

// T2 box and the HTOF segments behind it, set in k18::Initialize() from
// DCGEO "T2" (position) and DSIZE "T2" (full widths), HypTPC local [mm]
struct T2Box { Bool_t ok; Double_t x1, x2, y1, y2; };
T2Box g_t2_box = {false, 0., 0., 0., 0.};
// HTOF downstream face (azimuth 180 deg, z = +348.6 mm local, just in front
// of T2): segments 18-21 from -x to +x, 71 mm pitch (70 mm + clearance),
// as offline TPCAnalyzer::ExtrapolateToHTOF. Segment numbers 0-33 as the
// online HTOF histograms (0-based: HTOF_TDC_20_*, HitPat bin 20).
const Int_t    kHtofDownSeg0  = 18;
const Int_t    kHtofDownNSeg  = 4;
const Double_t kHtofSegWidth  = 71.;
std::vector<Int_t> g_t2_htof_segs; // segments overlapping the T2 x range

Bool_t
InT2Box(Double_t x, Double_t y)
{
  return g_t2_box.ok && g_t2_box.x1 < x && x < g_t2_box.x2
    && g_t2_box.y1 < y && y < g_t2_box.y2;
}

// BAC (DCGEO + DSIZE) and target boxes for the in-box fraction
// eff. = InBox / Tracks (no detector reference: BAC is a Cherenkov veto)
T2Box g_bac_box = {false, 0., 0., 0., 0.};
T2Box g_tgt_box = {false, 0., 0., 0., 0.};
// E72 cylindrical target seen from the beam: 2*radius x length [mm], as the
// outline on the Beam canvas (offline TPCPadHelper.hh)
const Double_t kTargetWX = 80.;
const Double_t kTargetWY = 100.;

Bool_t
InBox(const T2Box& b, Double_t x, Double_t y)
{
  return b.ok && b.x1 < x && x < b.x2 && b.y1 < y && y < b.y2;
}

// box of full widths (wx, wy) on the DCGEO entry; wx <= 0 takes DSIZE
T2Box
MakeBox(const TString& name, Double_t wx = -1., Double_t wy = -1.)
{
  try {
    const ThreeVector pos = gGeom.Global2LocalPos("HypTPC",
                                                  gGeom.GetGlobalPosition(name));
    if(wx <= 0.){
      const ThreeVector size = DetSizeMan::GetInstance().GetSize(name);
      wx = size.x();
      wy = size.y();
    }
    std::cout << "#D K18 " << name << " box x = [" << pos.x() - 0.5*wx << ", "
              << pos.x() + 0.5*wx << "], y = [" << pos.y() - 0.5*wy << ", "
              << pos.y() + 0.5*wy << "] mm" << std::endl;
    return {true, pos.x() - 0.5*wx, pos.x() + 0.5*wx,
            pos.y() - 0.5*wy, pos.y() + 0.5*wy};
  } catch(const std::exception& e) {
    std::cout << "#W K18 no " << name << " box (DCGEO / DSIZE), eff. disabled : "
              << e.what() << std::endl;
    return {false, 0., 0., 0., 0.};
  }
}

void
SetupT2Box()
{
  g_t2_box.ok = false;
  g_t2_htof_segs.clear();
  try {
    const ThreeVector pos = gGeom.Global2LocalPos("HypTPC",
                                                  gGeom.GetGlobalPosition("T2"));
    const ThreeVector size = DetSizeMan::GetInstance().GetSize("T2");
    g_t2_box = {true, pos.x() - 0.5*size.x(), pos.x() + 0.5*size.x(),
                pos.y() - 0.5*size.y(), pos.y() + 0.5*size.y()};
  } catch(const std::exception& e) {
    std::cout << "#W K18 no T2 box (DCGEO / DSIZE), T2 counts disabled : "
              << e.what() << std::endl;
    return;
  }
  /*
  for(Int_t i=0; i<kHtofDownNSeg; ++i){
    const Double_t x1 = (i - 0.5*kHtofDownNSeg) * kHtofSegWidth;
    const Double_t x2 = x1 + kHtofSegWidth;
    if(x1 < g_t2_box.x2 && g_t2_box.x1 < x2)
      g_t2_htof_segs.push_back(kHtofDownSeg0 + i);
  }
  */
  g_t2_htof_segs.push_back(20);
  std::cout << "#D K18 T2 box x = [" << g_t2_box.x1 << ", " << g_t2_box.x2
            << "], y = [" << g_t2_box.y1 << ", " << g_t2_box.y2
            << "] mm, HTOF segments :";
  for(auto seg : g_t2_htof_segs) std::cout << " " << seg;
  if(g_t2_htof_segs.empty()) std::cout << " none (T2 outside the HTOF face)";
  std::cout << std::endl;
}

// T2 efficiency text on the canvas (same style as the Multiplicity eff.):
// T2 detector eff. InBox&HTOF&T2 / InBox&HTOF from the count histogram
// (T2 hit = leading TDC in the T2_TDC window)
void
UpdateT2Text(TText* text, const TH1* h)
{
  if(!text || !h) return;
  const Double_t n_both = h->GetBinContent(4);
  const Double_t n_t2   = h->GetBinContent(5);
  text->SetTitle(Form("eff. %.3f", n_both > 0. ? n_t2/n_both : 0.));
}

// z of a DCGEO entry in the HypTPC local frame, NaN if missing
Double_t
LocalZFromGeom(const TString& name)
{
  try {
    const ThreeVector pos = gGeom.Global2LocalPos("HypTPC",
                                                  gGeom.GetGlobalPosition(name));
    std::cout << "#D K18 " << name << " plane from DCGEO : HypTPC local (x, y, z) = ("
              << pos.x() << ", " << pos.y() << ", " << pos.z() << ") mm"
              << std::endl;
    return pos.z();
  } catch(const std::exception& e) {
    std::cout << "#W K18 no " << name << " in DCGEO, " << name
              << " histograms will stay empty : " << e.what() << std::endl;
    return std::nan("");
  }
}

struct RKState { Double_t x, y, z, u, v; };

void
RKKick(Double_t u, Double_t v, Double_t q, const ThreeVector& B,
       Double_t& kx, Double_t& ky)
{
  const Double_t fac = std::sqrt(1.+u*u+v*v);
  kx = fac*(u*v*B.x() - (1.+u*u)*B.y() + v*B.z())*q*kCHLB;
  ky = fac*((1.+v*v)*B.x() - u*v*B.y() - u*B.z())*q*kCHLB;
}

RKState
RKStep(const RKState& s, Double_t q)
{
  const Double_t dr = kRKStepZ/std::sqrt(1.+s.u*s.u+s.v*s.v);
  const ThreeVector Z1(s.x, s.y, s.z);
  Double_t k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y;
  RKKick(s.u, s.v, q, gField.GetField(Z1), k1x, k1y);
  const ThreeVector B2 = gField.GetField(
    Z1 + ThreeVector(0.5*dr*s.u + 0.125*dr*dr*k1x,
                     0.5*dr*s.v + 0.125*dr*dr*k1y, 0.5*dr));
  RKKick(s.u + 0.5*dr*k1x, s.v + 0.5*dr*k1y, q, B2, k2x, k2y);
  RKKick(s.u + 0.5*dr*k2x, s.v + 0.5*dr*k2y, q, B2, k3x, k3y);
  const ThreeVector B4 = gField.GetField(
    Z1 + ThreeVector(dr*s.u + 0.5*dr*dr*k3x,
                     dr*s.v + 0.5*dr*dr*k3y, dr));
  RKKick(s.u + dr*k3x, s.v + dr*k3y, q, B4, k4x, k4y);
  RKState n;
  n.z = s.z + dr;
  n.x = s.x + dr*s.u + dr*dr/6.*(k1x + k2x + k3x);
  n.y = s.y + dr*s.v + dr*dr/6.*(k1y + k2y + k3y);
  n.u = s.u + dr/6.*(k1x + 2.*(k2x + k3x) + k4x);
  n.v = s.v + dr/6.*(k1y + 2.*(k2y + k3y) + k4y);
  return n;
}

// propagate to the planes z_planes (ascending); false if not all reached
Bool_t
RKPropagate(RKState s, Double_t q, const std::vector<Double_t>& z_planes,
            std::vector<RKState>& out)
{
  out.clear();
  std::size_t ivp = 0;
  for(Int_t istep=0; istep<kRKMaxStep && ivp<z_planes.size(); ++istep){
    const RKState n = RKStep(s, q);
    while(ivp<z_planes.size() &&
          (s.z - z_planes[ivp])*(n.z - z_planes[ivp]) <= 0.){
      const Double_t f = (z_planes[ivp] - s.z)/(n.z - s.z);
      out.push_back({s.x + f*(n.x - s.x), s.y + f*(n.y - s.y), z_planes[ivp],
                     s.u + f*(n.u - s.u), s.v + f*(n.v - s.v)});
      ++ivp;
    }
    s = n;
  }
  return ivp == z_planes.size();
}

}

namespace k18
{
//_____________________________________________________________________________
void
Initialize(const std::string& d5_file)
{
  g_d5_file = d5_file;
  kZVI = 0.5*(gGeom.GetLocalZ("BLC1a-U1") + gGeom.GetLocalZ("BLC1b-VP2"));
  kZVO = 0.5*(gGeom.GetLocalZ("BLC2a-U1") + gGeom.GetLocalZ("BLC2b-VP2"));
  std::cout << "#D K18 matrix planes from DCGEO : VI z = " << kZVI
            << " mm, VO z = " << kZVO << " mm" << std::endl;
  g_z_target = LocalZFromGeom("Target");
  g_z_t2     = LocalZFromGeom("T2");
  g_z_bac    = LocalZFromGeom("BAC");
  SetupT2Box();
  g_bac_box = MakeBox("BAC");
  g_tgt_box = MakeBox("Target", kTargetWX, kTargetWY);
  UpdateD5Field();
}

//_____________________________________________________________________________
void
ReloadD5Field()
{
  UpdateD5Field();
}

//_____________________________________________________________________________
std::vector<BeamTrack>
FitBeamTracks(const DCAnalyzer& dcAna)
{
  // re-read the D5 field at every run change
  static const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
  static UInt_t s_run = gUnpacker.get_run_number();
  if(s_run != gUnpacker.get_run_number()){
    s_run = gUnpacker.get_run_number();
    UpdateD5Field();
  }

  std::vector<BeamTrack> tracks;
  const auto& contIn  = dcAna.GetBcInTrackContainer();
  const auto& contOut = dcAna.GetBcOutTrackContainer();
  if(!gTM.IsReady() || contIn.empty() || contOut.empty()) return tracks;

  // combinatorial fit of all BcIn-BcOut pairs
  const Int_t nIn  = TMath::Min((Int_t)contIn.size(),  kMaxNTrack);
  const Int_t nOut = TMath::Min((Int_t)contOut.size(), kMaxNTrack);
  std::vector<K18Fit> cands;
  for(Int_t i=0; i<nIn; ++i){
    for(Int_t j=0; j<nOut; ++j){
      K18Fit fit;
      fit.ibcin  = i;
      fit.ibcout = j;
      if(DoK18Fit(contIn[i], contOut[j], fit))
        cands.push_back(fit);
    }
  }
  std::sort(cands.begin(), cands.end(),
            [](const K18Fit& a, const K18Fit& b){ return a.chisqr < b.chisqr; });

  // greedy exclusive pairing: a local track joins only one K18 track
  std::vector<Bool_t> used1(nIn, false), used2(nOut, false);
  const Double_t p0 = gTM.GetCentralMomentum();
  for(const auto& fit : cands){
    if(used1[fit.ibcin] || used2[fit.ibcout]) continue;
    used1[fit.ibcin]  = true;
    used2[fit.ibcout] = true;
    BeamTrack t;
    t.delta  = fit.par[4];
    t.p      = p0*(1. + fit.par[4]);
    t.chisqr = fit.chisqr;
    for(Int_t i=0; i<4; ++i) t.vo[i] = fit.vo[i];
    tracks.push_back(t);
  }
  return tracks;
}

//_____________________________________________________________________________
// RK start state of a beam track in the HypTPC local frame
Bool_t
StartState(const BeamTrack& t, RKState& start, Double_t& q)
{
  if(!gField.IsReady()) return false;
  if(t.p <= 0. || !std::isfinite(t.p)) return false;
  const ThreeVector pos = gGeom.Global2LocalPos("HypTPC",
    ThreeVector(t.vo[0], t.vo[2], kZVO));
  const ThreeVector dir = gGeom.Global2LocalDir("HypTPC",
    ThreeVector(t.vo[1], t.vo[3], 1.));
  if(TMath::Abs(dir.z()) < 1.e-12) return false;
  start = {pos.x(), pos.y(), pos.z(), dir.x()/dir.z(), dir.y()/dir.z()};
  q = kShsMapChargeSign*BeamCharge()/t.p;
  return true;
}

//_____________________________________________________________________________
Bool_t
TraceBeam(const BeamTrack& t, Double_t z_end, std::vector<ThreeVector>& points)
{
  points.clear();
  RKState s;
  Double_t q;
  if(!StartState(t, s, q)) return false;
  points.emplace_back(s.x, s.y, s.z);
  for(Int_t istep=0; istep<kRKMaxStep; ++istep){
    const RKState n = RKStep(s, q);
    if(!std::isfinite(n.z) || n.z <= s.z) return false;
    if(n.z >= z_end){
      const Double_t f = (z_end - s.z)/(n.z - s.z);
      points.emplace_back(s.x + f*(n.x - s.x), s.y + f*(n.y - s.y), z_end);
      return true;
    }
    points.emplace_back(n.x, n.y, n.z);
    s = n;
  }
  return false;
}

//_____________________________________________________________________________
Double_t
LocalZ(const TString& name)
{
  return LocalZFromGeom(name);
}

//_____________________________________________________________________________
void
Fill(const DCAnalyzer& dcAna, std::vector<TH1*>& hptr_array, Int_t particle,
     const std::vector<Int_t>& htof_segs, Bool_t t2_hit)
{
  static const Int_t multi_id  = HistMaker::getSequentialID(kK18Tracking, 0, kMulti,   1);
  static const Int_t mom_id    = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy1, 1);
  static const Int_t mom_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy1, 2);
  static const Int_t mom_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy1, 3);
  static const Int_t delta_id  = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy2, 1);
  static const Int_t chisqr_id = HistMaker::getSequentialID(kK18Tracking, 0, kChisqr,  1);
  static const Int_t rk_status_id = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 2);
  static const Int_t tgt_xy_id  = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 1);
  static const Int_t t2_xy_id   = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 2);
  static const Int_t tgt_x_id   = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat, 1);
  static const Int_t tgt_y_id   = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat, 2);
  static const Int_t t2_x_id    = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat, 3);
  static const Int_t t2_y_id    = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat, 4);
  // 2D per particle (beam TOF)
  static const Int_t tgt_xy_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 3);
  static const Int_t tgt_xy_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 4);
  static const Int_t t2_xy_pi_id  = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 5);
  static const Int_t t2_xy_k_id   = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 6);
  // BAC XY: all, Pi, K (no 1D)
  static const Int_t bac_xy_id    = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 7);
  static const Int_t bac_xy_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 8);
  static const Int_t bac_xy_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kHitPat2D, 9);
  // in-box counts (Tracks, InBox): BAC Pi / K, Target Pi / K
  static const Int_t bac_cnt_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 5);
  static const Int_t bac_cnt_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 6);
  static const Int_t tgt_cnt_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 7);
  static const Int_t tgt_cnt_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 8);
  // T2 counts per particle: Tracks, InBox, HTOF, InBox&HTOF, InBox&HTOF&T2
  static const Int_t t2_cnt_pi_id = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 3);
  static const Int_t t2_cnt_k_id  = HistMaker::getSequentialID(kK18Tracking, 0, kMulti, 4);

  // planes in ascending z, with the histograms to fill at each
  // box: in-box fraction for Pi / K (count ids, -1 without)
  struct Plane { Double_t z; Int_t xy_id, x_id, y_id, xy_pi_id, xy_k_id;
                 const T2Box* box = nullptr; Int_t cnt_pi_id = -1, cnt_k_id = -1;
                 TString name = ""; };
  static const std::vector<Plane> planes = []{
    std::vector<Plane> v;
    if(std::isfinite(g_z_target))
      v.push_back({g_z_target, tgt_xy_id, tgt_x_id, tgt_y_id, tgt_xy_pi_id, tgt_xy_k_id,
                   &g_tgt_box, tgt_cnt_pi_id, tgt_cnt_k_id, "Target"});
    if(std::isfinite(g_z_t2))
      v.push_back({g_z_t2, t2_xy_id, t2_x_id, t2_y_id, t2_xy_pi_id, t2_xy_k_id});
    if(std::isfinite(g_z_bac))
      v.push_back({g_z_bac, bac_xy_id, -1, -1, bac_xy_pi_id, bac_xy_k_id,
                   &g_bac_box, bac_cnt_pi_id, bac_cnt_k_id, "BAC"});
    std::sort(v.begin(), v.end(),
              [](const Plane& a, const Plane& b){ return a.z < b.z; });
    return v;
  }();
  static const std::vector<Double_t> plane_z = []{
    std::vector<Double_t> z;
    for(const auto& pl : planes) z.push_back(pl.z);
    return z;
  }();

  // beam really at T2: a U&D hit (HTOF_TDC window) in a segment behind T2
  Bool_t htof_hit = false;
  for(auto seg : htof_segs)
    if(std::find(g_t2_htof_segs.begin(), g_t2_htof_segs.end(), seg)
       != g_t2_htof_segs.end()) htof_hit = true;
  const Int_t t2_cnt_id = particle == kPi ? t2_cnt_pi_id
    : particle == kKaon ? t2_cnt_k_id : -1;

  const auto tracks = FitBeamTracks(dcAna);
  if(!gTM.IsReady() || dcAna.GetBcInTrackContainer().empty()
     || dcAna.GetBcOutTrackContainer().empty()) return;
  for(std::size_t it=0; it<tracks.size(); ++it){
    const auto& t = tracks[it];
    hptr_array[mom_id]->Fill(t.p);
    if(particle == kPi)   hptr_array[mom_pi_id]->Fill(t.p);
    if(particle == kKaon) hptr_array[mom_k_id]->Fill(t.p);
    hptr_array[delta_id]->Fill(t.delta);
    hptr_array[chisqr_id]->Fill(t.chisqr);

    // beam RK from VO to the target and T2
    RKState start;
    Double_t q;
    if(planes.empty() || !StartState(t, start, q)) continue;
    std::vector<RKState> vp;
    const Bool_t passed = RKPropagate(start, q, plane_z, vp);
    hptr_array[rk_status_id]->Fill(passed ? 1 : 0);
    if(!passed) continue;
    for(std::size_t i=0; i<planes.size(); ++i){
      hptr_array[planes[i].xy_id]->Fill(vp[i].x, vp[i].y);
      if(planes[i].x_id >= 0) hptr_array[planes[i].x_id]->Fill(vp[i].x);
      if(planes[i].y_id >= 0) hptr_array[planes[i].y_id]->Fill(vp[i].y);
      if(particle == kPi)   hptr_array[planes[i].xy_pi_id]->Fill(vp[i].x, vp[i].y);
      if(particle == kKaon) hptr_array[planes[i].xy_k_id]->Fill(vp[i].x, vp[i].y);
      // BAC / Target in-box fraction with the best track of the event
      const Int_t cnt_id = particle == kPi ? planes[i].cnt_pi_id
        : particle == kKaon ? planes[i].cnt_k_id : -1;
      if(it == 0 && cnt_id >= 0 && planes[i].box && planes[i].box->ok){
        TH1* h = hptr_array[cnt_id];
        h->Fill(0);
        if(InBox(*planes[i].box, vp[i].x, vp[i].y)) h->Fill(1);
        if(auto text = BoxEffText(planes[i].name, particle)){
          const Double_t n_trk = h->GetBinContent(1);
          text->SetTitle(Form("eff. %.3f",
                              n_trk > 0. ? h->GetBinContent(2)/n_trk : 0.));
        }
      }
      // T2 counts with the best track of the event (one beam per event)
      if(it == 0 && t2_cnt_id >= 0 && planes[i].z == g_z_t2 && g_t2_box.ok){
        const Bool_t in_box = InT2Box(vp[i].x, vp[i].y);
        TH1* h = hptr_array[t2_cnt_id];
        h->Fill(0);
        if(in_box)             h->Fill(1);
        if(htof_hit)           h->Fill(2);
        if(in_box && htof_hit) h->Fill(3);
        if(in_box && htof_hit && t2_hit) h->Fill(4);
        UpdateT2Text(T2CountText(particle), h);
      }
    }
  }
  hptr_array[multi_id]->Fill((Int_t)tracks.size());
}

//_____________________________________________________________________________
Double_t
T2CenterX()
{
  return g_t2_box.ok ? 0.5*(g_t2_box.x1 + g_t2_box.x2) : 0.;
}

//_____________________________________________________________________________
TText*
T2CountText(Int_t particle)
{
  auto make = [](const char* name){
    auto t = new TText(0.150, 0.800, "eff. 0.000");
    t->SetName(name);
    t->SetNDC();
    t->SetTextSize(0.080);
    return t;
  };
  static TText* text_pi = make("Beam_T2_Eff_Pi_Text");
  static TText* text_k  = make("Beam_T2_Eff_K_Text");
  if(particle == kPi)   return text_pi;
  if(particle == kKaon) return text_k;
  return nullptr;
}

//_____________________________________________________________________________
TText*
BoxEffText(const TString& name, Int_t particle)
{
  auto make = [](const char* tname){
    auto t = new TText(0.150, 0.800, "eff. 0.000");
    t->SetName(tname);
    t->SetNDC();
    t->SetTextSize(0.080);
    return t;
  };
  static TText* bac_pi = make("Beam_BAC_Eff_Pi_Text");
  static TText* bac_k  = make("Beam_BAC_Eff_K_Text");
  static TText* tgt_pi = make("Beam_Target_Eff_Pi_Text");
  static TText* tgt_k  = make("Beam_Target_Eff_K_Text");
  if(name == "BAC")
    return particle == kPi ? bac_pi : particle == kKaon ? bac_k : nullptr;
  if(name == "Target")
    return particle == kPi ? tgt_pi : particle == kKaon ? tgt_k : nullptr;
  return nullptr;
}

//_____________________________________________________________________________
TText*
T2SegText(Int_t particle)
{
  auto make = [](const char* name){
    TString segs;
    for(auto seg : g_t2_htof_segs)
      segs += Form("%s%d", segs.IsNull() ? "" : ",", seg);
    if(segs.IsNull()) segs = "none";
    auto t = new TText(0.150, 0.720, Form("HTOF seg %s", segs.Data()));
    t->SetName(name);
    t->SetNDC();
    t->SetTextSize(0.060);
    return t;
  };
  static TText* text_pi = make("Beam_T2_Seg_Pi_Text");
  static TText* text_k  = make("Beam_T2_Seg_K_Text");
  if(particle == kPi)   return text_pi;
  if(particle == kKaon) return text_k;
  return nullptr;
}

//_____________________________________________________________________________
// one point outside any frame, so that an unused line draws nothing
void
HideGraph(TGraph* gr)
{
  gr->Set(1);
  gr->SetPoint(0, -1.e6, -1.e6);
}

//_____________________________________________________________________________
std::vector<TGraph*>&
EventDisplayGraphs()
{
  static std::vector<TGraph*> graphs = []{
    std::vector<TGraph*> g;
    for(Int_t i=0; i<kMaxDisplayTracks; ++i){
      auto gr = new TGraph();
      gr->SetName(Form("BeamRK_%d", i));
      gr->SetLineColor(kRed);
      gr->SetLineWidth(2);
      HideGraph(gr);
      g.push_back(gr);
    }
    return g;
  }();
  return graphs;
}

//_____________________________________________________________________________
void
UpdateEventDisplay(const DCAnalyzer& dcAna, TH1* h_target_x)
{
  // RK end plane: SCH from DCGEO (HypTPC local), as the event display frame
  static const Double_t z_sch = LocalZFromGeom("SCH");
  auto& graphs = EventDisplayGraphs();
  // unused lines keep one point far outside the frame: a 0-point TGraph
  // gives "TGraphPainter::PaintGraph: illegal number of points (0)"
  for(auto& gr : graphs) HideGraph(gr);
  if(!std::isfinite(z_sch)) return;
  const auto tracks = FitBeamTracks(dcAna);
  std::vector<ThreeVector> points;
  for(std::size_t i=0; i<tracks.size() && i<graphs.size(); ++i){
    if(!TraceBeam(tracks[i], z_sch, points)) continue;
    // event display: horizontal = z, vertical = x [mm]
    graphs[i]->Set(points.size());
    for(std::size_t j=0; j<points.size(); ++j)
      graphs[i]->SetPoint(j, points[j].z(), points[j].x());
    // x at the target plane, interpolated between the RK points
    if(!h_target_x || !std::isfinite(g_z_target)) continue;
    for(std::size_t j=1; j<points.size(); ++j){
      const auto& a = points[j-1];
      const auto& b = points[j];
      if(a.z() <= g_z_target && g_z_target <= b.z() && b.z() > a.z()){
        const Double_t f = (g_z_target - a.z())/(b.z() - a.z());
        h_target_x->Fill(a.x() + f*(b.x() - a.x()));
        break;
      }
    }
  }
}
}

//______________________________________________________________________________
void
ConfMan::initializeK18Tracking()
{
  k18::Initialize(name_file_["D5FIELD:"]);
}
