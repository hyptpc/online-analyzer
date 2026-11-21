// -*- C++ -*-

#include "K18Tracking.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TH1.h>
#include <TMath.h>
#include <TMinuit.h>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCGeomMan.hh"
#include "DCLocalTrack.hh"
#include "DCLTrackHit.hh"
#include "HistMaker.hh"
#include "TransferMatrixMan.hh"

#include <UnpackerManager.hh>

namespace
{
auto& gTM = TransferMatrixMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();

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
  return true;
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
  UpdateD5Field();
}

//_____________________________________________________________________________
void
ReloadD5Field()
{
  UpdateD5Field();
}

//_____________________________________________________________________________
void
Fill(const DCAnalyzer& dcAna, std::vector<TH1*>& hptr_array)
{
  // re-read the D5 field at every run change
  static const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
  static Int_t s_run = gUnpacker.get_run_number();
  if(s_run != gUnpacker.get_run_number()){
    s_run = gUnpacker.get_run_number();
    UpdateD5Field();
  }

  const auto& contIn  = dcAna.GetBcInTrackContainer();
  const auto& contOut = dcAna.GetBcOutTrackContainer();
  if(!gTM.IsReady() || contIn.empty() || contOut.empty()) return;

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

  static const Int_t multi_id  = HistMaker::getSequentialID(kK18Tracking, 0, kMulti,   1);
  static const Int_t mom_id    = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy1, 1);
  static const Int_t delta_id  = HistMaker::getSequentialID(kK18Tracking, 0, kEnergy2, 1);
  static const Int_t chisqr_id = HistMaker::getSequentialID(kK18Tracking, 0, kChisqr,  1);

  // greedy exclusive pairing: a local track joins only one K18 track
  std::vector<Bool_t> used1(nIn, false), used2(nOut, false);
  const Double_t p0 = gTM.GetCentralMomentum();
  Int_t ntK18 = 0;
  for(const auto& fit : cands){
    if(used1[fit.ibcin] || used2[fit.ibcout]) continue;
    used1[fit.ibcin]  = true;
    used2[fit.ibcout] = true;
    ++ntK18;
    hptr_array[mom_id]->Fill(p0*(1. + fit.par[4]));
    hptr_array[delta_id]->Fill(fit.par[4]);
    hptr_array[chisqr_id]->Fill(fit.chisqr);
  }
  hptr_array[multi_id]->Fill(ntK18);
}
}

//______________________________________________________________________________
void
ConfMan::initializeK18Tracking()
{
  k18::Initialize(name_file_["D5FIELD:"]);
}
