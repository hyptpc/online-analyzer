// -*- C++ -*-

#include "EventAnalyzer.hh"

#include <iostream>
#include <string>
#include <vector>

#include <TMath.h>
#include <TString.h>

#include "DCAnalyzer.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DetectorID.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#include <DAQNode.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
const auto& gUConf = hddaq::unpacker::GConfig::get_instance();
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
EventAnalyzer::EventAnalyzer()
{
}

//_____________________________________________________________________________
EventAnalyzer::~EventAnalyzer()
{
}


//_____________________________________________________________________________
void
EventAnalyzer::DCRawHit(const TString& dcname, const RawData& rawData)
{
  static const auto& digit_info = gUConf.get_digit_info();
  for (const auto& name_str : DCNameList.at(dcname)) {
    const auto name = name_str.Data();
    Int_t detector_id = digit_info.get_device_id(name);
    Int_t nplane = digit_info.get_n_plane(detector_id);
    for(Int_t plane=0; plane<nplane; ++plane){
      const auto& cont = rawData.GetDCRawHC(detector_id, plane);
      Int_t multi = 0;
      Int_t cmulti = 0;
      for(Int_t i=0, n=cont.size(); i<n; ++i){
        auto hit = cont[i];
        if(!hit) continue;
        auto wire = hit->WireId();
        Bool_t is_good = false;
        for(Int_t j=0, m=hit->GetTdcSize(); j<m; ++j){
          if(m != hit->GetTrailingSize()) break;
          auto l = hit->GetTdc(j);
          auto t = hit->GetTrailing(j);
          auto tot = l - t;
          if(gUser.IsInRange(Form("%s_TOT", name), tot) &&
             gUser.IsInRange(Form("%s_TDC", name), l))
            is_good = true;
          for(const auto& totcut: std::vector<TString>{"", "C"}){
            if(totcut == "C" && !gUser.IsInRange(Form("%s_TOT", name), tot))
              continue;
            auto c = totcut.Data();
          }
        }
        if(hit->GetTdcSize() > 0){
          ++multi;
        }
        if(is_good){
          ++cmulti;
        }
      }
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::DCHit(const TString& dcname, const DCAnalyzer& dcAna)
{
  static const auto& digit_info = gUConf.get_digit_info();
  // DC
  for (const auto& name_str : DCNameList.at(dcname)) {
    const auto name = name_str.Data();
    auto detector_id = digit_info.get_device_id(name);
    Int_t nplane = digit_info.get_n_plane(detector_id);
    for(Int_t plane=0; plane<nplane; ++plane){
      Int_t plane_idx = plane;
      if (detector_id == DetIdBLC1b || detector_id == DetIdBLC2b) plane_idx = plane + nplane; 
      auto hc1 = dcname == "BcIn" ?
                 dcAna.GetBcInHC(plane_idx) : dcAna.GetBcOutHC(plane_idx);
      Int_t multi = 0;
      if(plane%2 == 0) { //pair plane hit pattern start
        Int_t wire1, wire2;
        for(const auto &hit1: hc1){
          wire1 = hit1->GetWire();
          auto hc2 = dcname == "BcIn" ?
                     dcAna.GetBcInHC(plane_idx) : dcAna.GetBcOutHC(plane_idx);
          for(const auto &hit2: hc2){
            wire2 = hit2->GetWire();
          }
        }
      } //pair plane hit pattern end

      for(const auto& hit: hc1){
        auto wire = hit->GetWire();	
        Bool_t is_good = false;
        for(Int_t j=0, m=hit->GetDriftTimeSize(); j<m; ++j){
          if(!hit->IsGood(j)) continue;
          auto dt  = hit->GetDriftTime(j);
          auto dl  = hit->GetDriftLength(j);
          auto tot = hit->GetTot(j);
          is_good = true;
        }
        if(is_good){
          ++multi;
        }
      }
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::BcInTracking(DCAnalyzer& dcAna)
{
  // static const auto& digit_info = gUConf.get_digit_info();

  for (const auto& track : dcAna.GetBcInTrackContainer()) {
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t x0 = track->GetX0();
    Double_t y0 = track->GetY0();
    Double_t u0 = track->GetU0();
    Double_t v0 = track->GetV0();
    // Double_t theta = track->GetTheta();

    for (const auto& lthit : track->GetHitArray()) {
      const auto hit  = lthit->GetHit();
      const auto name = hit->GetRawHit()->DetectorName().Data();
      auto plane      = hit->PlaneId();
      auto wire       = lthit->GetWire();
      auto dt         = lthit->DriftTime();
      auto dl         = lthit->DriftLength();
      auto res  = lthit->GetResidual();
      auto wp   = lthit->GetWirePosition();
      auto pos  = lthit->GetLocalHitPos();
      auto sign = (pos - wp > 0.) ? 1 : -1;
    }
  }
}

//_____________________________________________________________________________
void
EventAnalyzer::BcOutTracking(DCAnalyzer& dcAna)
{
  // static const auto& digit_info = gUConf.get_digit_info();
  for (const auto& track : dcAna.GetBcOutTrackContainer()) {
    Int_t nh = track->GetNHit();
    Double_t chisqr = track->GetChiSquare();
    Double_t x0 = track->GetX0();
    Double_t y0 = track->GetY0();
    Double_t u0 = track->GetU0();
    Double_t v0 = track->GetV0();
    // Double_t theta = track->GetTheta();

    for (const auto& lthit : track->GetHitArray()) {
      const auto hit  = lthit->GetHit();
      const auto name = hit->GetRawHit()->DetectorName().Data();
      auto plane      = hit->PlaneId();
      auto wire       = lthit->GetWire();
      auto dt         = lthit->DriftTime();
      auto dl         = lthit->DriftLength();
      auto res  = lthit->GetResidual();
      auto wp   = lthit->GetWirePosition();
      auto pos  = lthit->GetLocalHitPos();
      auto sign = (pos - wp > 0.) ? 1 : -1;
    }
  }
}

