// -*- C++ -*-

#include <iostream>
#include <string>
#include <vector>

#include <TGFileBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TString.h>

#include <Unpacker.hh>
#include <UnpackerManager.hh>

#include "user_analyzer.hh"

#include "BH2Filter.hh"
#include "BH2Hit.hh"
#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCTdcCalibMan.hh"
#include "DetectorID.hh"
#include "EventAnalyzer.hh"
#include "EventDisplay.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FieldMan.hh"
#include "HistMaker.hh"
#include "HodoAnalyzer.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "MacroBuilder.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#define BH2FILTER 1

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace
  {
    const UnpackerManager& gUnpacker  = GUnpacker::get_instance();
    const HistMaker&       gHist      = HistMaker::getInstance();
    const DCGeomMan&       gGeom      = DCGeomMan::GetInstance();
          BH2Filter&       gBH2Filter = BH2Filter::GetInstance();
          EventDisplay&    gEvDisp    = EventDisplay::GetInstance();
    const UserParamMan&    gUser      = UserParamMan::GetInstance();

    std::vector<TH1*> hptr_array;

  }

//____________________________________________________________________________
Int_t
process_begin( const std::vector<std::string>& argv )
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  gConfMan.InitializeParameter<BH2Filter>("BH2FLT");
  gConfMan.InitializeParameter<DCGeomMan>("DCGEOM");
  gConfMan.InitializeParameter<DCTdcCalibMan>("TDCCALIB");
  gConfMan.InitializeParameter<DCDriftParamMan>("DRFTPM");
  gConfMan.InitializeParameter<HodoParamMan>("HDPRM");
  gConfMan.InitializeParameter<HodoPHCMan>("HDPHC");
  gConfMan.InitializeParameter<FieldMan>("KURAMA");
  gConfMan.InitializeParameter<UserParamMan>("USER");
  gConfMan.InitializeParameter<EventDisplay>();
  if( !gConfMan.IsGood() ) return -1;

  gBH2Filter.Print();
  // gBH2Filter.SetVerbose();

  gEvDisp.SetStyle();

  return 0;
}

//____________________________________________________________________________
Int_t
process_end( void )
{
  return 0;
}

//____________________________________________________________________________
Int_t
process_event( void )
{
  //  static const Double_t MaxMultiHitBcOut  = gUser.GetParameter("MaxMultiHitBcOut");
  static const Double_t MaxMultiHitSdcIn  = gUser.GetParameter("MaxMultiHitSdcIn");
  static const Double_t MaxMultiHitSdcOut = gUser.GetParameter("MaxMultiHitSdcOut");

  EventAnalyzer event;
  event.DecodeRawData();
  event.DecodeDCAnalyzer();
  event.DecodeHodoAnalyzer();
  // event.TimeCutBFT();

  const RawData*      const rawData = event.GetRawData();
  const HodoAnalyzer* const hodoAna = event.GetHodoAnalyzer();
  const DCAnalyzer*   const dcAna   = event.GetDCAnalyzer();

  gEvDisp.DrawText();

  // Trigger Flag
  {
    const HodoRHitContainer &cont = rawData->GetTrigRawHC();
    std::vector<Int_t> trigflag(NumOfSegTFlag);
    Int_t nh = cont.size();
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId()+1;
      Int_t tdc = hit->GetTdc1();
      trigflag[seg-1] = tdc;
    }
    gEvDisp.DrawTrigger( trigflag );
    if (trigflag[trigger::kSpillOnEnd] || trigflag[trigger::kSpillOffEnd])
      return 0;
  }

  // BH1
  {
    const HodoRHitContainer &cont = rawData->GetBH1RawHC();
    Int_t nh = cont.size();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      if( !hit ) continue;
      Int_t Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if( Tu>0 || Td>0 ) hit_flag = true;
    }
    gEvDisp.PrintHit("BH1", 0.86, hit_flag);
  }

  // BFT
  std::vector<Double_t> BftXCont;
  {
    Int_t ncl = hodoAna->GetNClustersBFT();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<ncl; ++i ){
      FiberCluster *cl = hodoAna->GetClusterBFT(i);
      if( !cl ) continue;
      Double_t pos = cl->MeanPosition();
      BftXCont.push_back( pos );
      hit_flag = true;
    }
    gEvDisp.PrintHit("BFT", 0.83, hit_flag);
  }

  // BH2
  {
    const HodoRHitContainer &cont = rawData->GetBH2RawHC();
    Int_t nh=cont.size();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      Int_t seg = hit->SegmentId();
      if( !hit ) continue;
      Int_t Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if( Tu>0 || Td>0 ){
	hit_flag = true;
	gEvDisp.DrawHitHodoscope( gGeom.GetDetectorId("BH2"), seg, Tu, Td );
      }
    }
    gEvDisp.PrintHit("BH2", 0.80, hit_flag);
  }

  // SFT
  {
    for( Int_t l=0; l<NumOfLayersSFT; ++l ){
      Int_t ncl = hodoAna->GetNClustersSFT(l);
      Bool_t hit_flag = false;
      for( Int_t i=0; i<ncl; ++i ){
	FiberCluster *cl = hodoAna->GetClusterSFT(l,i);
	if( !cl ) continue;
	Int_t seg   = (Int_t)cl->MeanSeg();
	Int_t plane = cl->PlaneId();
	hit_flag = true;
	gEvDisp.DrawHitWire( plane+7, seg );
      }
      gEvDisp.PrintHit("SFT", 0.77, hit_flag);
    }
  }

  // SAC
  {
    const HodoRHitContainer &cont = rawData->GetSACRawHC();
    Int_t nh=cont.size();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      if( !hit ) continue;
      Int_t Tu=hit->GetTdcUp(), Td=hit->GetTdcDown();
      if( Tu>0 || Td>0 ) hit_flag = true;
    }
    gEvDisp.PrintHit("SAC", 0.74, hit_flag);
  }

  // SCH
  {
    Int_t nhSch = hodoAna->GetNHitsSCH();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<nhSch; ++i ){
      FiberHit *hit = hodoAna->GetHitSCH(i);
      if( !hit ) continue;
      Int_t mh  = hit->GetNumOfHit();
      Int_t seg = hit->SegmentId();
      bool flag = false;
      for( Int_t m=0; m<mh; ++m ){
        // Double_t leading = hit->GetLeading(m);
        // if( MinTdcSCH<leading && leading<MaxTdcSCH ){
          flag = true;
        // }
      }
      if( flag ){
        hit_flag = true;
        gEvDisp.DrawHitHodoscope( gGeom.GetDetectorId("SCH"), seg );
      }
    }
    gEvDisp.PrintHit("SCH", 0.71, hit_flag);
  }

  // TOF
  {
    const HodoRHitContainer &cont = rawData->GetTOFRawHC();
    Int_t nh = cont.size();
    Bool_t hit_flag = false;
    for( Int_t i=0; i<nh; ++i ){
      HodoRawHit *hit = cont[i];
      if( !hit ) continue;
      Int_t seg = hit->SegmentId();
      Int_t Tu = hit->GetTdcUp(), Td = hit->GetTdcDown();
      if( Tu>0 || Td>0 ){
        hit_flag = true;
        gEvDisp.DrawHitHodoscope( gGeom.GetDetectorId("TOF"), seg, Tu, Td );
      }
    }
    gEvDisp.PrintHit("TOF", 0.68, hit_flag);
  }
  // if( nhTof==0 ) return true;

  // BcOut
  // Double_t multi_BcOut = 0.;
  // for( Int_t layer=1; layer<=NumOfLayersBcOut; ++layer ){
  //   const DCHitContainer &cont = dcAna->GetBcOutHC(layer);
  //   Int_t nhIn=cont.size();
  //   for( Int_t i=0; i<nhIn; ++i ){
  //     DCHit*   hit  = cont[i];
  //     Double_t wire = hit->GetWire();
  //     Int_t      mhit = hit->GetTdcSize();
  //     bool     goodFlag = false;
  //     ++multi_BcOut;
  //     for (Int_t j=0; j<mhit; j++) {
  // 	if (hit->IsWithinRange(j)) {
  // 	  goodFlag = true;
  // 	  break;
  // 	}
  //     }
  //     if( goodFlag ){
  // 	hddaq::cout << "BcOut : " << wire << std::endl;
  // 	gEvDisp.DrawHitWire( layer+112, Int_t(wire) );
  //     }
  //   }
  // }
  // multi_BcOut /= (Double_t)NumOfLayersBcOut;
  // if( multi_BcOut > MaxMultiHitBcOut ) return true;

  // SdcIn
  Double_t multi_SdcIn = 0.;
  for( int layer=1; layer<=NumOfLayersSdcIn; ++layer ){
    const DCHitContainer &contIn = dcAna->GetSdcInHC(layer);
    int nhIn=contIn.size();
    if ( nhIn > MaxMultiHitSdcIn )
      continue;
    for( int i=0; i<nhIn; ++i ){
      DCHit  *hit  = contIn[i];
      Double_t  wire = hit->GetWire();
      int     mhit = hit->GetTdcSize();
      bool    goodFlag = false;
      ++multi_SdcIn;
      for (int j=0; j<mhit; j++) {
	if (hit->IsWithinRange(j)) {
	  goodFlag = true;
	  break;
	}
      }
      if( goodFlag )
	gEvDisp.DrawHitWire( layer, int(wire) );
    }
  }
  // multi_SdcIn /= (Double_t)NumOfLayersSdcIn;
  // if( multi_SdcIn > MaxMultiHitSdcIn ) return true;

  // SdcOut
  Double_t multi_SdcOut = 0.;
  for( int layer=1; layer<=NumOfLayersSdcOut; ++layer ){
    const DCHitContainer &contOut = dcAna->GetSdcOutHC(layer);
    int nhOut = contOut.size();
    if ( nhOut > MaxMultiHitSdcOut )
      continue;
    for( int i=0; i<nhOut; ++i ){
      DCHit  *hit  = contOut[i];
      Double_t  wire = hit->GetWire();
      ++multi_SdcOut;
      gEvDisp.DrawHitWire( layer+30, int(wire) );
    }
  }
  // multi_SdcOut /= (Double_t)NumOfLayersSdcOut;
  // if( multi_SdcOut > MaxMultiHitSdcOut ) return true;

  event.ApplyBH2Filter();
  event.TrackSearchBcOut();

  Int_t ntBcOut = dcAna->GetNtracksBcOut();
  for( Int_t it=0; it<ntBcOut; ++it ){
    DCLocalTrack *tp = dcAna->GetTrackBcOut( it );
    if( tp ) gEvDisp.DrawBcOutLocalTrack( tp );
  }

  // if( ntBcOut==0 ) return true;

  event.TrackSearchSdcIn();
  Int_t ntSdcIn = dcAna->GetNtracksSdcIn();
  for( Int_t it=0; it<ntSdcIn; ++it ){
    DCLocalTrack *tp = dcAna->GetTrackSdcIn( it );
    if( tp ) gEvDisp.DrawSdcInLocalTrack( tp );
  }

  event.TrackSearchSdcOut();
  Int_t ntSdcOut = dcAna->GetNtracksSdcOut();
  for( Int_t it=0; it<ntSdcOut; ++it ){
    DCLocalTrack *tp = dcAna->GetTrackSdcOut( it );
    if( tp ) gEvDisp.DrawSdcOutLocalTrack( tp );
  }

  event.TrackSearchKurama();
  //  Int_t ntKurama = dcAna->GetNTracksKurama();

  gEvDisp.Update();

  gEvDisp.EndOfEvent();

  return 0;
}

}
