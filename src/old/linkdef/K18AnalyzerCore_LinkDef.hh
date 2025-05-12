// -*- C++ -*-

#include <vector>

#if defined(__MAKECINT__) || defined(__MAKECLING__)

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace hddaq;
// #pragma link C++ namespace hddaq::gui;
#pragma link C++ namespace hddaq::unpacker;
#pragma link C++ namespace analyzer;
#pragma link C++ enum DetectorType;
#pragma link C++ enum SubDetectorType;
#pragma link C++ enum DataType;
#pragma link C++ enum ParticleType;

#pragma link C++ class BH2Cluster+;
#pragma link C++ class BH2Filter+;
#pragma link C++ class BH2Hit+;
#pragma link C++ class DCAnalyzer+;
#pragma link C++ class DCHit+;
#pragma link C++ class DCPairHitCluster+;
#pragma link C++ class DCRawHit+;
#pragma link C++ class DCLocalTrack+;
#pragma link C++ class DCLTrackHit+;
#pragma link C++ class FiberCluster+;
#pragma link C++ class FiberHit+;
#pragma link C++ class FieldElements+;
#pragma link C++ class FLHit+;
#pragma link C++ class HodoAnalyzer+;
#pragma link C++ class HodoCluster+;
#pragma link C++ class Hodo1Hit+;
#pragma link C++ class Hodo2Hit+;
#pragma link C++ class HodoRawHit+;
#pragma link C++ class K18TrackD2U+;
#pragma link C++ class KuramaFieldMap+;
#pragma link C++ class KuramaTrack+;
#pragma link C++ class RawData+;
// #pragma link C++ class RungeKutta+;
// #pragma link C++ class SsdCluster+;
#pragma link C++ class TrackHit+;

#pragma link C++ class ConfMan+;
#pragma link C++ class DCDriftParamMan+;
#pragma link C++ class DCGeomMan+;
#pragma link C++ class DCGeomRecord+;
#pragma link C++ class DCTdcCalibMan+;
#pragma link C++ class EMCParamMan+;
#pragma link C++ class EventDisplay+;
#pragma link C++ class Exception+;
#pragma link C++ class FieldMan+;
#pragma link C++ class GeAdcCalibMan+;
#pragma link C++ class HistMaker+;
#pragma link C++ class HodoParamMan+;
#pragma link C++ class HodoPHCMan+;
#pragma link C++ class HttpServer+;
#pragma link C++ class K18TransMatrix+;
#pragma link C++ class MatrixParamMan+;
#pragma link C++ class MsTParamMan+;
#pragma link C++ class ScalerAnalyzer+;
#pragma link C++ class UserParamMan+;

#pragma link C++ class analyzer::EventAnalyzer+;
#pragma link C++ class debug::ObjectCounter+;
#pragma link C++ class debug::Timer+;

#pragma link C++ class std::vector< std::vector<Int_t> >+;
#pragma link C++ class std::vector< std::vector<Double_t> >+;
#pragma link C++ class std::vector< TString >+;
#pragma link C++ class std::vector< std::vector< std::vector<Double_t> > >+;

#endif
