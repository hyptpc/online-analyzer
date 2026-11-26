// -*- C++ -*-

#include "TpcPadMapWeb.hh"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TGraph.h>
#include <TH1D.h>
#include <TH2Poly.h>
#include <TNamed.h>
#include <TSystem.h>

#include "HttpServer.hh"

namespace
{
const char* kFolder = "/TPCPadMap";

struct PadMap
{
  TH2Poly*                 poly = nullptr;
  TH1D*                    h    = nullptr;
  std::vector<TH2PolyBin*> bins;
};
std::vector<PadMap> g_maps;

//_____________________________________________________________________________
// polygons of the TH2Poly bins, in bin order, as
// {"n":N,"x":[[x0,..],..],"y":[[y0,..],..]} [mm], rounded to 0.01 mm
std::string
GeometryJson(TH2Poly* poly)
{
  std::ostringstream os, ox, oy;
  ox.setf(std::ios::fixed); ox.precision(2);
  oy.setf(std::ios::fixed); oy.precision(2);
  Int_t n = 0;
  for(auto obj : *poly->GetBins()){
    auto g = dynamic_cast<TGraph*>(static_cast<TH2PolyBin*>(obj)->GetPolygon());
    ox << (n ? ",[" : "[");
    oy << (n ? ",[" : "[");
    for(Int_t i=0, np=(g ? g->GetN() : 0); i<np; ++i){
      ox << (i ? "," : "") << g->GetX()[i];
      oy << (i ? "," : "") << g->GetY()[i];
    }
    ox << "]";
    oy << "]";
    ++n;
  }
  os << "{\"n\":" << n << ",\"x\":[" << ox.str() << "],\"y\":["
     << oy.str() << "]}";
  return os.str();
}

//_____________________________________________________________________________
// the page is web/tpc_padmap.html next to the analyzer sources
TString
WebDir()
{
  TString dir = gSystem->DirName(__FILE__);
  dir += "/../web";
  gSystem->ExpandPathName(dir);
  return dir;
}
}

namespace tpcweb
{
//_____________________________________________________________________________
void
Initialize(TH1* adc2d, TH1* tdc2d, TH1* hit2d, TH1* thre2d)
{
  auto& gHttp = HttpServer::GetInstance();
  const char* names[] = { "TPC_PadADC", "TPC_PadTDC", "TPC_PadHit",
                          "TPC_PadThreshold" };
  TH1* src[] = { adc2d, tdc2d, hit2d, thre2d };
  TH2Poly* geom = nullptr;
  for(Int_t i=0; i<4; ++i){
    auto poly = dynamic_cast<TH2Poly*>(src[i]);
    if(!poly || !poly->GetBins()) continue;
    PadMap m;
    m.poly = poly;
    for(auto obj : *poly->GetBins())
      m.bins.push_back(static_cast<TH2PolyBin*>(obj));
    const Int_t n = m.bins.size();
    m.h = new TH1D(names[i], poly->GetTitle(), n, 0, n);
    m.h->SetDirectory(nullptr);
    m.h->SetMinimum(poly->GetMinimumStored());
    m.h->SetMaximum(poly->GetMaximumStored());
    gHttp.Register(m.h, kFolder);
    g_maps.push_back(m);
    if(!geom) geom = poly;
  }
  if(!geom) return;
  // the four pad maps share the same bins (HistMaker::createTPC)
  // the polygons as JSON text in the title
  auto str = new TNamed("TPC_PadGeom", GeometryJson(geom).c_str());
  gHttp.Register(str, kFolder);
  gHttp.Hide(TString(kFolder) + "/TPC_PadGeom");

  const TString dir = WebDir();
  if(gSystem->AccessPathName(dir + "/tpc_padmap.html")){
    std::cout << "#W tpcweb::Initialize() " << dir
              << "/tpc_padmap.html not found" << std::endl;
    return;
  }
  gHttp.AddLocation("tpcmap/", dir);
  std::cout << "#D tpcweb::Initialize() light TPC pad maps at"
            << " /tpcmap/tpc_padmap.html" << std::endl;
}

//_____________________________________________________________________________
void
Update()
{
  for(auto& m : g_maps){
    Double_t* a = m.h->GetArray();
    for(std::size_t i=0, n=m.bins.size(); i<n; ++i)
      a[i+1] = m.bins[i]->GetContent();
    m.h->SetTitle(m.poly->GetTitle());
  }
}
}
