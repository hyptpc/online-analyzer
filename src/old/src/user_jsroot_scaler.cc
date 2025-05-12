// -*- C++ -*-

// Author: Shuhei Hayakawa

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGFileBrowser.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>

#include <THttpServer.h>
#include <TKey.h>
#include <TSystem.h>
#include <TTimeStamp.h>

#include <DAQNode.hh>
#include <filesystem_util.hh>
#include <Unpacker.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <std_ostream.hh>

#include "user_analyzer.hh"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "HttpServer.hh"
#include "ScalerAnalyzer.hh"

#define CFT 0

namespace analyzer
{
using namespace hddaq::unpacker;
using namespace hddaq;

namespace
{
HttpServer& gHttp = HttpServer::GetInstance();
ScalerAnalyzer scaler_on;
ScalerAnalyzer scaler_off;
const std::chrono::milliseconds flush_interval(100);
}

//____________________________________________________________________________
Int_t
process_begin(const std::vector<std::string>& argv)
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  if(!gConfMan.IsGood()) return -1;

  gHttp.SetPort(9092);
  gHttp.Open();
  gHttp.SetItemField("/", "_monitoring", "100");
  gHttp.SetItemField("/", "_layout", "vert3");
  gHttp.SetItemField("/", "_drawitem", "[ScalerOn,ScalerOff,Tag]");
  gHttp.CreateItem("/ScalerOn", "DAQ Scaler On");
  gHttp.SetItemField("/ScalerOn", "_kind", "Text");
  gHttp.CreateItem("/ScalerOff", "DAQ Scaler Off");
  gHttp.SetItemField("/ScalerOff", "_kind", "Text");
  gHttp.CreateItem("/Tag", "Tag Check");
  gHttp.SetItemField("/Tag", "_kind", "Text");
  std::stringstream ss;
  ss << "<div style='color: white; background-color: black;"
     << "width: 100%; height: 100%;'>"
     << "Tag cheker is running" << "</div>";
  gHttp.SetItemField("/Tag", "value", ss.str().c_str());
  gHttp.Hide("/Reset");

  scaler_on.SetFlag(ScalerAnalyzer::kSeparateComma);
  scaler_on.SetFlag(ScalerAnalyzer::kSpillBySpill);
  scaler_on.SetFlag(ScalerAnalyzer::kSpillOn);

  scaler_off.SetFlag(ScalerAnalyzer::kSeparateComma);
  scaler_off.SetFlag(ScalerAnalyzer::kSpillBySpill);
  scaler_off.SetFlag(ScalerAnalyzer::kSpillOff);

  //////////////////// Set Channels
  // ScalerAnalylzer::Set(Int_t column,
  //                       Int_t raw,
  //                       ScalerInfo(name, module, channel));
  // scaler information is defined from here.
  // please do not use a white space character.
  {
    Int_t c = ScalerAnalyzer::kLeft;
    Int_t r = 0;
    scaler_on.Set(c, r++, ScalerInfo("BH1",        0, 16));
    scaler_on.Set(c, r++, ScalerInfo("BH1-SUM",   -1, -1));
    scaler_on.Set(c, r++, ScalerInfo("BH1-01",     1,  0));
    scaler_on.Set(c, r++, ScalerInfo("BH1-02",     1,  1));
    scaler_on.Set(c, r++, ScalerInfo("BH1-03",     1,  2));
    scaler_on.Set(c, r++, ScalerInfo("BH1-04",     1,  3));
    scaler_on.Set(c, r++, ScalerInfo("BH1-05",     1,  4));
    scaler_on.Set(c, r++, ScalerInfo("BH1-06",     1,  5));
    scaler_on.Set(c, r++, ScalerInfo("BH1-07",     1,  6));
    scaler_on.Set(c, r++, ScalerInfo("BH1-08",     1,  7));
    scaler_on.Set(c, r++, ScalerInfo("BH1-09",     1,  8));
    scaler_on.Set(c, r++, ScalerInfo("BH1-10",     1,  9));
    scaler_on.Set(c, r++, ScalerInfo("BH1-11",     1, 10));
    scaler_on.Set(c, r++, ScalerInfo("BH2",        0, 17));
    scaler_on.Set(c, r++, ScalerInfo("BH2-SUM",   -1, -1));
    scaler_on.Set(c, r++, ScalerInfo("BH2-01",     0, 64));
    scaler_on.Set(c, r++, ScalerInfo("BH2-02",     0, 65));
    scaler_on.Set(c, r++, ScalerInfo("BH2-03",     0, 66));
    scaler_on.Set(c, r++, ScalerInfo("BH2-04",     0, 67));
    scaler_on.Set(c, r++, ScalerInfo("BH2-05",     0, 68));
    scaler_on.Set(c, r++, ScalerInfo("BH2-06",     0, 69));
    scaler_on.Set(c, r++, ScalerInfo("BH2-07",     0, 70));
    scaler_on.Set(c, r++, ScalerInfo("BH2-08",     0, 71));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-1",    2, 2));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-2",    2, 3));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-3",    2, 4));
    scaler_on.Set( c, r++, ScalerInfo( "AFT-01",     0, 80 ) );
    scaler_on.Set( c, r++, ScalerInfo( "AFT-02",     0, 81 ) );
    scaler_on.Set( c, r++, ScalerInfo( "AFT-03",     0, 82 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kCenter;
    Int_t r = 0;
    scaler_on.Set(c, r++, ScalerInfo("10M-Clock",    0,  0));
    scaler_on.Set(c, r++, ScalerInfo("TM",           0,  9));
    scaler_on.Set(c, r++, ScalerInfo("SY",           0, 10));
    scaler_on.Set(c, r++, ScalerInfo("K-Beam",       0, 35));
    scaler_on.Set(c, r++, ScalerInfo("Pi-Beam",      0, 39));
    scaler_on.Set(c, r++, ScalerInfo("Beam",         0, 36));
    scaler_on.Set(c, r++, ScalerInfo("BH1-1/100-PS", 1, 11));
    scaler_on.Set(c, r++, ScalerInfo("BH1-1/1e5-PS", 1, 12));
    //    scaler_on.Set(c, r++, ScalerInfo("TOF-24",       0, 29));
    scaler_on.Set(c, r++, ScalerInfo("Mtx2D-1",      0, 32));
    scaler_on.Set(c, r++, ScalerInfo("Mtx2D-2",      0, 33));
    scaler_on.Set(c, r++, ScalerInfo("Mtx3D",        0, 34));
    // scaler_on.Set(c, r++, ScalerInfo("Other1",       0, 25));
    // scaler_on.Set(c, r++, ScalerInfo("Other2",       0, 26));
    scaler_on.Set(c, r++, ScalerInfo("Other3",       0, 27));
    scaler_on.Set(c, r++, ScalerInfo("Other4",       0, 28));
    // scaler_on.Set(c, r++, ScalerInfo("LSO1",         2, 48));
    // scaler_on.Set(c, r++, ScalerInfo("LSO2",         2, 49));
    // scaler_on.Set(c, r++, ScalerInfo("LSO1xGe13",    2, 50));
    // scaler_on.Set(c, r++, ScalerInfo("LSO2xGe24",    2, 51));
    // scaler_on.Set(c, r++, ScalerInfo("GeCoin1",      2, 52));
    // scaler_on.Set(c, r++, ScalerInfo("GeCoin2",      2, 53));
    // scaler_on.Set(c, r++, ScalerInfo("LSO1High",     2, 54));
    // scaler_on.Set(c, r++, ScalerInfo("LSO2High",     2, 55));
    // scaler_on.Set(c, r++, ScalerInfo("GeHigh1",      2, 56));
    // scaler_on.Set(c, r++, ScalerInfo("GeHigh2",      2, 57));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-A",       0, 35));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-B",       0, 36));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-C",       0, 37));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-D",       0, 38));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-E",       0, 39));
    scaler_on.Set(c, r++, ScalerInfo("BEAM-F",       0, 40));
    scaler_on.Set(c, r++, ScalerInfo("BAC",        0, 18));
    scaler_on.Set(c, r++, ScalerInfo("TOF",        0, 22));
    scaler_on.Set(c, r++, ScalerInfo("AC1",        0, 23));
    scaler_on.Set(c, r++, ScalerInfo("WC",         0, 24));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-4",    2, 5));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-5",    2, 6));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-6",    2, 7));
    scaler_on.Set( c, r++, ScalerInfo( "AFT-04",     0, 83 ) );
    scaler_on.Set( c, r++, ScalerInfo( "AFT-05",     0, 84 ) );
    scaler_on.Set( c, r++, ScalerInfo( "AFT-06",     0, 85 ) );
  }

  {
    Int_t c = ScalerAnalyzer::kRight;
    Int_t r = 0;
    scaler_on.Set(c, r++, ScalerInfo("Spill",        -1, -1));
    scaler_on.Set(c, r++, ScalerInfo("Real-Time",     0,  1));
    scaler_on.Set(c, r++, ScalerInfo("Live-Time",     0,  2));
    scaler_on.Set(c, r++, ScalerInfo("L1-Req",        0,  3));
    scaler_on.Set(c, r++, ScalerInfo("L1-Acc",        0,  4));
    // scaler_on.Set(c, r++, ScalerInfo("MstClr",        0,  5));
    // scaler_on.Set(c, r++, ScalerInfo("Clear",         0,  6));
    scaler_on.Set(c, r++, ScalerInfo("L2-Req",        0,  7));
    scaler_on.Set(c, r++, ScalerInfo("L2-Acc",        0,  8));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-A",        0, 41));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-B",        0, 42));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-C",        0, 43));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-D",        0, 44));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-E",        0, 45));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-F",        0, 46));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-A-PS",     0, 48));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-B-PS",     0, 49));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-C-PS",     0, 50));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-D-PS",     0, 51));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-E-PS",     0, 52));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-F-PS",     0, 53));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-PSOR-A",   0, 54));
    scaler_on.Set(c, r++, ScalerInfo("TRIG-PSOR-B",   0, 55));
    scaler_on.Set(c, r++, ScalerInfo("Clock-PS",      0, 56));
    scaler_on.Set(c, r++, ScalerInfo("Reserve2-PS",   0, 57));
    scaler_on.Set(c, r++, ScalerInfo("Level1-PS",     0, 58));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-7",    2, 8));
    scaler_on.Set(c, r++, ScalerInfo("V792gate-8",    2, 9));
    scaler_on.Set( c, r++, ScalerInfo( "AFT-07",      0, 86));
    scaler_on.Set( c, r++, ScalerInfo( "AFT-08",      0, 87));
    scaler_on.Set( c, r++, ScalerInfo( "AFT-09",      0, 88));
  }

  for(Int_t i=0; i<ScalerAnalyzer::MaxColumn; ++i){
    for(Int_t j=0; j<ScalerAnalyzer::MaxRow; ++j){
      scaler_off.Set(i, j, scaler_on.GetScalerInfo(i, j));
    }
  }

  scaler_on.PrintFlags();
  scaler_off.PrintFlags();

  return 0;
}

//____________________________________________________________________________
Int_t
process_end()
{
  return 0;
}

//____________________________________________________________________________
Int_t
process_event()
{
  static auto& gUnpacker = GUnpacker::get_instance();
  static auto* root = gUnpacker.get_root();
  static auto& gConfig = GConfig::get_instance();
  static const TString tout = gConfig.get_control_param("tout");

  Int_t run_number = root->get_run_number();
  Int_t event_number = gUnpacker.get_event_number();

  std::stringstream ss;

  // Tag
  ss.str("");
  ss << "<div style='color: white; background-color: black;"
     << "width: 100%; height: 100%;'>";
  ss << "RUN " << run_number << "   Event " << event_number
     << "<br>";
  if(!gUnpacker.is_good()){
    std::ifstream ifs(tout);
    if(ifs.good()){
      TString buf;
      while(!ifs.eof()){
	buf.ReadLine(ifs, false);
	if(buf.Contains("!") && !buf.Contains("............!"))
	  buf = "<font color='yellow'>" + buf + "</font>";
	// if(buf.Contains("bUuSELselthdg"))
	//   ss << TString(' ', 24);
	ss << buf << "<br>";
      }
      ifs.close();
      tag_summary.seekp(0, std::ios_base::beg);
    } else {
      ss << Form("Failed to read %s", tout.Data());
    }
    ss << "</div>";
    gHttp.SetItemField("/Tag", "value", ss.str().c_str());
  }

  const Int_t MaxDispRow = 29;

  // Scaler Spill On
  auto now = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::system_clock::now().time_since_epoch());
  static auto prev_flush = now;
  Bool_t flush_flag = (now - prev_flush) > flush_interval;
  prev_flush = now;

  if(scaler_on.Decode()){
    if(flush_flag && !scaler_on.IsSpillEnd())
      // if(!scaler_on.IsSpillEnd())
      return 0;

    ss.str("");
    ss << "<div style='color: white; background-color: black;"
       << "width: 100%; height: 100%;'>";
    ss << "<table border=\"0\" width=\"700\" cellpadding=\"0\">";
    TString end_mark = scaler_on.IsSpillEnd() ? "Spill End" : "";
    ss << "<tr><td width=\"100\">RUN</td><td align=\"right\" width=\"100\">"
       << scaler_on.SeparateComma(run_number) << "</td><td width=\"100\">"
       << " : Event Number" << "</td><td align=\"right\">"
       << scaler_on.SeparateComma(event_number) << "</td>"
       << "<td width=\"100\">" << " : " << "</td>"
       << "<td align=\"right\" width=\"100\">" << end_mark << "</td>"
       << "<tr><td></td><td></td><td></td></tr>";
    for(Int_t j=0; j<MaxDispRow; ++j){
      ss << "<tr>";
      for(Int_t i=0; i<ScalerAnalyzer::MaxColumn; ++i){
	// for(Int_t i=0; i<ScalerAnalyzer::MaxColumn; ++i){
        TString n = scaler_on.GetScalerName(i, j);
        // if(n.Contains("n/a"))
        //   continue;
	ss << "<td>";
	if(i != 0)
	  ss << " : ";
	ss << n << "</td>"
	   << "<td align=\"right\">"
           << scaler_on.SeparateComma(scaler_on.Get(i, j)) << "</td>";
      }
      ss << "</tr>";
    }
    ss << "<tr><td></td><td></td><td></td></tr>"
       << "<tr><td>K-Beam/TM</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Fraction("K-Beam", "TM"))
       << "</td>"
       << "<td> : Live/Real</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Fraction("Live-Time", "Real-Time"))
       << "</td>"
       << "<td> : DAQ-Eff</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Fraction("L1-Acc", "L1-Req"))
       << "</td>"
       << "</tr></tr>"
       << "<td>L1Req/K-Beam</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Fraction("L1-Req", "K-Beam"))
       << "</td>"
       << "<td> : L2-Eff</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Fraction("L2-Acc", "L1-Acc"))
       << "</td>"
       << "<td> : Duty-Factor</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_on.Duty())
       << "</td>"
       << "</tr>";
    ss << "</table>";
    ss << "</div>";
    gHttp.SetItemField("/ScalerOn", "value", ss.str().c_str());
  }

  // Scaler Spill Off
  if(scaler_off.Decode()){
    if(flush_flag && !scaler_off.IsSpillEnd())
      // if(!scaler_off.IsSpillEnd())
      return 0;

    ss.str("");
    ss << "<div style='color: white; background-color: black;"
       << "width: 100%; height: 100%;'>";
    ss << "<table border=\"0\" width=\"700\" cellpadding=\"0\">";
    TString end_mark = scaler_off.IsSpillEnd() ? "Spill End" : "";
    ss << "<tr><td width=\"100\">RUN</td><td align=\"right\" width=\"100\">"
       << scaler_off.SeparateComma(run_number) << "</td><td width=\"100\">"
       << " : Event Number" << "</td><td align=\"right\">"
       << scaler_off.SeparateComma(event_number) << "</td>"
       << "<td width=\"100\">" << " : " << "</td>"
       << "<td align=\"right\" width=\"100\">" << end_mark << "</td>"
       << "<tr><td></td><td></td><td></td></tr>";
    for(Int_t j=0; j<MaxDispRow; ++j){
      // for(Int_t j=0; j<ScalerAnalyzer::MaxRow; ++j){
      ss << "<tr>";
      for(Int_t i=0; i<ScalerAnalyzer::MaxColumn; ++i){
        TString n = scaler_off.GetScalerName(i, j);
        // if(n.Contains("n/a"))
        //   continue;
	ss << "<td>";
	if(i != 0)
	  ss << " : ";
	ss << n << "</td>"
	   << "<td align=\"right\">"
           << scaler_off.SeparateComma(scaler_off.Get(i, j)) << "</td>";
      }
      ss << "</tr>";
    }
    ss << "<tr><td></td><td></td><td></td></tr>"
       << "<tr><td>K-Beam/TM</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Fraction("K-Beam", "TM"))
       << "</td>"
       << "<td> : Live/Real</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Fraction("Live-Time", "Real-Time"))
       << "</td>"
       << "<td> : DAQ-Eff</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Fraction("L1-Acc", "L1-Req"))
       << "</td>"
       << "</tr></tr>"
       << "<td>L1Req/K-Beam</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Fraction("L1-Req", "K-Beam"))
       << "</td>"
       << "<td> : L2-Eff</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Fraction("L2-Acc", "L1-Acc"))
       << "</td>"
       << "<td> : Duty-Factor</td>"
       << "<td align=\"right\">"
       << Form("%.6f", scaler_off.Duty())
       << "</td>"
       << "</tr>";
    ss << "</table>";
    ss << "</div>";
    gHttp.SetItemField("/ScalerOff", "value", ss.str().c_str());
  }

  if (gUnpacker.is_good()){
    ss.str("");
    ss << "<div style='color: white; background-color: black;"
       << "width: 100%; height: 100%;'>"
       << "Tag cheker is running" << "</div>";
    gHttp.SetItemField("/Tag", "value", ss.str().c_str());
  }

  if(scaler_on.IsSpillEnd() || scaler_off.IsSpillEnd()){
    gSystem->Sleep(150);
    gSystem->ProcessEvents();
  }

  if (!gUnpacker.is_good()){
    static const TString host(gSystem->Getenv("HOSTNAME"));
    static auto prev_spill = scaler_on.Get("Spill");
    auto        curr_spill = scaler_on.Get("Spill");
    if(host.Contains("k18term4") &&
       event_number > 1 && prev_spill != curr_spill){
      std::cout << "exec tagslip sound!" << std::endl;
      gSystem->Exec("ssh k18epics.monitor.k18net \"aplay ~/sound/tagslip.wav\" &");
    }
    prev_spill = curr_spill;
  }
  gSystem->ProcessEvents();
  return 0;
}
}
