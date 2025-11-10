// -*- C++ -*-

#include "MacroBuilder.hh"

#include <string>

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMacro.h>
#include <TString.h>
#include <TText.h>
#include <TArc.h>
#include <TLine.h>
#include <TPolyLine.h>

#include "DetectorID.hh"
#include "Main.hh"
#include "HistHelper.hh"
#include "HistMaker.hh"

namespace analyzer
{

namespace macro
{

//____________________________________________________________________________
TObject*
Get( TString name )
{
  std::string process = Main::getInstance().getArgv().at(0);
  Int_t n = process.find("bin");
  TString path = process.substr(0, n)+"src/analyzer/macro/";

  if( name.Contains(".C") )
    name.ReplaceAll(".C","");

  return new TMacro(path+name+".C");
}

}

//____________________________________________________________________________
// For HttpServer
namespace http
{

namespace
{
enum eUorD { kU, kD, kUorD };
enum eAorT { kA, kT, kAorT };
}

//____________________________________________________________________________
TCanvas*
BAC()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(3, 2);
  for( Int_t i=0; i<6; ++i ){
    c1->cd(i+1)->SetLogy();
    if (i<5) {
      TH1 *h = GHist::get( HistMaker::getUniqueID(kBAC, 0, kADC, i + 1 ) );
      if( !h ) continue;
      // h->SetMinimum(0.1);
      h->Draw();
      // h->GetXaxis()->SetRangeUser(xmin,xmax);
      h = GHist::get( HistMaker::getUniqueID(kBAC, 0, kADCwTDC, i + 1 ) );
      if( !h ) continue;
      h->SetLineColor(kRed+1);
      h->Draw("same");
    } else {
      TH1 *h = GHist::get(HistMaker::getUniqueID(kBAC, 0, kTDC, 5));
      if( !h ) continue;
      h->Draw();
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
KVC1()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 4);
  for(Int_t ch=0; ch<4; ++ch){
    for(Int_t seg=0; seg<4; ++seg){
      c1->cd(seg+ch*4+1)->SetLogy();
      if(ch<3){
	TH1 *h = GHist::get( HistMaker::getUniqueID(kKVC1, ch, kADC, seg+1));
	if(!h) continue;
	h->Draw();
	h = GHist::get( HistMaker::getUniqueID(kKVC1, ch, kADCwTDC, seg+1));
	if(!h) continue;
	h->SetLineColor(kRed+1);
	h->Draw("same");
      } else {
	TH1 *h = GHist::get(HistMaker::getUniqueID(kKVC1, 2, kTDC, seg+1));
	if( !h ) continue;
	h->Draw();
      }
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
KVC2()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  const Int_t n_seg = 4;
  const Int_t n_ch = 5;
  c1->Divide(n_seg, n_ch+1);
  for(Int_t ch=0; ch<n_ch+1; ++ch){
    for(Int_t seg=0; seg<n_seg; ++seg){
      c1->cd(seg+ch*n_seg+1)->SetLogy();
      if(ch<n_ch){
	TH1 *h = GHist::get( HistMaker::getUniqueID(kKVC2, ch, kADC, seg+1));
	if(!h) continue;
	h->Draw();
	h = GHist::get( HistMaker::getUniqueID(kKVC2, ch, kADCwTDC, seg+1));
	if(!h) continue;
	h->SetLineColor(kRed+1);
	h->Draw("same");
      } else {
	TH1 *h = GHist::get(HistMaker::getUniqueID(kKVC2, 4, kTDC, seg+1));
	if( !h ) continue;
	h->Draw();
      }
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
KVC()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  const Int_t n_seg = 8;
  const Int_t n_ch = 5;
  c1->Divide(n_seg, n_ch+1);
  for(Int_t ch=0; ch<n_ch+1; ++ch){
    for(Int_t seg=0; seg<n_seg; ++seg){
      c1->cd(seg+ch*n_seg+1)->SetLogy();
      if(ch<n_ch){
	TH1 *h = GHist::get( HistMaker::getUniqueID(kKVC, ch, kADC, seg+1));
	if(!h) continue;
	h->Draw();
	h = GHist::get( HistMaker::getUniqueID(kKVC, ch, kADCwTDC, seg+1));
	if(!h) continue;
	h->SetLineColor(kRed+1);
	h->Draw("same");
      } else {
	TH1 *h = GHist::get(HistMaker::getUniqueID(kKVC, 4, kTDC, seg+1));
	if( !h ) continue;
	h->Draw();
      }
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
SAC()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 3);
  for( Int_t i=0; i<10; ++i ){
    c1->cd(i+1)->SetLogy();
    if (i<9) {
      TH1 *h = GHist::get( HistMaker::getUniqueID(kSAC, 0, kADC, i + 1 ) );
      if( !h ) continue;
      h->Draw();
      h = GHist::get( HistMaker::getUniqueID(kSAC, 0, kADCwTDC, i + 1 ) );
      if( !h ) continue;
      h->SetLineColor(kRed+1);
      h->Draw("same");
    } else {
      TH1 *h = GHist::get(HistMaker::getUniqueID(kSAC, 0, kTDC, 9));
      if( !h ) continue;
      h->Draw();
    }
  }

  return c1;
}

//____________________________________________________________________________
TCanvas*
BH2ADC(Int_t ch)
{
  TString UorD = (ch == 0) ? "U" : "D";
  TCanvas *c1 = new TCanvas(__func__ + UorD, __func__ + UorD);
  c1->Divide(5, 3);
  const Int_t n_seg = 15;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kBH2, ch, kADC, seg+1));
    if(!h) continue;
    h->Draw();
    h = GHist::get( HistMaker::getUniqueID(kBH2, ch, kADCwTDC, seg+1));
    if(!h) continue;
    h->SetLineColor(kRed+1);
    h->Draw("same");
  }
  // for(Int_t ch=0; ch<3; ++ch){
  //   c1->cd(ch+4)->SetLogy();
  //   TH1 *h = GHist::get(HistMaker::getUniqueID(kBH2, ch, kTDC, 1));
  //   if( !h ) continue;
  //   h->Draw();
  // }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BH2TDC(Int_t ch)
{
  TString UorD = (ch == 0) ? "U" : (ch == 1 ) ? "D" : "MT";
  TCanvas *c1 = new TCanvas(__func__ + UorD, __func__ + UorD);
  c1->Divide(5, 3);
  const Int_t n_seg = 15;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kBH2, ch, kTDC, seg+1));
    if(!h) continue;
    h->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BH2ADCU()
{
  return BH2ADC(0);
}

//____________________________________________________________________________
TCanvas*
BH2ADCD()
{
  return BH2ADC(1);
}

//____________________________________________________________________________
TCanvas*
BH2TDCU()
{
  return BH2TDC(0);
}

//____________________________________________________________________________
TCanvas*
BH2TDCD()
{
  return BH2TDC(1);
}

//____________________________________________________________________________
TCanvas*
BH2TDCMT()
{
  return BH2TDC(2);
}

//____________________________________________________________________________
TCanvas*
HTOFADC(Int_t ch)
{
  TString UorD = (ch == 0) ? "U" : (ch == 1 ) ? "D" : "SUM";
  TCanvas *c1 = new TCanvas(__func__ + UorD, __func__ + UorD);
  c1->Divide(6, 6);
  const Int_t n_seg = 34;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kHTOF, ch, kADC, seg+1));
    if(!h) continue;
    h->Draw();
    h = GHist::get( HistMaker::getUniqueID(kHTOF, ch, kADCwTDC, seg+1));
    if(!h) continue;
    h->SetLineColor(kRed+1);
    h->Draw("same");
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
HTOFTDC(Int_t ch)
{
  TString UorD = (ch == 0) ? "U" : (ch == 1 ) ? "D" : "HT";
  TCanvas *c1 = new TCanvas(__func__ + UorD, __func__ + UorD);
  c1->Divide(6, 6);
  const Int_t n_seg = 34;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kHTOF, ch, kTDC, seg+1));
    if(!h) continue;
    h->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
HTOFADCU()
{
  return HTOFADC(0);
}

//____________________________________________________________________________
TCanvas*
HTOFADCD()
{
  return HTOFADC(1);
}

//____________________________________________________________________________
TCanvas*
HTOFADCSUM()
{
  return HTOFADC(2);
}

//____________________________________________________________________________
TCanvas*
HTOFTDCU()
{
  return HTOFTDC(0);
}

//____________________________________________________________________________
TCanvas*
HTOFTDCD()
{
  return HTOFTDC(1);
}

//____________________________________________________________________________
TCanvas*
HTOFTDCHT()
{
  return HTOFTDC(2);
}

//____________________________________________________________________________
TCanvas*
HTOFThreshold()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(6, 6);
  const Int_t n_seg = 34;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1);
    auto h = GHist::get( HistMaker::getUniqueID(kHTOF, 0, kThreshold, seg+1));
    if(!h) continue;
    h->SetMinimum(0.0);
    h->SetMaximum(1.5);
    h->SetLineWidth(2);
    h->Draw("hist");
  }
  return c1;
}

//_____________________________________________________________________________
TCanvas*
HTOFTDC2()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 5);
  Int_t n_seg = 4;
  Int_t n_ch = 5;
  for(Int_t ch=0; ch<n_ch; ++ch){
    for(Int_t i=0; i<n_seg; ++i){
      Int_t seg = i+18;
      c1->cd(i+ch*n_seg+1)->SetLogy();
      auto h = GHist::get(HistMaker::getUniqueID(kHTOF, ch, kTDC, seg+1));
      if(!h) continue;
      h->Draw();
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
T1()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(2, 1);
  c1->cd(1)->SetLogy();
  TH1 *h = GHist::get( HistMaker::getUniqueID(kT1, 0, kADC, 1 ) );
  if( !h ) return c1;
  h->Draw();
  h = GHist::get( HistMaker::getUniqueID(kT1, 0, kADCwTDC, 1 ) );
  if( !h ) return c1;
  h->SetLineColor(kRed+1);
  h->Draw("same");

  c1->cd(2)->SetLogy();
  h = GHist::get( HistMaker::getUniqueID(kT1, 0, kTDC, 1 ) );
  if( !h ) return c1;
  h->Draw();
  
  return c1;
}

//____________________________________________________________________________
TCanvas*
CVCADC()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 4);
  const Int_t n_seg = 8;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kCVC, 0, kADC, seg+1));
    if(!h) continue;
    h->Draw();
    h = GHist::get( HistMaker::getUniqueID(kCVC, 0, kADCwTDC, seg+1));
    if(!h) continue;
    h->SetLineColor(kRed+1);
    h->Draw("same");

    c1->cd(seg+n_seg+1)->SetLogy();
    h = GHist::get(HistMaker::getUniqueID(kCVC, 1, kADC, seg+1));
    if(!h) continue;
    h->Draw();
    h = GHist::get( HistMaker::getUniqueID(kCVC, 1, kADCwTDC, seg+1));
    if(!h) continue;
    h->SetLineColor(kRed+1);
    h->Draw("same");
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
CVCTDC()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 6);
  const Int_t n_seg = 8;
  for(Int_t seg=0; seg<n_seg; ++seg){
    c1->cd(seg+1)->SetLogy();
    auto h = GHist::get(HistMaker::getUniqueID(kCVC, 0, kTDC, seg+1));
    if(!h) continue;
    h->Draw();

    c1->cd(seg+n_seg+1)->SetLogy();
    h = GHist::get(HistMaker::getUniqueID(kCVC, 1, kTDC, seg+1));
    if(!h) continue;
    h->Draw();

    c1->cd(seg+2*n_seg+1)->SetLogy();
    h = GHist::get(HistMaker::getUniqueID(kCVC, 2, kTDC, seg+1));
    if(!h) continue;
    h->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
SAC3()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(2, 1);
  c1->cd(1)->SetLogy();
  TH1 *h = GHist::get( HistMaker::getUniqueID(kSAC3, 0, kADC, 1 ) );
  if( !h ) return c1;
  h->Draw();
  h = GHist::get( HistMaker::getUniqueID(kSAC3, 0, kADCwTDC, 1 ) );
  if( !h ) return c1;
  h->SetLineColor(kRed+1);
  h->Draw("same");

  c1->cd(2)->SetLogy();
  h = GHist::get( HistMaker::getUniqueID(kSAC3, 0, kTDC, 1 ) );
  if( !h ) return c1;
  h->Draw();
  
  return c1;
}

//____________________________________________________________________________
TCanvas*
SFV()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(3, 2);
  for( Int_t i=0; i<5; ++i ){
    c1->cd(i+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(kSFV, 0, kTDC, i + 1 ) );
    if( !h ) continue;
    h->Draw();
  }  
  return c1;
}

//____________________________________________________________________________
TCanvas*
COBO()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 2);
  for( Int_t i=0; i<8; ++i ){
    c1->cd(i+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(kCOBO, 0, kTDC, i + 1 ) );
    if( !h ) continue;
    h->Draw();
  }  
  return c1;
}


//_____________________________________________________________________________
TCanvas*
BVHTDCTOT()
{
  auto c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 2);
  for (Int_t seg=0; seg<NumOfSegBVH; ++seg) {
    c1->cd(seg+1); // ->SetLogy();
    auto h1 = GHist::get(HistMaker::getUniqueID(kBVH, 0, kTDC, seg+1));
    if(!h1) continue;
    h1->Draw();
    c1->cd(seg+1+NumOfSegBVH); // ->SetLogy();
    auto h2 = GHist::get(HistMaker::getUniqueID(kBVH, 0, kADC, seg+1));
    if(!h2) continue;
    h2->Draw();
  }
  return c1;
}

//_____________________________________________________________________________
TCanvas*
T1T2()
{
  auto c1 = new TCanvas(__func__, __func__);
  c1->Divide(4, 2);
  const Int_t n_ch = 2;
  std::vector<Int_t> device_id = { kT1, kT2 };
  for (Int_t i=0, n=device_id.size(); i<n; ++i) {
    for (Int_t ch=0; ch<n_ch; ++ch) {
      c1->cd(i*n_ch+ch+1)->SetLogy();
      auto h1 = GHist::get(HistMaker::getUniqueID(device_id[i], 0, kADC, ch+1));
      if(!h1) continue;
      h1->Draw();
      auto h2 = GHist::get(HistMaker::getUniqueID(device_id[i], 0, kADCwTDC, ch+1));
      if(!h2) continue;
      h2->SetLineColor(kRed+1);
      h2->Draw("same");
    }
  }
  for (Int_t i=0, n=device_id.size(); i<n; ++i) {
    for (Int_t ch=0; ch<n_ch; ++ch) {
      c1->cd(i*n_ch+ch+5)->SetLogy();
      auto h1 = GHist::get(HistMaker::getUniqueID(device_id[i], 0, kTDC, ch+1));
      if(!h1) continue;
      h1->Draw();
    }
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
Multiplicity()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  std::vector<UInt_t> hid_list = {
    HistMaker::getUniqueID(kBHT, 0, kMulti, 0),
    HistMaker::getUniqueID(kT0, 0, kMulti, 0),
    HistMaker::getUniqueID(kBH2, 0, kMulti, 0),
    HistMaker::getUniqueID(kBAC, 0, kMulti, 0),
    HistMaker::getUniqueID(kHTOF, 0, kMulti, 0),
    HistMaker::getUniqueID(kKVC, 0, kMulti, 0),
    HistMaker::getUniqueID(kKVC, 0, kMulti, 1),
    HistMaker::getUniqueID(kCVC, 0, kMulti, 0),
  };
  c1->Divide(4, 2);
  for(Int_t i=0, n=hid_list.size(); i<n; ++i){
    c1->cd(i+1);
    auto h1 = GHist::get(hid_list[i]);
    h1->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
HitPat()
{
  TCanvas *c1 = new TCanvas(__func__, __func__);
  std::vector<UInt_t> hid_list = {
    HistMaker::getUniqueID(kBHT, 0, kHitPat, 0),
    HistMaker::getUniqueID(kT0, 0, kHitPat, 0),
    HistMaker::getUniqueID(kBH2, 0, kHitPat, 0),
    HistMaker::getUniqueID(kBAC, 0, kHitPat, 0),
    HistMaker::getUniqueID(kHTOF, 0, kHitPat, 0),
    HistMaker::getUniqueID(kHTOF, 0, kHitPat, 1),
    HistMaker::getUniqueID(kHTOF, 0, kHitPat, 2),
    HistMaker::getUniqueID(kHTOF, 0, kHitPat, 3),
    HistMaker::getUniqueID(kKVC, 0, kHitPat, 0),
    HistMaker::getUniqueID(kT1, 0, kHitPat, 0),
    HistMaker::getUniqueID(kCVC, 0, kHitPat, 0),
    HistMaker::getUniqueID(kSAC3, 0, kHitPat, 0),
    HistMaker::getUniqueID(kSFV, 0, kHitPat, 0),
  };
  c1->Divide(4, 4);
  for(Int_t i=0, n=hid_list.size(); i<n; ++i){
    c1->cd(i+1);
    auto h1 = GHist::get(hid_list[i]);
    h1->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BLDCTDC( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kTDC, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
    if( det==kCDC2 ) continue;
    h = GHist::get( HistMaker::getUniqueID(det, 0, kTDC2D, i+1) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->Draw("hist same");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BHTTOT( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, init, kTOT, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BHTTDCvsTOT( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, init, kADC2D, i+1) );
    if( !h ) continue;
    gPad->SetLogz();
    h->Draw("col");
    //    h->GetXaxis()->SetRangeUser(500,1500);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCTOT( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kTOT, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BLDCTDCvsTOT( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kADC2D, i+1) );
    if( !h ) continue;
    gPad->SetLogz();
    h->Draw("col");
    h->GetXaxis()->SetRangeUser(500,1500);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCHitPat( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kHitPat, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BLDCResid( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kResid, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCXYProf( DetectorType det, std::string strDet,int offs)
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(3,3);
  double range=100;
  //  if(det==kSDC) range=5;
  for( Int_t i=0; i<9; ++i ){
    c1->cd(i+1);
    TH2 *h = (TH2*)GHist::get( HistMaker::getUniqueID(det, 0, kProf, i+1+offs) );
    if( !h ) continue;
    if(i%3==0){
      h->Draw("colz");
      h->GetXaxis()->SetRangeUser(-range,range);
      h->GetYaxis()->SetRangeUser(-range,range);
      if(det==kBPC){
	TArc arc;
	arc.SetLineColor(2);
	arc.SetLineWidth(2);
	arc.SetFillStyle(0);
	arc.DrawArc(0,4,34);
	// arc.DrawArc(0,0,34);
      }
    }else{
      h->Draw();
      h->GetXaxis()->SetRangeUser(-range,range);
    }
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCProf( DetectorType det, std::string strDet)
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(2,2);
  double range=100;
  //  if(det==kSDC) range=5;
  int ilist[4]={1,31,32,33};
  for( Int_t i=0; i<4; ++i ){
    c1->cd(i+1);
    TH2 *h = (TH2*)GHist::get( HistMaker::getUniqueID(det, 0, kProf, ilist[i]) );
    if( !h ) continue;
    h->Draw("colz");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BeamAxis()
{
  TString cname=__func__;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(3,3);
  double range=100;
  int detlist[2]={kBLC1a,kBLC2a};
  for( Int_t i=0; i<2; ++i ){
    c1->cd(3*i+1);
    TH2 *h2 = (TH2*)GHist::get( HistMaker::getUniqueID(detlist[i], 0, kProf, 1) );
    if( !h2 ) continue;
    h2->Draw("colz");
    TH1* h1 = (TH1*)GHist::get( HistMaker::getUniqueID(detlist[i], 0, kProf, 2) );
    if( !h1 ) continue;
    c1->cd(3*i+2);
    h1->Draw();
    h1 = (TH1*)GHist::get( HistMaker::getUniqueID(detlist[i], 0, kProf, 10) );
    if( !h1 ) continue;
    c1->cd(3*i+3);
    h1->Draw();
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
BLDCCorr(int ichm, std::string strDet)
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(2,2);
  double range=100;    
  for( Int_t i=0; i<4; ++i ){
    c1->cd(i+1);
    TH2 *h2 = (TH2*)GHist::get( HistMaker::getUniqueID(kAna, 0, 40, ichm*10+1+i) );
    if( !h2 ) continue;
    h2->Draw("colz");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCMulti( DetectorType det, std::string strDet, int nlayers, int init, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  TText tex;
  c1->Divide(nx,ny);
  for( Int_t i=init; i<init+nlayers; ++i ){
    c1->cd(i-init+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kMulti, i+1) );
    if( !h ) continue;
    //    h->GetXaxis()->SetRangeUser(0,2000);
    h->SetMinimum(0.1);
    h->Draw("hist");
    h->SetMinimum(0.1);
    h = GHist::get( HistMaker::getUniqueID(det, 0, kMulti2D, i+1) );
    if( !h ) continue;
    //    h->GetXaxis()->SetRangeUser(0,2000);
    h->SetLineColor(2);
    h->Draw("histsame");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCEff( DetectorType det, std::string strDet )
{
  TString cname=__func__;
  TCanvas *c1 = new TCanvas(cname, cname);
  TText tex;
  c1->Divide(2,2);
  c1->cd(1);
  TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kEff, 0) );
  if( h )   h->Draw("hist");
  return c1;
}
//____________________________________________________________________________
TCanvas*
MHTDCMeanTime( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kCTime, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
MHTDCTDC( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, subDet, kTDC, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
    h = GHist::get( HistMaker::getUniqueID(det, subDet, kTDC2D, i+1) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->Draw("hist same");
  }
  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPC( void )
{
  std::vector<Int_t> id = {
    HistMaker::getUniqueID( kTPC, 0, kADC ),
    HistMaker::getUniqueID( kTPC, 0, kTDC ),
    HistMaker::getUniqueID( kTPC, 0, kPede ),
    HistMaker::getUniqueID( kTPC, 0, kMulti )
  };

  auto c1 = new TCanvas( __func__, __func__ );
  c1->Divide( 2, 2 );
  for( Int_t i=0, n=id.size(); i<n; ++i ){
    c1->cd( i+1 );
    auto h = GHist::get( id[i] );
    if( h ) h->Draw( "colz" );
  }
  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPC2D( void )
{
  std::vector<Int_t> id = {
    HistMaker::getUniqueID( kTPC, 0, kADC2D, 3 ),
    HistMaker::getUniqueID( kTPC, 0, kADC2D, 4 ),
    HistMaker::getUniqueID( kTPC, 0, kFADC ),
    HistMaker::getUniqueID( kTPC, 2, kTDC )
  };

  auto c1 = new TCanvas( __func__, __func__ );
  c1->Divide( 2, 2 );
  for( Int_t i=0, n=id.size(); i<n; ++i ){
    c1->cd( i+1 );
    if( i==1 || i==2 ) c1->cd( i+1 )->SetLogz();
    auto h = GHist::get( id[i] );
    if( h ) h->Draw( "colz" );
  }
  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPC3D( void )
{
  std::vector<Int_t> id = {
    HistMaker::getUniqueID( kTPC, 2, kTDC2D ),
    HistMaker::getUniqueID( kTPC, 3, kTDC2D ),
    HistMaker::getUniqueID( kTPC, 0, kADC2D ),
    HistMaker::getUniqueID( kTPC, 0, kTDC2D ),
    HistMaker::getUniqueID( kTPC, 1, kTDC2D ),
    HistMaker::getUniqueID( kTPC, 0, kADC2D, 4 )
  };

  auto c1 = new TCanvas( __func__, __func__ );
  c1->Divide( 3, 2 );
  for( Int_t i=0, n=id.size(); i<n; ++i ){
    c1->cd( i+1 );
    auto h = GHist::get( id[i] );
    if( h ) h->Draw( "colz" );
  }
  return c1;
}
//_____________________________________________________________________________
TCanvas*
TPCADCPAD( void )
{
  auto id = HistMaker::getUniqueID( kTPC, 0, kADC2D );
  auto c1 = new TCanvas( __func__, __func__ );
  c1->cd()->SetLogz();
  auto h = GHist::get( id );
  if( h ) h->Draw( "colz" );

  Double_t l = (500./2.)/(1+sqrt(2.));
  Double_t px[9]={-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),
    l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),
    -l*(1+sqrt(2.))};
  Double_t py[9]={l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,
    -l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,
    l};
  TPolyLine* pLine = new TPolyLine( 9, px, py );
  pLine->SetLineColor(1);
  pLine->SetFillColorAlpha(kWhite, 0);
  pLine->Draw();

  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPCHTOFPAD( void )
{
  auto id = HistMaker::getUniqueID( kTPC, 0, kADC2D, 3 );
  auto id_htof = HistMaker::getUniqueID( kHTOF, 0, kHitPat, 2 );
  auto c1 = new TCanvas( __func__, __func__ );
  c1->cd();
  
  auto h_htof = GHist::get( id_htof );
  if( h_htof ){
    h_htof->SetMaximum( 200 );
    h_htof->Draw( "colz same" );
  }
  else 
  {std::cout << " no h_tof " << std::endl; getchar();}
  
  auto h = GHist::get( id );
  if( h ){
    h->SetLineWidth( 0 );
    h->GetXaxis()->SetRangeUser(-400,400);
    h->GetYaxis()->SetRangeUser(-400,400);
    h->SetMaximum( 200 );
    h->Draw( "col same" );
  }

  Double_t l = (500./2.)/(1+sqrt(2.));
  Double_t px[9]={-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),
    l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),
    -l*(1+sqrt(2.))};
  Double_t py[9]={l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,
    -l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,
    l};
  TPolyLine* pLine = new TPolyLine( 9, px, py );
  pLine->SetLineColor(1);
  pLine->SetFillColorAlpha(kWhite, 0);
  pLine->Draw();

  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPCTDCPAD( void )
{
  auto id = HistMaker::getUniqueID( kTPC, 0, kADC2D, 3 );
  auto c1 = new TCanvas( __func__, __func__ );
  c1->cd();
  auto h = GHist::get( id );
  if( h ) h->Draw( "colz" );

  Double_t l = (500./2.)/(1+sqrt(2.));
  Double_t px[9]={-l*(1+sqrt(2.)),-l,l,l*(1+sqrt(2.)),
    l*(1+sqrt(2.)),l,-l,-l*(1+sqrt(2.)),
    -l*(1+sqrt(2.))};
  Double_t py[9]={l,l*(1+sqrt(2.)),l*(1+sqrt(2.)),l,
    -l,-l*(1+sqrt(2.)),-l*(1+sqrt(2.)),-l,
    l};
  TPolyLine* pLine = new TPolyLine( 9, px, py );
  pLine->SetLineColor(1);
  pLine->SetFillColorAlpha(kWhite, 0);
  pLine->Draw();

  return c1;
}

//_____________________________________________________________________________
TCanvas*
TPCFADC( void )
{
  auto id = HistMaker::getUniqueID( kTPC, 0, kFADC );
  auto c1 = new TCanvas( __func__, __func__ );
  c1->cd()->SetLogz();
  auto h = GHist::get( id );
  if( h ) h->Draw( "colz" );
  return c1;
}

//____________________________________________________________________________
TCanvas*
TriggerFlagMHTDCTDC( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, subDet, kTDC, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
    h = GHist::get( HistMaker::getUniqueID(det, subDet, kTDC2D, i+1) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->Draw("hist same");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
SDDMHTDC( DetectorType det, std::string strDet, int subDet, int nlayers, int nx, int ny )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  const int nports=nlayers;
  const int nunits=4;
  int ipad=0;
  int type[2]={kTDC,kTDC2D};
  if(subDet==1) type[0]=kResetL,type[1]=kResetT;
  for( Int_t iport=0; iport<nports; ++iport ){
    for( Int_t iunit=0; iunit<nunits; ++iunit ){
      c1->cd(ipad+1);
      int index=iunit+iport*nunits;
      TH1 *h = GHist::get( HistMaker::getUniqueID(det, index,type[0], 1) );
      if( !h ) continue;
      h->Draw("hist");
      //    h->GetXaxis()->SetRangeUser(0,2000);
      h = GHist::get( HistMaker::getUniqueID(det, index, type[1], 1) );
      if( !h ) continue;
      h->SetLineColor(2);
      h->Draw("hist same");
      ipad++;
    }
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
BLDCWIRE( DetectorType det, std::string strDet, int layer, int nwires, int nx, int ny )
{
  TString cname=Form("%s_%s_layer%d",__func__,strDet.data(),layer);
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nwires; ++i ){
    c1->cd(i+1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, layer+1, kTDC, i+1) );
    if( !h ) continue;
    h->Draw("hist");
    //    h->GetXaxis()->SetRangeUser(0,2000);
    // h = GHist::get( HistMaker::getUniqueID(det, layer+1, kTDC2D, i+1) );
    // if( !h ) continue;
    // h->SetLineColor(2);
    // h->Draw("hist same");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
QDC( DetectorType det, std::string strDet, int subDet, int begin, int nlayers, int nx, int ny, double xmin, double xmax )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  for( Int_t i=0; i<nlayers; ++i ){
    c1->cd(i+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, subDet, kADC, begin + i + 1 ) );
    if( !h ) continue;
    h->SetMinimum(0.1);
    h->Draw("hist");
    // h->GetXaxis()->SetRangeUser(xmin,xmax);
    h = GHist::get( HistMaker::getUniqueID(det, subDet, kADCwTDC, begin + i + 1 ) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->Draw("hist same");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
PbG( DetectorType det, std::string strDet, int subDet, int begin, int nlayers, int nx, int ny, double xmin, double xmax )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(nx,ny);
  int ipad=0;
  for( Int_t i=0; i<nlayers; ++i ){
    ipad++;
    if(ipad==17||ipad==24||ipad==31) ipad+=3;
    c1->cd(ipad)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, subDet, kADC, begin + i + 1 ) );
    if( !h ) continue;
    h->SetMinimum(0.1);
    h->Draw("hist");
    h->GetXaxis()->SetRangeUser(xmin,xmax);
    h = GHist::get( HistMaker::getUniqueID(det, subDet, kADCwTDC, begin + i + 1 ) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->Draw("hist same");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
SDDADC( DetectorType det, std::string strDet, int iport, double xmin, double xmax, int rebin )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  int nsdds=8;
  int nunits=4;
  c1->Divide(nsdds,nunits);
  int ipad=0;
  for( Int_t iunit=0; iunit<nunits; ++iunit ){    
    for( Int_t i=0; i<nsdds; ++i ){    
      c1->cd(ipad+1)->SetLogy();
      int index=iunit+nunits*iport;
      TH1 *h = GHist::get( HistMaker::getUniqueID(det, index, kADC, i+1 ) );
      if( !h ) continue;
      h->Rebin(rebin);
      h->Draw("hist");
      h->SetMinimum(0.1);
      h->GetXaxis()->SetRangeUser(xmin,xmax);
      h = GHist::get( HistMaker::getUniqueID(det, index, kADCwTDC,  i + 1 ) );
      if( !h ) continue;
      h->SetLineColor(2);
      h->Rebin(rebin);
      h->Draw("hist same");
      ipad++;
    }
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
Check( int type, std::string strType)
{
  TString name=Form("Check_%s",strType.data());
  TCanvas *c1 = new TCanvas(name,name);
  c1->Divide(3,3);
  int detid[9]={kT0,kDEF,kCDH,kCDH,kCDH,kPbF2,kPbG,kVeto,kBTC};
  int seg[9]={2,2,0,16,32,14,16,0,0};
  int ud[9]={0,0,0,0,0,0,0,0,0};
  for(int i=0;i<9;i++){
    c1->cd(i+1)->SetLogy();
    TH1 *h = GHist::get( HistMaker::getUniqueID(detid[i], ud[i], kADC, seg[i] + 1 ) );
    if( !h ) continue;
    //    h->GetXaxis()->SetRangeUser(xmin,xmax);
    h->SetMinimum(0.1);
    h->Draw("hist");
    h = GHist::get( HistMaker::getUniqueID(detid[i], ud[i], kADCwTDC, seg[i] + 1 ) );
    if( !h ) continue;
    h->SetLineColor(2);
    h->SetMinimum(0.1);
    h->Draw("hist same");
  }
  return c1;
}

//____________________________________________________________________________
TCanvas*
Ana( int type )
{
  if(type==1){
    TCanvas *c1 = new TCanvas("BLC2a-T0","BLC2a-T0");
    TH1 *h;
    for(int i=0;i<5;i++){
      h = GHist::get( HistMaker::getUniqueID(kAna, 0, 2, i+1) );
      if( h )  { 
	h->SetLineColor(i+1);
	if(i==0) h->Draw("box");
	else     h->Draw("boxsame");
      }
    }
    return c1;
  }
  if(type==2){
    TCanvas *c1 = new TCanvas("BLC2b-T0","BLC2b-T0");
    TH1 *h;
    for(int i=0;i<5;i++){
      h = GHist::get( HistMaker::getUniqueID(kAna, 0, 3, i+1) );
      if( h )  { 
	h->SetLineColor(i+1);
	if(i==0) h->Draw("box");
	else     h->Draw("boxsame");
      }
    }
    return c1;
  }
  if(type==3){
    TCanvas *c1 = new TCanvas("BPC1-DEF","BLC1-DEF");
    TH1 *h;
    for(int i=0;i<5;i++){
      h = GHist::get( HistMaker::getUniqueID(kAna, 0, 4, i+1) );
      if( h )  { 
	h->SetLineColor(i+1);
	if(i==0) h->Draw("box");
	else     h->Draw("boxsame");
      }
    }
    return c1;
  }
  if(type==4){
    TCanvas *c1 = new TCanvas("BPC2-DEF","BPC2-DEF");
    TH1 *h;
    for(int i=0;i<5;i++){
      h = GHist::get( HistMaker::getUniqueID(kAna, 0, 5, i+1) );
      if( h )  { 
	h->SetLineColor(i+1);
	if(i==0) h->Draw("box");
	else     h->Draw("boxsame");
      }
    }
    return c1;
  }
  if(type==10){
    TCanvas *c1 = new TCanvas("BLDC_WireCorr","BLDC_WireCorr");
    c1->Divide(3,2);
    TH1 *h;
    const int nchm=4;
    DetectorType chm[nchm]={kBLC1a,kBLC1b, kBLC2a, kBLC2b};
    for(int i=0;i<nchm;i++){
      c1->cd(i+1);
      h = GHist::get( HistMaker::getUniqueID(chm[i], 0, kWireCorr, 0) );
      if( h )  { 
	h->Draw("box");
      }
    }
    return c1;
  }
}
//____________________________________________________________________________
TCanvas*
TOF( int type, std::string strTrig )
{
  TString name=Form("TOF_%s",strTrig.data());
  TCanvas *c1 = new TCanvas(name,name);
  TH1 *h;
  TString labels[6]={"All","ACHit","Kaon gate","Pion gate","Proton gate","Deuteron gate"};
  int color[6]={1,2,3,4,6,7};
  int width[6]={1,2,1,1,1,1};
  c1->Divide(1,2);
  TLegend *leg=new TLegend(0.1,0.6,0.3,0.9);
  for( int j=0;j<5;j++){
    h = GHist::get( HistMaker::getUniqueID(kAna, 0, 1, 20*type+j+1) ); // all
    h->SetLineColor(color[j]); 
    h->SetLineWidth(width[j]); 
    leg->AddEntry(h,Form("#color[%d]{%s}",color[j],labels[j].Data()),"l");
    if(j==0){
      h->SetAxisRange(-40,20);
      c1->cd(1)->SetLogy();
      h->Draw();
      c1->cd(2);
      h->Draw();
    }else{
      c1->cd(1)->SetLogy();
      h->Draw("same");
      c1->cd(2);
      h->Draw("same");
    }
  }
  leg->Draw();
  return c1;
}
//____________________________________________________________________________
TCanvas*
TOF_Btrg( int type, std::string strTrig )
{
  TString name=Form("TOF_Btrg_%s",strTrig.data());
  TCanvas *c1 = new TCanvas(name,name);
  TH1 *h;
  TString labels[4]={"Btrg All","Pion","Kaon","Proton"};
  int color[4]={1,2,3,4};
  int width[4]={1,1,2,2};
  int ch[4]={7,11,12,13};
  c1->Divide(1,2);
  TLegend *leg=new TLegend(0.1,0.6,0.3,0.9);
  for( int j=0;j<4;j++){
    h = GHist::get( HistMaker::getUniqueID(kAna, 0, 1, 20*type+ch[j]) ); // all
    h->SetLineColor(color[j]); 
    h->SetLineWidth(width[j]); 
    leg->AddEntry(h,Form("#color[%d]{%s}",j+1,labels[j].Data()),"l");
    if(j==0){
      h->SetAxisRange(-40,20);
      c1->cd(1)->SetLogy();
      h->Draw();
      c1->cd(2);
      h->Draw();
    }else{
      c1->cd(1)->SetLogy();
      h->Draw("same");
      c1->cd(2);
      h->Draw("same");
    }
  }
  leg->Draw();
  return c1;
}
//____________________________________________________________________________
TCanvas*
TOF2( int type, std::string strTrig )
{
  TString name=Form("TOF2_%s",strTrig.data());
  TCanvas *c1 = new TCanvas(name,name);
  TH1 *h;
  TString labels[4]={"Beam trigger","Pion trigger","Kaon trigger","Proton trigger"};
  c1->Divide(2,2);
  for( int j=0;j<4;j++){
    h = GHist::get( HistMaker::getUniqueID(kAna, 0, 1, 20*type+j+7) ); // all
    c1->cd(j+1)->SetGrid();
    h->Draw();
    h->SetAxisRange(-40,20);
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
TOF2D( int type, std::string strTrig )
{
  TString name=Form("TOF2D_%s",strTrig.data());
  TCanvas *c1 = new TCanvas(name,name);
  TH1 *h;
  int ch[4]={1,3,4,5};
  c1->Divide(2,2);
  for( int j=0;j<4;j++){
    h = GHist::get( HistMaker::getUniqueID(kAna, 0, type, ch[j]) ); // all
    c1->cd(j+1)->SetLogz();
    gPad->SetGrid();
    h->Draw("colz");
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
AC(int type, std::string strTrig)
{
  TString name=Form("AC_%s",strTrig.data());
  TCanvas *c1 = new TCanvas(name,name);
  TString labels[4]={"All","ACHit","TOF Kaon","TOF Pion"};
  TH1 *h;
  c1->SetLogy();
  int color[4]={1,2,4,3};
  int width[4]={1,1,2,2};
  TLegend *leg=new TLegend(0.6,0.6,0.9,0.9);
  for(int i=0;i<4;i++){
    h = GHist::get( HistMaker::getUniqueID(kAna, 0, type, i+1) );
    if( !h ) continue;
    h->SetLineColor(color[i]);
    h->SetLineWidth(width[i]);
    leg->AddEntry(h,Form("#color[%d]{%s}",color[i],labels[i].Data()),"l");
    if(i==0) h->Draw();
    else{
      h->Draw("same");
    }
    leg->Draw();
  }
  return c1;
}
//____________________________________________________________________________
TCanvas*
MHTDCHitPatMulti( DetectorType det, std::string strDet )
{
  TString cname=__func__+strDet;
  TCanvas *c1 = new TCanvas(cname, cname);
  c1->Divide(2,2);
  {
    c1->cd(1);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kHitPat, 1) );
    if( h ){
      h->Draw("hist");
      //    h->GetXaxis()->SetRangeUser(0,2000);
    }
    h = GHist::get( HistMaker::getUniqueID(det, 0, kHitPat, 2) );
    if( h ){
      h->SetLineColor(2);
      h->Draw("histsame");
      //    h->GetXaxis()->SetRangeUser(0,2000);
    }
  }
  {
    c1->cd(2);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kMulti, 1) );
    if( h ){
      h->Draw("hist");
      //    h->GetXaxis()->SetRangeUser(0,2000);
    }
    h = GHist::get( HistMaker::getUniqueID(det, 0, kMulti, 2) );
    if( h ){
      h->SetLineColor(2);
      h->Draw("histsame");
      //    h->GetXaxis()->SetRangeUser(0,2000);
    }
  }
  if(det!=kTriggerFlag){
    c1->cd(3);
    TH1 *h = GHist::get( HistMaker::getUniqueID(det, 0, kHitPat2D, 0) );
    if( h ){
      h->Draw("box");
      //    h->GetXaxis()->SetRangeUser(0,2000);
    }
  }
  return c1;
}
//____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
void
UpdateBLDCEfficiency()
{
  const int nchm=4;
  static TCanvas *c1[nchm];
  TString chm_name[nchm]={"BLC1a","BLC1b","BLC2a","BLC2b"};
  int layer[nchm]={8,8,8,8};
  int tmpn=0;
  for(int i=0;i<nchm;i++) tmpn+=layer[i];
  static std::vector<TText*> tex(tmpn);
  tmpn=0;
  for(int ichm=0;ichm<nchm;ichm++){
    c1[ichm] = (TCanvas*)gROOT->FindObject("BLDCMulti_"+chm_name[ichm]);
    if(!c1[ichm]){
      std::cout<<"no canvas for "<<chm_name[ichm]<<std::endl;
      return;
    }
    std::vector<TString> name;
    TString dname=chm_name[ichm];
    for(int i=0;i<layer[ichm];i++){
      name.push_back(dname+Form("_Multi_%d",i));
    }
    for(Int_t i=0, n=name.size(); i<n; ++i){
      c1[ichm]->cd(i+1);
      TH1 *h = (TH1*)gPad->FindObject(name[i]);
      if(!h) continue;
      Double_t zero = h->GetBinContent(1);
      Double_t all  = h->GetEntries();
      Double_t eff  = 1. - zero/all;
      int itex=i+tmpn;
      if(tex[itex]) delete tex[itex];
      tex[itex] = new TText;
      tex[itex]->SetNDC();
      tex[itex]->SetTextSize(0.100);
      tex[itex]->SetText(0.500,0.600,Form("eff. %.3f",eff));
      tex[itex]->Draw();
    }
    tmpn+=layer[ichm];
  } 
}
//_____________________________________________________________________________
void
UpdateCounterEfficiency()
{
  std::vector<TString> ch_name = {
    "BHT_Multi_0", "T0_Multi_0", "BH2_Multi_0",
    "BAC_Multi_0", "HTOF_Multi_0", "KVC_Multi_0",
    "KVC_Multi_1", "CVC_Multi_0"
  };
  static std::vector<TText*> tex(ch_name.size());
  static auto c1 = (TCanvas*)gROOT->FindObject("Multiplicity");
  if(!c1) return;
  for(int i=0, n=ch_name.size(); i<n; ++i){
    c1->cd(i+1);
    TH1 *h = (TH1*)gPad->FindObject(ch_name[i]);
    if(!h) continue;
    Double_t zero = h->GetBinContent(1);
    Double_t all  = h->GetEntries();
    Double_t eff  = 1. - zero/all;
    if(tex[i]) delete tex[i];
    tex[i] = new TText;
    tex[i]->SetNDC();
    tex[i]->SetTextSize(0.080);
    tex[i]->SetText(0.490,0.700,Form("eff. %.3f",eff));
    tex[i]->Draw();
  }
} 
}
}

