// Author: Tomonori Takahashi

#include "HistHelper.hh"

#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH2Poly.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TString.h>

#include "lexical_cast.hh"

ClassImp(GHist)

//______________________________________________________________________________
GHist::GHist( void )
  : m_index()
{
}

//______________________________________________________________________________
GHist::~GHist( void )
{
}
//______________________________________________________________________________
TH1*
GHist::D1( const TString& name,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup )
{
  return new TH1D( name, title, nbinsx, xlow, xup );
}

//______________________________________________________________________________
TH1*
GHist::D1( Int_t id,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if( !g_h.m_index[id].IsNull() ){
    std::cerr << "#E HistHelper::D1" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH1D( name, title, nbinsx, xlow, xup );
}

//______________________________________________________________________________
TH1*
GHist::I1( const TString& name,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup )
{
  return new TH1I( name, title, nbinsx, xlow, xup );
}

//______________________________________________________________________________
TH1*
GHist::I1( Int_t id,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if("" != g_h.m_index[id]){
    std::cerr << "#E HistHelper::I1" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH1I( name, title, nbinsx, xlow, xup );
}

//______________________________________________________________________________
Int_t
GHist::fill1( const TString& name,
	      Double_t x,
	      Double_t w )
{
  TH1* h = GHist::get(name);
  if (h)
    return h->Fill(x, w);
  return 0;
}

//______________________________________________________________________________
Int_t
GHist::fill1( Int_t id,
	      Double_t x,
	      Double_t w )
{
  TH1* h = GHist::get(id);
  if (h)
    return h->Fill(x, w);
  return 0;
}

//______________________________________________________________________________
Int_t
GHist::fill1( TH1* h,
	      Double_t x,
	      Double_t w )
{

  if (h)
    {
//       std::cout << "#D H1::fill() \n"
// 		<< h->GetName() << " x = " << x << " w = " << w
// 		<< std::endl;
      return h->Fill(x, w);
    }
  else
    std::cerr << "#E GHist::fill()\n"
	      << " got null pointer" << std::endl;
  return 0;
}

//______________________________________________________________________________
TH2*
GHist::D2( const TString& name,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup,
	   Int_t nbinsy,
	   Double_t ylow,
	   Double_t yup )
{
  return new TH2D( name, title,
		   nbinsx, xlow, xup,
		   nbinsy, ylow, yup );
}

//______________________________________________________________________________
TH2*
GHist::D2( Int_t id,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup,
	   Int_t nbinsy,
	   Double_t ylow,
	   Double_t yup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if( !g_h.m_index[id].IsNull() ){
    std::cerr << "#E HistHelper::I1" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH2D( name, title,
		   nbinsx, xlow, xup,
		   nbinsy, ylow, yup );
}

//______________________________________________________________________________
TH2*
GHist::I2( const TString& name,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup,
	   Int_t nbinsy,
	   Double_t ylow,
	   Double_t yup )
{
  return new TH2I( name, title,
		   nbinsx, xlow, xup,
		   nbinsy, ylow, yup );
}

//______________________________________________________________________________
TH2*
GHist::I2( Int_t id,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup,
	   Int_t nbinsy,
	   Double_t ylow,
	   Double_t yup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if( !g_h.m_index[id].IsNull() ){
    std::cerr << "#E HistHelper::I2" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH2I( name, title,
		   nbinsx, xlow, xup,
		   nbinsy, ylow, yup );
}

//______________________________________________________________________________
TH2*
GHist::P2( const TString& name,
	   const TString& title,
	   Double_t xlow,
	   Double_t xup,
	   Double_t ylow,
	   Double_t yup )
{
  return new TH2Poly( name, title,
                      xlow, xup,
                      ylow, yup );
}

//______________________________________________________________________________
TH2*
GHist::P2( Int_t id,
	   const TString& title,
	   Double_t xlow,
	   Double_t xup,
	   Double_t ylow,
	   Double_t yup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if( !g_h.m_index[id].IsNull() ){
    std::cerr << "#E HistHelper::P2" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH2Poly( name, title,
                      xlow, xup,
                      ylow, yup );
}

//______________________________________________________________________________
TH3*
GHist::D3( Int_t id,
	   const TString& title,
	   Int_t nbinsx,
	   Double_t xlow,
	   Double_t xup,
	   Int_t nbinsy,
	   Double_t ylow,
	   Double_t yup,
	   Int_t nbinsz,
	   Double_t zlow,
	   Double_t zup )
{
  TString name = title;
  GHist& g_h = GHist::getInstance();
  if( !g_h.m_index[id].IsNull() ){
    std::cerr << "#E HistHelper::D3" << std::endl;
    std::cerr << " The histogram already exists with the same ID" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  g_h.m_index[id] = name;
  return new TH3D( name, title,
		   nbinsx, xlow, xup,
		   nbinsy, ylow, yup,
		   nbinsz, zlow, zup);
}

//______________________________________________________________________________
Int_t
GHist::fill2( const TString& name,
	      Double_t x,
	      Double_t y,
	      Double_t w )
{
  TH2* h = dynamic_cast<TH2*>(GHist::get(name));
  if (h)
    return h->Fill(x, y, w);
  return 0;
}

//______________________________________________________________________________
Int_t
GHist::fill2( Int_t id,
	      Double_t x,
	      Double_t y,
	      Double_t w )
{
  TH2* h = dynamic_cast<TH2*>(GHist::get(id));
  if (h)
    return h->Fill(x, y, w);
  return 0;
}

//______________________________________________________________________________
Int_t
GHist::fill2( TH2* h,
	      Double_t x,
	      Double_t y,
	      Double_t w )
{
  if (h)
    return h->Fill(x, y, w);
  return 0;
}

//______________________________________________________________________________
//______________________________________________________________________________
TH1*
GHist::get( const TString& name )
{
  TObject *ptr = gROOT->FindObject(name);
  if( !ptr ){
    std::cerr << "#E HistHelper::get" << std::endl;
    std::cerr << " NULL pointer is returned" << std::endl;
    std::cerr << " name : " << name << std::endl;
    return NULL;
  }
  return dynamic_cast<TH1*>(ptr);
}

//______________________________________________________________________________
TH1*
GHist::get(Int_t id)
{
  TH1 *ptr = GHist::get(GHist::getInstance().m_index[id]);
  if( !ptr ){
    std::cerr << "#E HistHelper::get" << std::endl;
    std::cerr << " NULL pointer is returned" << std::endl;
    std::cerr << " ID : " << id << std::endl;
    return NULL;
  }
  return ptr;
}

//______________________________________________________________________________
GHist&
GHist::getInstance( void )
{
  static GHist g_instance;
  return g_instance;
}
