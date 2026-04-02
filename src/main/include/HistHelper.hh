// Author: Tomonori Takahashi

#ifndef ANALYZER_HISTOGRAM_HELPER_H
#define ANALYZER_HISTOGRAM_HELPER_H

#include <string>
#include <vector>
#include <map>
#include <Rtypes.h>
#include <TString.h>

class TH1;
class TH2;
class TH3;

class GHist
{

protected:
  std::map<Int_t, TString> m_index;

public:
  static GHist& getInstance( void );
  virtual ~GHist( void );

  // methods for TH1 --------------------------------------------------
  static TH1* D1(const TString& name,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xup);
  static TH1* D1(Int_t id,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xup);
  static TH1* I1(const TString& name,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xup);
  static TH1* I1(Int_t id,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xup);
  static Int_t  fill1(const TString& name,
		    Double_t x,
		    Double_t w=1.0);
  static Int_t  fill1(Int_t id,
		    Double_t x,
		    Double_t w=1.0);
  static Int_t  fill1(TH1* h,
		    Double_t x,
		    Double_t w=1.0);

  // methods for TH2 --------------------------------------------------
  static TH2* D2(const TString& name,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xip,
		 Int_t nbinsy,
		 Double_t ylow,
		 Double_t yup);
  static TH2* D2(Int_t id,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xip,
		 Int_t nbinsy,
		 Double_t ylow,
		 Double_t yup);
  static TH2* I2(const TString& name,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xip,
		 Int_t nbinsy,
		 Double_t ylow,
		 Double_t yup);
  static TH2* I2(Int_t id,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xip,
		 Int_t nbinsy,
		 Double_t ylow,
		 Double_t yup);
  static TH2* P2(const TString& name,
		 const TString& title,
		 Double_t xlow,
		 Double_t xip,
		 Double_t ylow,
		 Double_t yup);
  static TH2* P2(Int_t id,
		 const TString& title,
		 Double_t xlow,
		 Double_t xip,
		 Double_t ylow,
		 Double_t yup);
  static TH3* D3(Int_t id,
		 const TString& title,
		 Int_t nbinsx,
		 Double_t xlow,
		 Double_t xip,
		 Int_t nbinsy,
		 Double_t ylow,
		 Double_t yup,
		 Int_t nbinsz,
		 Double_t zlow,
		 Double_t zup);
  static Int_t  fill2(const TString& name,
		    Double_t x,
		    Double_t y,
		    Double_t w=1.0);
  static Int_t  fill2(Int_t id,
		    Double_t x,
		    Double_t y,
		    Double_t w=1.0);
  static Int_t  fill2(TH2* h,
		    Double_t x,
		    Double_t y,
		    Double_t w=1.0);

  // Common mothod ------------------------------------------------------
  // No function which return TH2 pointer is prepared.
  // If the actual class of the pointer is TH2,
  // dynamic_cast in the user source code is need.
  static TH1* get(const TString& name);
  static TH1* get(Int_t id);

private:
  GHist();
  GHist(const GHist&);
  GHist& operator=(const GHist&);

  ClassDef(GHist, 0)

};

#endif
