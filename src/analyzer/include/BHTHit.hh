// -*- C++ -*-

#ifndef BHT_HIT_HH
#define BHT_HIT_HH

#include <string>

#include <TObject.h>
#include <TString.h>

#include <std_ostream.hh>

#include "Hodo2Hit.hh"

//______________________________________________________________________________
class BHTHit : public Hodo2Hit
{
public:
  explicit  BHTHit( HodoRawHit *object, int index=0 );
  virtual  ~BHTHit( void );

private:
  BHTHit( const BHTHit& object );
  BHTHit& operator =( const BHTHit& object );

public:
  Bool_t   Calculate( void ) override;
  Double_t GetWidthTop( Int_t n=0 ) const { return Double_t(GetAUp(n)); }
  Double_t GetWidthBottom( Int_t n=0 ) const { return Double_t(GetADown(n)); }
  Double_t GetTotTop( Int_t n=0 )      const { return GetWidthTop(n);      }
  Double_t GetTotBottom( Int_t n=0 )   const { return GetWidthBottom(n);   }

  virtual Bool_t ReCalc( Bool_t allpyRecursively=false )
  { return BHTHit::Calculate(); }

  //  ClassDef(BHTHit,0);
};

#endif
