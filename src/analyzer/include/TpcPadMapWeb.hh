// -*- C++ -*-

#ifndef TPC_PAD_MAP_WEB_HH
#define TPC_PAD_MAP_WEB_HH

class TH1;


namespace tpcweb
{
void Initialize(TH1* adc2d, TH1* tdc2d, TH1* hit2d, TH1* thre2d);
void Update();
}

#endif
