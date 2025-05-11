// -*- C++ -*-
#include "BHTHit.hh"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#define PHC 0
//ClassImp(BHTHit);

namespace
{
  const std::string& class_name("BHTHit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
}

//______________________________________________________________________________
BHTHit::BHTHit( HodoRawHit *object, int index)
  : Hodo2Hit(object,index)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
BHTHit::~BHTHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
Bool_t
BHTHit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
//  if( m_status ){
  if( m_is_calculated ){
    hddaq::cerr << func_name << " already calculated" << std::endl;
    return false;
  }

  if( m_raw->GetNumOfTdcHits()!=2 )
    return false;
  if( !gHodo.IsReady() ){
    hddaq::cerr << func_name << " HodoParamMan must be initialized" << std::endl;
    return false;
  }
#if PHC
  if( !gPHC.IsReady() ){
    hddaq::cerr << func_name << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }
#endif
  // Detector information
  Int_t cid  = m_raw->DetectorId();
  Int_t plid = m_raw->PlaneId();
  Int_t seg  = m_raw->SegmentId();


  for(int i=0;i<m_raw->GetSizeTdcUp();i++){ 
    int tdc1 = m_raw->GetTdc1(i);
    if(m_raw->GetSizeTdcDown()<i+1) break;
    int tdc2 = m_raw->GetTdc2(i);
    if( tdc1<0 || tdc2<0 ){
      if(i==0) return false;
      else break;
    }   
    if(m_raw->GetSizeAdcUp()<i+1) break;
    if(m_raw->GetSizeAdcDown()<i+1) break;
    double time1 = 0., time2 = 0.;
    if( !gHodo.GetTime( cid, plid, seg, 0, tdc1, time1 ) ||
	!gHodo.GetTime( cid, plid, seg, 1, tdc2, time2 ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "U/D, " << tdc1 << "/" << tdc2 << ", " << "time"
		  << ")" << std::endl;
      return false;
    }
    m_t1.push_back( time1 );
    m_t2.push_back( time2 );
    m_a1.push_back( tdc1-m_raw->GetAdcUp(i) );
    m_a2.push_back( tdc2-m_raw->GetAdcDown(i) );
    m_index++;
  }
  m_is_calculated = true;
  return true;
}
