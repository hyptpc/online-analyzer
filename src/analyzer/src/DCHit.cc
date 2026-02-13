/**
 *  file: DCHit.cc
 *  date: 2017.04.10
 *
 */

#include "DCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "BLDCWireMapMan.hh"
#include "XTMapMan.hh"
//#include "DCDriftParamMan.hh"
#include "DetectorID.hh"
//#include "DCGeomMan.hh"
#include "DCParameters.hh"
#include "DCTdcCalibMan.hh"
#include "DCLTrackHit.hh"
#include "DebugCounter.hh"
#include "MathTools.hh"
//#include "RootHelper.hh"

#include "TRandom.h"

namespace
{
  const std::string& class_name("DCHit");
  //  const DCGeomMan&       gGeom  = DCGeomMan::GetInstance();
  const BLDCWireMapMan&       gGeom  = BLDCWireMapMan::GetInstance();
  const DCTdcCalibMan&   gTdc   = DCTdcCalibMan::GetInstance();
  const XTMapMan& gXt = XTMapMan::GetInstance();
  const bool SelectTDC1st  = true;
}

//______________________________________________________________________________
DCHit::DCHit( void )
  : m_layer(-1), m_wire(-1),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(9999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer )
  : m_layer( layer ), m_wire(-1),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(-999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer, double wire )
  : m_layer(layer), m_wire(wire),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(-999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int cid,int layer, double wire )
  : m_cid(cid), m_layer(layer), m_wire(wire),
    m_wpos(-9999.), m_angle(0.),
    m_cluster_size(0.),
    m_zero_suppressed(false),
    m_time_corrected(false),
    m_good_waveform(false),
    m_pedestal(-999),
    m_peak_height(-999),
    m_peak_position(-999),
    m_deviation(-999.),
    m_amplitude(-999.),
    m_peak_time(-999.),
    m_adc_sum(-999.),
    m_de(-999.),
    m_rms(-999.),
    m_chisqr(-999.),
    m_belong_kaon(false)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::~DCHit( void )
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
DCHit::SetTdcVal( int tdc )
{
  m_tdc.push_back(tdc);
  m_belong_track.push_back(false);
  m_dl_range.push_back(false);
}

//______________________________________________________________________________
void
DCHit::SetAdcVal(int adc )
{
  m_adc.push_back(adc);
}

//______________________________________________________________________________
void
DCHit::ClearRegisteredHits( void )
{
  int n = m_register_container.size();
  for(int i=0; i<n; ++i){
    delete m_register_container[i];
  }
}

//______________________________________________________________________________
bool
DCHit::CalcDCObservables( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gGeom.IsReady() ||
      !gTdc.IsReady() ||
      !gXt.IsReady() )
    return false;

  m_wpos  = gGeom.CalcWirePosition( m_cid, m_layer,m_wire ) *10;
  m_angle = gGeom.GetTiltAngle( m_cid, m_layer );
  m_z     = gGeom.GetLocalZ( m_cid, m_layer )*10;
  //  std::cout<<func_name<<"  "<<m_cid<<"  "<<m_layer<<"  "<<m_wire<<"  "<<m_wpos<<std::endl;
  bool status = true;
  int  nhtdc = m_tdc.size();
  for ( int i=0; i<nhtdc; ++i ) {
    Double_t ctime;
    if( !gTdc.GetTime( m_cid, m_layer, m_wire, m_tdc[i], ctime ) )
      return false;
    double dtime=ctime;
    double dlength=gXt.CalcDriftLength( m_cid, m_layer, m_wire, dtime );
    //    std::cout<<i<<" / "<<nhtdc<<"  tdc:"<<m_tdc[i]<<"  dt:"<<ctime<<" ,dl:"<<dlength<<std::endl;
    m_dt.push_back( dtime );
    m_dl.push_back( dlength );
    m_dl_range[i] = false;
    switch(m_cid){
    case DetIdBLC1a:
      if( 0.01<m_dl[i] && m_dl[i]<3.99 )
  	m_dl_range[i] = true;
      break;      
    case DetIdBLC1b:
      if( 0.01<m_dl[i] && m_dl[i]<3.99 )
  	m_dl_range[i] = true;
      break;      
    case DetIdBLC2a:
      if( 0.01<m_dl[i] && m_dl[i]<2.49 )
  	m_dl_range[i] = true;
      break;      
    case DetIdBLC2b:
      if( 0.01<m_dl[i] && m_dl[i]<2.49 )
  	m_dl_range[i] = true;
      break;      
    case DetIdBPC1:
      if( 0.01<m_dl[i] && m_dl[i]<3.59 )
  	m_dl_range[i] = true;
      break;      
    case DetIdBPC2:
      if( 0.01<m_dl[i] && m_dl[i]<2.99 )
  	m_dl_range[i] = true;
      break;      
    default:
      hddaq::cout << "#E " << func_name << " "
  		  << "invalid counter id : " << m_cid << std::endl;
      status = false;
      break;
    }
  }
  
  int    tdc1st = 0;
  double dl1st  = -1.;
  double dt1st  = -1.;
  for( int i=0; i<nhtdc; ++i ){
    if( tdc1st < m_tdc[i]
	&& m_dl_range[i] ){
      tdc1st = m_tdc[i];
      dl1st  = m_dl[i];
      dt1st  = m_dt[i];
    }
  }
  if( tdc1st>0 ){
    m_tdc.clear();
    m_dt.clear();
    m_dl.clear();
    m_tdc.push_back( tdc1st );
    m_dt.push_back( dt1st );
    m_dl.push_back( dl1st );
  } else {
    m_tdc.clear();
    m_dt.clear();
    m_dl.clear();
  }
  //  std::cout<<"Fill hits"<<std::endl;
  // std::cout<<func_name<<"  "<<m_cid<<"  "<<m_layer+1<<"  "<<m_wire+1<<"  wpos:"<<m_wpos<<"  ,angle: "<<m_angle<<"  ,zpos: "<<m_z<<std::endl;
  // std::cout<<"STATUS: "<<status<<std::endl;
  
  return status;
}

//______________________________________________________________________________
bool
DCHit::DoTimeCorrection( double offset )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_time_corrected ){
    Print(func_name+" already corrected!");
    return false;
  }

#if 0
  Print("Before Correction");
#endif

  DoubleVec ctime;
  int nh = m_time.size();
  for(int i=0; i<nh; ++i){
    ctime.push_back( m_time[i] + offset );
  }

  m_time = ctime;
  m_peak_time    += offset;
  m_time_corrected = true;

#if 0
  Print("After Correction");
#endif

  return true;
}

//______________________________________________________________________________
double
DCHit::GetResolution( void ) const
{
  return 0.1;
  //  return gGeom.GetResolution(m_layer);
}

//______________________________________________________________________________
void
DCHit::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int w = 16;
  ost << "#D " << func_name << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
      << std::setw(w) << std::left << "wpos"  << m_wpos  << std::endl
      << std::setw(w) << std::left << "angle" << m_angle << std::endl
      << std::setw(w) << std::left << "z"     << m_z     << std::endl
      << std::setw(w) << std::left << "kaon"  << m_belong_kaon << std::endl;

  ost << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy( m_tdc.begin(), m_tdc.end(),
	     std::ostream_iterator<int>(ost, " ") );
  ost << std::endl;

  if(!m_is_ssd){
    ost << std::endl << std::setw(w) << std::left
	<< "trailing" << m_trailing.size() << " : ";
    std::copy(m_trailing.begin(), m_trailing.end(),
	      std::ostream_iterator<int>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "drift time" << m_dt.size() << " : ";
    std::copy(m_dt.begin(), m_dt.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "drift length" << m_dl.size() << " : ";
    std::copy(m_dl.begin(), m_dl.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "trailing time" << m_trailing_time.size() << " : ";
    std::copy(m_trailing_time.begin(), m_trailing_time.end(),
	      std::ostream_iterator<double>(ost, " "));
    ost << std::endl;
    ost << std::setw(w) << std::left
	<< "belongTrack" << m_belong_track.size() << " : ";
    std::copy(m_belong_track.begin(), m_belong_track.end(),
	      std::ostream_iterator<bool>(ost, " "));
    ost << std::endl << std::setw(w) << std::left
	<< "dlRange" << m_dl_range.size() << " : ";
    std::copy(m_dl_range.begin(), m_dl_range.end(),
	      std::ostream_iterator<bool>(ost, " "));
    ost << std::endl;
  }
}
