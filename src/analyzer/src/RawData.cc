/**
 *  file: RawData.cc
 *  date: 2017.04.10
 *
 */

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "DCRawHit.hh"
#include "HodoRawHit.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

#define OscillationCut 0
#define MsTPhase1      1

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("RawData");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const UserParamMan&    gUser     = UserParamMan::GetInstance();
  enum EUorD { kOneSide=1, kBothSide=2 };
  enum EDCDataType { kLeading, kTrailing, kNDCDataType };
#if OscillationCut
  const int  MaxMultiHitDC  = 16;
#endif
  
  //______________________________________________________________________________
  inline bool
  AddHodoRawHit( HodoRHitContainer& cont,
		 int id, int plane, int seg, int UorD, int AorT, int data, bool DEBUG=false )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    
    HodoRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      HodoRawHit *q = cont[i];
      if( q->DetectorId()==id &&
	  q->PlaneId()==plane &&
	  q->SegmentId()==seg ){
	if(DEBUG) std::cout<<"add data"<<std::endl;
	p=q; break;
      }
    }
    if( !p ){
      p = new HodoRawHit( id, plane, seg );
      cont.push_back(p);
      if(DEBUG) std::cout<<"========= new hodorawhit ======="<<std::endl;
    }
    if(DEBUG){
    //    if(id==DetIdBHT){
      std::cout   << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << AorT  << std::endl
		  << "UorD       = " << UorD  << std::endl
		  << "Data       = " << data  << std::endl;      
    }
    switch(AorT){
    case 0:
      if( UorD==0 ) p->SetAdcUp(data);
      else          p->SetAdcDown(data);
      break;
    case 1:
      if( UorD==0 ) p->SetTdcUp(data);
      else          p->SetTdcDown(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong AorT " << std::endl
		  << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << AorT  << std::endl
		  << "UorD       = " << UorD  << std::endl;
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  inline bool
  AddDCRawHit( DCRHitContainer& cont,
	       int plane, int wire, int tdc, int type=kLeading )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");
    
    DCRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      DCRawHit *q = cont[i];
      if( q->PlaneId()==plane &&
	  q->WireId()==wire ){
	p=q; break;
      }
    }
    if( !p ){
      p = new DCRawHit( plane, wire );
      cont.push_back(p);
    }

    switch(type){
    case kLeading:
      p->SetTdc(tdc);
      break;
    case kTrailing:
      p->SetTrailing(tdc);
      break;
    default:
      hddaq::cerr << func_name << " wrong data type " << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "WireId     = " << wire  << std::endl
		  << "DataType   = " << type  << std::endl;
      return false;
    }
    return true;
  }
  
  //______________________________________________________________________________
  inline void
  DecodeHodo( int id, int plane, int nseg, int nch, HodoRHitContainer& cont, bool DEBUG=false )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	for( int AorT=0; AorT<2; ++AorT ){
	  int nhit = gUnpacker.get_entries( id, plane, seg, UorD, AorT );
	  if( nhit<=0 ) continue;
	  for(int ihit=0;ihit<nhit;++ihit){
	    int data = gUnpacker.get( id, plane, seg, UorD, AorT, ihit );
	    AddHodoRawHit( cont, id, plane, seg, UorD, AorT, data, DEBUG );
	  }
	}
      }
    }
  }
  inline void
  DecodeBHT( int id, int plane, int nseg, int nch, HodoRHitContainer& cont, bool DEBUG=false )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	int nhit = gUnpacker.get_entries( id, plane, seg, UorD, kLeading ); //leading
	if( nhit<=0 ) continue;
	for(int ihit=0;ihit<nhit;++ihit){
	  int data = gUnpacker.get( id, plane, seg, UorD, kLeading, ihit );
	  AddHodoRawHit( cont, id, plane, seg, UorD, 1, data, DEBUG );
	}
	nhit = gUnpacker.get_entries( id, plane, seg, UorD, kTrailing ); //trailing
	if( nhit<=0 ) continue;
	for(int ihit=0;ihit<nhit;++ihit){
	  int data = gUnpacker.get( id, plane, seg, UorD, kTrailing, ihit );
	  AddHodoRawHit( cont, id, plane, seg, UorD, 0, data, DEBUG );
	}
      }
    }
  }
  inline void
  DecodeDC( int id, int layer, int nwire, DCRHitContainer& cont, int tdcmin, int tdcmax)
  {
    for(int wire=0; wire<nwire; ++wire){
      int nhit = gUnpacker.get_entries( id, layer, 0, wire, 0 );
#if OscillationCut
      if( nhit>MaxMultiHitDC ) continue;
#endif
      for(int i=0; i<nhit; i++ ){
	int data = gUnpacker.get( id, layer, 0, wire, 0, i);
	if( data<tdcmin || tdcmax<data ) continue;
	AddDCRawHit( cont, layer, wire, data );
      }
    }
  }
  
  //______________________________________________________________________________
  inline void
  DecodeHodo( int id, int nseg, int nch, HodoRHitContainer& cont )
  {
    DecodeHodo( id, 0, nseg, nch, cont );
  }
  
}

//______________________________________________________________________________
RawData::RawData( void )
  : m_is_decoded(false),
    m_BHDRawHC(),
    m_T0RawHC(),
    m_CVCRawHC(),
    m_NCRawHC(),
    m_BLC1aRawHC(8+1),
    m_BLC1bRawHC(8+1),
    m_BLC2aRawHC(8+1),
    m_BLC2bRawHC(8+1),
    m_BcOutRawHC(NumOfLayersBcOut+1)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
RawData::~RawData( void )
{
  ClearAll();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
RawData::ClearAll( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  del::ClearContainer( m_BHDRawHC );
  del::ClearContainer( m_T0RawHC );
  // del::ClearContainer( m_CVCRawHC );
  // del::ClearContainer( m_NCRawHC );
  del::ClearContainerAll( m_BLC1aRawHC );
  del::ClearContainerAll( m_BLC1bRawHC );
  del::ClearContainerAll( m_BLC2aRawHC );
  del::ClearContainerAll( m_BLC2bRawHC );
  del::ClearContainerAll( m_BcOutRawHC );
}

//______________________________________________________________________________
bool
RawData::DecodeHits( void )
{

  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_decoded ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();
  // DecodeHodo( DetIdCVC, 0, 10, kBothSide, m_CVCRawHC );
  // DecodeHodo( DetIdNC,  0,  6, kBothSide, m_NCRawHC );
  DecodeHodo( DetIdBHD, 0, 16, kBothSide, m_BHDRawHC );
  DecodeHodo( DetIdT0,  0, 5, kBothSide, m_T0RawHC );
  for(int i=0;i<8;i++){
    DecodeDC( DetIdBLC1a, i, 32, m_BLC1aRawHC[i], 0,2000 );
    DecodeDC( DetIdBLC1b, i, 32, m_BLC1bRawHC[i], 0,2000 );
    DecodeDC( DetIdBLC2a, i, 32, m_BLC2aRawHC[i], 0,2000 );
    DecodeDC( DetIdBLC2b, i, 32, m_BLC2bRawHC[i], 0,2000 );

    //BcOut
    DecodeDC( DetIdBLC2a, i, 32, m_BcOutRawHC[i],0,2000);
    DecodeDC( DetIdBLC2b, i, 32, m_BcOutRawHC[i+8],0,2000);
    
  }

  m_is_decoded = true;
  return true;
}

//______________________________________________________________________________
bool
RawData::DecodeCalibHits( void )
{
  // del::ClearContainer( m_VmeCalibRawHC );

  // for( int plane=0; plane<NumOfPlaneVmeCalib; ++plane ){
  //   DecodeHodo( DetIdVmeCalib, plane, NumOfSegVmeCalib,
  // 		kOneSide, m_VmeCalibRawHC );
  // }

  return true;
}
