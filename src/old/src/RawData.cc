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

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("RawData");
  const UnpackerManager& gUnpacker = GUnpacker::get_instance();
  const UserParamMan&    gUser     = UserParamMan::GetInstance();
  enum EUorD { kOneSide=1, kBothSide=2 };
  enum EHodoDataType { kHodoAdc, kHodoLeading, kHodoTrailing, kHodoOverflow, kHodoNDataType };
  enum EDCDataType   { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
#if OscillationCut
  const int  MaxMultiHitDC  = 16;
#endif

  //______________________________________________________________________________
  inline bool
  AddHodoRawHit( HodoRHitContainer& cont,
		 int id, int plane, int seg, int UorD, int type, int data )
  {
    static const std::string func_name("["+class_name+"::"+__func__+"()]");

    HodoRawHit *p = 0;
    for( std::size_t i=0, n=cont.size(); i<n; ++i ){
      HodoRawHit *q = cont[i];
      if( q->DetectorId()==id &&
	  q->PlaneId()==plane &&
	  q->SegmentId()==seg ){
	p=q; break;
      }
    }
    if( !p ){
      p = new HodoRawHit( id, plane, seg );
      cont.push_back(p);
    }

    switch(type){
    case kHodoAdc:
      if( UorD==0 ) p->SetAdcUp(data);
      else          p->SetAdcDown(data);
      break;
    case kHodoLeading:
      if( UorD==0 ) p->SetTdcUp(data);
      else          p->SetTdcDown(data);
      break;
    case kHodoTrailing:
      if( UorD==0 ) p->SetTdcTUp(data);
      else          p->SetTdcTDown(data);
      break;
    case kHodoOverflow:
      p->SetTdcOverflow(data);
      break;
    default:
      hddaq::cerr << func_name << " wrong data type " << std::endl
		  << "DetectorId = " << id    << std::endl
		  << "PlaneId    = " << plane << std::endl
		  << "SegId      = " << seg   << std::endl
		  << "AorT       = " << type  << std::endl
		  << "UorD       = " << UorD  << std::endl;
      return false;
    }
    return true;
  }

  //______________________________________________________________________________
  inline bool
  AddDCRawHit( DCRHitContainer& cont,
	       int plane, int wire, int data, int type=kDcLeading )
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
    case kDcLeading:
      p->SetTdc(data);
      break;
    case kDcTrailing:
      p->SetTrailing(data);
      break;
    case kDcOverflow:
      p->SetTdcOverflow(data);
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
  DecodeHodo( int id, int plane, int nseg, int nch, HodoRHitContainer& cont )
  {
    for( int seg=0; seg<nseg; ++seg ){
      for( int UorD=0; UorD<nch; ++UorD ){
	for( int AorT=0; AorT<2; ++AorT ){
	  int nhit = gUnpacker.get_entries( id, plane, seg, UorD, AorT );
	  if( nhit<=0 ) continue;
	  for(int m = 0; m<nhit; ++m){
	    int data = gUnpacker.get( id, plane, seg, UorD, AorT, m );
	    AddHodoRawHit( cont, id, plane, seg, UorD, AorT, data );
	  }
	}
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
    m_BH1RawHC(),
    m_BH2RawHC(),
    m_TOFRawHC(),
    m_BFTRawHC(NumOfPlaneBFT),
    m_SCHRawHC(),
#ifdef E40
    m_SACRawHC(),
    m_HtTOFRawHC(),
    m_LCRawHC(),
    m_SFTRawHC(NumOfPlaneSFT),
    m_CFTRawHC(NumOfPlaneCFT),
    m_BGORawHC(),
    m_BGOFadcRawHC(NumOfSegBGO),
    m_PiIDRawHC(),
    m_FHT1RawHC(2*NumOfLayersFHT1),
    m_FHT2RawHC(2*NumOfLayersFHT2),
#endif
    m_BcInRawHC(NumOfLayersBcIn+1),
    m_BcOutRawHC(NumOfLayersBcOut+1),
    m_SdcInRawHC(NumOfLayersSdcIn+1),
    m_SdcOutRawHC(NumOfLayersSdcOut+1),
    m_ScalerRawHC(),
    m_TrigRawHC(),
    m_VmeCalibRawHC(),
    m_FpgaBH2MtRawHC()
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

  del::ClearContainer( m_BH1RawHC );
  del::ClearContainer( m_BH2RawHC );
  del::ClearContainer( m_TOFRawHC );

  del::ClearContainerAll( m_BFTRawHC );

#ifdef E40
  del::ClearContainer( m_SACRawHC );
  del::ClearContainer( m_HtTOFRawHC );
  del::ClearContainer( m_LCRawHC );
  del::ClearContainerAll( m_SFTRawHC );
  del::ClearContainerAll( m_CFTRawHC );
  del::ClearContainer( m_BGORawHC );
  //del::ClearContainerAll( m_BGOFadcRawHC );
  for (int i=0; i<NumOfSegBGO; i++) {
    m_BGOFadcRawHC[i].clear();
    m_BGOFadcRawHC[i].shrink_to_fit();
  }
  del::ClearContainer( m_PiIDRawHC );
  del::ClearContainerAll( m_FHT1RawHC );
  del::ClearContainerAll( m_FHT2RawHC );
#endif

  del::ClearContainer( m_SCHRawHC );

  del::ClearContainerAll( m_BcInRawHC );
  del::ClearContainerAll( m_BcOutRawHC );
  del::ClearContainerAll( m_SdcInRawHC );
  del::ClearContainerAll( m_SdcOutRawHC );

  del::ClearContainer( m_ScalerRawHC );
  del::ClearContainer( m_TrigRawHC );
  del::ClearContainer( m_VmeCalibRawHC );

  del::ClearContainer( m_FpgaBH2MtRawHC );
}

//______________________________________________________________________________
bool
RawData::DecodeHits( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  static const double MinBC3_TDC  = gUser.GetParameter("BC3_TDC", 0);
  static const double MaxBC3_TDC  = gUser.GetParameter("BC3_TDC", 1);
  static const double MinBC4_TDC  = gUser.GetParameter("BC4_TDC", 0);
  static const double MaxBC4_TDC  = gUser.GetParameter("BC4_TDC", 1);
  static const double MinSDC1_TDC = gUser.GetParameter("SDC1_TDC", 0);
  static const double MaxSDC1_TDC = gUser.GetParameter("SDC1_TDC", 1);
  // static const double MinSDC2_TDC = gUser.GetParameter("SDC2_TDC", 0);
  // static const double MaxSDC2_TDC = gUser.GetParameter("SDC2_TDC", 1);
  static const double MinSDC3_TDC = gUser.GetParameter("SDC3_TDC", 0);
  static const double MaxSDC3_TDC = gUser.GetParameter("SDC3_TDC", 1);
  static const double MinSDC4_TDC = gUser.GetParameter("SDC4_TDC", 0);
  static const double MaxSDC4_TDC = gUser.GetParameter("SDC4_TDC", 1);

  if( m_is_decoded ){
    hddaq::cout << "#D " << func_name << " "
		<< "already decoded!" << std::endl;
    return false;
  }

  ClearAll();

  // BH1
  DecodeHodo( DetIdBH1, NumOfSegBH1, kBothSide, m_BH1RawHC );
  // BH2
  DecodeHodo( DetIdBH2, NumOfSegBH2, kBothSide, m_BH2RawHC );
  // TOF
  DecodeHodo( DetIdTOF, NumOfSegTOF, kBothSide, m_TOFRawHC );

  //BFT
  for( int plane=0; plane<NumOfPlaneBFT; ++plane ){
    for(int seg = 0; seg<NumOfSegBFT; ++seg){
      int nhit = gUnpacker.get_entries( DetIdBFT, plane, 0, seg, 0 );
      if( nhit>0 ){
	for(int i = 0; i<nhit; ++i){
	  int leading  = gUnpacker.get( DetIdBFT, plane, 0, seg, 0, i )  ;
	  int trailing = gUnpacker.get( DetIdBFT, plane, 0, seg, 1, i )  ;
	  AddHodoRawHit( m_BFTRawHC[plane], DetIdBFT, plane, seg , 0, kHodoLeading,  leading );
	  AddHodoRawHit( m_BFTRawHC[plane], DetIdBFT, plane, seg , 0, kHodoTrailing, trailing );
	}
      }
    }
  }

  //SCH
  for(int seg=0; seg<NumOfSegSCH; ++seg){
    int nhit = gUnpacker.get_entries( DetIdSCH, 0, seg, 0, 0 );
    if( nhit>0 ){
      for(int i = 0; i<nhit; ++i){
	int leading  = gUnpacker.get( DetIdSCH, 0, seg, 0, 0, i );
	int trailing = gUnpacker.get( DetIdSCH, 0, seg, 0, 1, i );
	AddHodoRawHit( m_SCHRawHC, DetIdSCH, 0, seg , 0, kHodoLeading,  leading );
	AddHodoRawHit( m_SCHRawHC, DetIdSCH, 0, seg , 0, kHodoTrailing, trailing );
      }
    }
  }

#ifdef E40
  // SAC
  DecodeHodo( DetIdSAC, NumOfSegSAC, kOneSide,  m_SACRawHC );
  // TOF-HT
  DecodeHodo( DetIdHtTOF,NumOfSegHtTOF,kOneSide,  m_HtTOFRawHC );
  // LC
  DecodeHodo( DetIdLC,  NumOfSegLC,  kOneSide,  m_LCRawHC );
  //SFT
  for( int plane=0; plane<NumOfPlaneSFT; ++plane ){
    int nseg = 0;
    switch ( plane ) {
    case 0: nseg = NumOfSegSFT_UV;  break;
    case 1: nseg = NumOfSegSFT_UV;  break;
    case 2: nseg = NumOfSegSFT_X; break;
    case 3: nseg = NumOfSegSFT_X; break;
    default: break;
    }
    for ( int seg = 0; seg < nseg; ++seg ) {
      for ( int LorT = 0; LorT < 2; ++LorT ) {
	int nhit = gUnpacker.get_entries( DetIdSFT, plane, 0, seg, LorT );
	if( nhit>0 ){
	  for(int i = 0; i<nhit; ++i){
	    int edge = gUnpacker.get( DetIdSFT, plane, 0, seg, LorT, i )  ;
	    AddHodoRawHit( m_SFTRawHC[plane], DetIdSFT, plane, seg, 0, kHodoLeading+LorT, edge );
	  }
	}
      }
    }
  }

  //CFT
  for( int plane=0; plane<NumOfPlaneCFT; ++plane ){
    //CFT TDC
    for(int seg = 0; seg<NumOfSegCFT[plane]; ++seg){
      int nhit_tdc      = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 0 );
      int nhit_trailing = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 1 );
      if( nhit_tdc>0 ){
	for(int i = 0; i<nhit_tdc; ++i){
	  int leading  = gUnpacker.get( DetIdCFT, plane, seg, 0, 0, i );
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, kHodoLeading, leading );
	}
      }
      if(nhit_trailing>0){
	for(int i = 0; i<nhit_trailing; ++i){
	  int trailing = gUnpacker.get( DetIdCFT, plane, seg, 0, 1, i )  ;
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, kHodoTrailing, trailing );
	}
      }
      //if( nhit_tdc==0 )continue; // w/ or w/o TDC,   comment out => pedestal

      //CFT ADC HI
      int nhit_adc_hi = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 2 );
      if( nhit_adc_hi>0 ){
	for(int i = 0; i<nhit_adc_hi; ++i){
	  int adc_hi  = gUnpacker.get( DetIdCFT, plane, seg, 0, 2, i )  ;
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 0, 0, adc_hi );//ADC Hi
	}
      }
      else continue;
      //CFT ADC LOW
      int nhit_adc_low = gUnpacker.get_entries( DetIdCFT, plane, seg, 0, 3 );
      if( nhit_adc_low>0 ){
	for(int i = 0; i<nhit_adc_low; ++i){
	  int adc_low = gUnpacker.get( DetIdCFT, plane, seg, 0, 3, i );
	  AddHodoRawHit( m_CFTRawHC[plane], DetIdCFT, plane, seg , 1, 0, adc_low );//ADC Low
	}
      }
      else continue;
    }
  }

  //BGO
  for(int seg=0; seg<NumOfSegBGO; ++seg){
    int ped      = 0;
    int integral = 0;

    int nhit_a = gUnpacker.get_entries( DetIdBGO, 0, seg, 0, 0 );
    if( nhit_a>0 ){
      for(int i = 0; i<nhit_a; ++i){
	unsigned int fadc = gUnpacker.get(DetIdBGO, 0, seg, 0, 0 ,i);

	m_BGOFadcRawHC[seg].push_back(fadc);

	if( fadc == 0xffff ){
	}else{
	  if( ped == 0 ){ped = fadc;} // 1st value is seemed to be pedestal
	  int delta = ped - fadc;
	  /*
	  std::cout << "BGO : seg=" << seg  << ", i=" << i << "/nhit=" << nhit
		    << ", fadc=" << fadc << ", pede=" << ped
		    << ", integral=" << integral << ", +=" << delta // (ped - fadc)
		    << std::endl;
	  */
	  integral += delta;
	}
      }
      AddHodoRawHit( m_BGORawHC, DetIdBGO, 0, seg , 0, kHodoAdc, integral );
    }

    //TDC
    unsigned int nhit_t = gUnpacker.get_entries(DetIdBGO, 0, seg, 0, 1);
    for(unsigned int m = 0; m<nhit_t; ++m){
      int tdc = gUnpacker.get(DetIdBGO, 0, seg, 0, kHodoLeading, m);
      AddHodoRawHit( m_BGORawHC, DetIdBGO, 0, seg , 0, kHodoLeading, tdc );
      //std::cout << "BGO TDC : seg=" << seg  << ", tdc=" << tdc	<< std::endl;
    }

  }

  //PiID counter
  for(int seg=0; seg<NumOfSegPiID; ++seg){
    int nhit_l = gUnpacker.get_entries( DetIdPiID, 0, seg, 0, 0 );
    int nhit_t = gUnpacker.get_entries( DetIdPiID, 0, seg, 0, 1 );
    if( nhit_l>0 ){
      for(int i = 0; i<nhit_l; ++i){
	int leading  = gUnpacker.get( DetIdPiID, 0, seg, 0, 0, i );
	AddHodoRawHit( m_PiIDRawHC, DetIdPiID, 0, seg , 0, kHodoLeading,  leading );
	//std::cout << "PiID TDC : seg=" << seg  << ", tdc=" << leading << std::endl;
      }
    }
    if( nhit_t>0 ){
      for(int j = 0; j<nhit_t; ++j){
	int trailing = gUnpacker.get( DetIdPiID, 0, seg, 0, 1, j );
	AddHodoRawHit( m_PiIDRawHC, DetIdPiID, 0, seg , 0, kHodoTrailing, trailing );
	//std::cout << "PiID TDC trailing : seg=" << seg  << ", tdc=" << trailing << std::endl;
      }
    }

  }


  //FHT1
  for( int layer=0; layer<NumOfLayersFHT1; ++layer ){
    for(int seg = 0; seg<MaxSegFHT1; ++seg){
      for(int UorD = 0; UorD<2; ++UorD){
	for(int LorT = 0; LorT<2; ++LorT){
	  int nhit = gUnpacker.get_entries( DetIdFHT1, layer, seg, UorD, LorT);
	  for(int i = 0; i<nhit; ++i){
	    int time  = gUnpacker.get( DetIdFHT1, layer, seg, UorD, LorT, i )  ;
	    AddHodoRawHit( m_FHT1RawHC[2*layer + UorD], DetIdFHT1, layer, seg , UorD, kHodoLeading+LorT, time );
	  }// multihit
	}// LorT
      }// UorD
    }// seg
  }// layer

  //FHT2
  for( int layer=0; layer<NumOfLayersFHT2; ++layer ){
    for(int seg = 0; seg<MaxSegFHT2; ++seg){
      for(int UorD = 0; UorD<2; ++UorD){
	for(int LorT = 0; LorT<2; ++LorT){
	  int nhit = gUnpacker.get_entries( DetIdFHT2, layer, seg, UorD, LorT);
	  for(int i = 0; i<nhit; ++i){
	    int time  = gUnpacker.get( DetIdFHT2, layer, seg, UorD, LorT, i )  ;
	    AddHodoRawHit( m_FHT2RawHC[2*layer + UorD], DetIdFHT2, layer, seg , UorD, kHodoLeading+LorT, time );
	  }// multihit
	}// LorT
      }// UorD
    }// seg
  }// layer
#endif

  // BC3&BC4 MWDC
  for(int plane=0; plane<NumOfLayersBcOut; ++plane ){
    if( plane<NumOfLayersBc ){
      for(int wire=0; wire<MaxWireBC3; ++wire){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdBC3, plane, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdBC3, plane, 0, wire, lt, i);
	    if( data<MinBC3_TDC || MaxBC3_TDC<data ) continue;
	    AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data, lt );
	  }
	}
      }
    }
    else{
      for(int wire=0; wire<MaxWireBC4; ++wire){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdBC4, plane-NumOfLayersBc, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data =  gUnpacker.get( DetIdBC4, plane-NumOfLayersBc, 0, wire, lt, i );
	    if( data<MinBC4_TDC || MaxBC4_TDC<data ) continue;
	    AddDCRawHit( m_BcOutRawHC[plane+1], plane+PlMinBcOut, wire+1, data, lt );
	  }
	}
      }
    }
  }

  // SdcIn (SDC1)
  for( int plane=0; plane<NumOfLayersSDC1; ++plane ){
    for( int wire=0; wire<MaxWireSDC1; ++wire ){
      int nhit = gUnpacker.get_entries( DetIdSDC1, plane, 0, wire, 0 );
#if OscillationCut
      if( nhit>MaxMultiHitDC ) continue;
#endif
      for(int i=0; i<nhit; i++ ){
	int data = gUnpacker.get( DetIdSDC1, plane, 0, wire, 0, i ) ;
	if( data<MinSDC1_TDC || MaxSDC1_TDC<data ) continue;
	// AddDCRawHit( m_SdcInRawHC[plane+1], plane+PlMinSdcIn, wire+1, data );
	AddDCRawHit( m_SdcInRawHC[plane+1], plane+1, wire+1, data );
      }
    }
  }

  // SdcOut (SDC3&SDC4)
  for( int plane=0; plane<NumOfLayersSDC3+NumOfLayersSDC4; ++plane ){
    if( plane<NumOfLayersSDC3 ){
      for( int wire=0; wire<MaxWireSDC3; ++wire ){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdSDC3, plane, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC3, plane, 0, wire, lt, i );
	    if( lt == 0 && ( data<MinSDC3_TDC || MaxSDC3_TDC<data ) ) continue;
	    if( lt == 1 && data<MinSDC3_TDC ) continue;
	    //	    if((plane == 0 || plane == 1) && 53 < wire && wire < 65) continue;
	    //	    if((plane == 2 || plane == 3) && 61 < wire && wire < 68) continue;
	    AddDCRawHit( m_SdcOutRawHC[plane+1], plane+PlMinSdcOut, wire+1, data , lt);
	  }// for(i)
	}// for(lt)
      }// for(wire)
    }// if(SDC3/4)
    else{
      int MaxWireSDC4;
      if( plane==NumOfLayersSDC3 || plane==(NumOfLayersSDC3+1) )
	MaxWireSDC4 = MaxWireSDC4Y;
      else
	MaxWireSDC4 = MaxWireSDC4X;
      for( int wire=0; wire<MaxWireSDC4; ++wire ){
	for(int lt = 0; lt<2; ++lt){
	  int nhit = gUnpacker.get_entries( DetIdSDC4, plane-NumOfLayersSDC3, 0, wire, lt );
#if OscillationCut
	  if( nhit>MaxMultiHitDC ) continue;
#endif
	  for(int i=0; i<nhit; i++ ){
	    int data = gUnpacker.get( DetIdSDC4, plane-NumOfLayersSDC3, 0, wire, lt ,i );
	    if( lt == 0 && ( data<MinSDC4_TDC || MaxSDC4_TDC<data ) ) continue;
	    if( lt == 1 && data<MinSDC4_TDC ) continue;
	    //	    if((plane == 4 || plane == 5) && 30 < wire && wire < 38) continue;
	    //	    if((plane == 6 || plane == 7) && 44 < wire && wire < 54) continue;
	    AddDCRawHit( m_SdcOutRawHC[plane+1],  plane+PlMinSdcOut, wire+1, data , lt);
	  }// for(i)
	}// for(lt)
      }//for(wire)
    }// if(SDC3/4)
  }// for(plane)

  // Scaler
  for( int l = 0; l<NumOfScaler; ++l){
    for( int seg=0; seg<NumOfSegScaler; ++seg ){
      int nhit = gUnpacker.get_entries( DetIdScaler, l, 0, seg, 0 );
      if( nhit>0 ){
	int data = gUnpacker.get( DetIdScaler, l, 0, seg, 0 );
	AddHodoRawHit( m_ScalerRawHC, DetIdScaler, l, seg, 0, 0, data );
      }
    }
  }// for(l)

  // trigger Flag ---------------------------------------------------------------
  DecodeHodo( DetIdTrig, NumOfSegTrig, kOneSide, m_TrigRawHC );

  {
    // BH2Mt
    static const int type_mt = gUnpacker.get_data_id("BH2", "fpga_meantime");
    for(int seg = 0; seg<NumOfSegBH2; ++seg){
      int mhit = gUnpacker.get_entries( DetIdBH2, 0, seg, 0, type_mt );
      for(int m = 0; m<mhit; ++m){
	int data = gUnpacker.get( DetIdBH2, 0, seg, 0, type_mt , m);
	AddHodoRawHit( m_FpgaBH2MtRawHC, DetIdFpgaBH2Mt, 0, seg, 0, kHodoLeading, data );
      }// for(m)
    }// for(seg)
  }

  m_is_decoded = true;
  return true;
}

//______________________________________________________________________________
bool
RawData::DecodeCalibHits( void )
{
  del::ClearContainer( m_VmeCalibRawHC );

  for( int plane=0; plane<NumOfPlaneVmeCalib; ++plane ){
    DecodeHodo( DetIdVmeCalib, plane, NumOfSegVmeCalib,
  		kOneSide, m_VmeCalibRawHC );
  }

  return true;
}
