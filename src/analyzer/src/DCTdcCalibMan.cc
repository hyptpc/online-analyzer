// -*- C++ -*-

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCTdcCalibMan.hh"
#include "FuncName.hh"

//ClassImp(DCTdcCalibMan);

//______________________________________________________________________________
namespace
{
  inline UInt_t MakeKey( Int_t det, Int_t plane, Double_t wire )
  {
    return (det<<20) | (plane<<10) | Int_t(wire);
  }
}

//______________________________________________________________________________
struct DCTdcCalMap
{
  DCTdcCalMap( Double_t q0, Double_t q1 )
    : p0(q0), p1(q1)
  {}
  Double_t p0, p1;
};

// initialize DCTdcCalibMan --------------------------------------------------
void
ConfMan::initializeDCTdcCalibMan()
{
  if(name_file_["DCTDC:"] != ""){
    DCTdcCalibMan& gDCTdcCalib = DCTdcCalibMan::GetInstance();
    gDCTdcCalib.SetFileName(name_file_["DCTDC:"]);
    flag_[kIsGood] = gDCTdcCalib.Initialize();
  }else{
    std::cout << "#E ConfMan::"
	      << " File path does not exist in " << name_file_["DCTDC:"] 
	      << std::endl;
    flag_.reset(kIsGood);
  }
}
// initialize DCTdcCalibMan --------------------------------------------------
//______________________________________________________________________________
DCTdcCalibMan::DCTdcCalibMan( void )
  : TObject(),
    m_is_ready(false),
    m_file_name(),
    m_map()
{
}

//______________________________________________________________________________
DCTdcCalibMan::~DCTdcCalibMan( void )
{
  ClearElements();
}

//______________________________________________________________________________
void
DCTdcCalibMan::ClearElements( void )
{
  std::map <UInt_t, DCTdcCalMap *>::iterator itr;
  for( itr=m_map.begin(); itr!=m_map.end(); ++itr ){
    delete itr->second;
    itr->second = 0;
  }
  m_map.clear();
}

//______________________________________________________________________________
Bool_t
DCTdcCalibMan::Initialize( void )
{
  Int_t cid, pid, wid;
  Double_t p0, p1;

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line[0] == '#' ) continue;
    if(TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream iss( line );
    if( iss >> cid >> pid >> wid >> p0 >> p1 ){
      UInt_t key = MakeKey( cid, pid, wid );
      DCTdcCalMap *p = new DCTdcCalMap( p0, p1 );
      if(p){
	if(m_map[key]) delete m_map[key];
	m_map[key]=p;
      }
      else{
	std::cerr << FUNC_NAME << ": new fail. "
		  << " Det  =" << std::setw(3) << std::dec << cid
		  << " Plane=" << std::setw(3) << std::dec << pid
		  << " Wire =" << std::setw(3) << std::dec << wid
		  << std::endl;
      }
    }
    else{
      std::cerr << FUNC_NAME << ": Bad format => "
      		<< line << std::endl;
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
DCTdcCalMap*
DCTdcCalibMan::GetMap( Int_t DetId, Int_t PlaneId, Double_t WireId ) const
{
  UInt_t key = MakeKey( DetId, PlaneId, WireId );
  std::map<UInt_t, DCTdcCalMap*>::const_iterator itr = m_map.find(key);
  if( itr != m_map.end() )
    return itr->second;
  else
    return nullptr;
}

//______________________________________________________________________________
Bool_t
DCTdcCalibMan::GetTime( Int_t DetId, Int_t PlaneId, Double_t WireId,
			Int_t tdc, Double_t & time ) const
{
  DCTdcCalMap *p = GetMap(DetId,PlaneId,WireId);
  if(p){
    time=(p->p0)+(p->p1)*tdc;
    return true;
  }
  p = GetMap(DetId,PlaneId,0);
  if(p){
    time=(p->p0)+(p->p1)*tdc;
    return true;
  }
  else{
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " DetId  =" << std::setw(3) << std::dec << DetId
		<< " PlaneId=" << std::setw(3) << std::dec << PlaneId
		<< " WireId=" << std::setw(3) << std::dec << WireId
		<< std::endl;
    // return false;
  }
}

//______________________________________________________________________________
Bool_t
DCTdcCalibMan::GetTdc( Int_t DetId, Int_t PlaneId, Double_t WireId,
		       Double_t time, Int_t &tdc ) const
{
  DCTdcCalMap *p = GetMap( DetId, PlaneId, WireId );
  if(p){
    tdc = (Int_t)((time-(p->p0))/(p->p1));
    return true;
  }
  p = GetMap( DetId, PlaneId, 0 );
  if(p){
    tdc = (Int_t)((time-(p->p0))/(p->p1));
    return true;
  }
  else{
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " DetId  =" << std::setw(3) << std::dec << DetId
		<< " PlaneId=" << std::setw(3) << std::dec << PlaneId
		<< " WireId=" << std::setw(3) << std::dec << WireId
		<< std::endl;
    // return false;
  }
}
