// -*- C++ -*-

#include "MatrixParamMan.hh"

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "ThreeVector.hh"

namespace
{
}

//______________________________________________________________________________
MatrixParamMan::MatrixParamMan( void )
  : m_is_ready( false ),
    m_file_name_2d1(),
    m_file_name_2d2(),
    m_file_name_3d()
{
  // 2D
  m_enable_2d1.resize( NumOfSegTOF_Mtx );
  m_enable_2d2.resize( NumOfSegTOF_Mtx );
  for( int i=0; i<NumOfSegTOF_Mtx; ++i ){
    m_enable_2d1[i].resize( NumOfSegSCH );
    m_enable_2d2[i].resize( NumOfSegSCH );
  }
  // 3D
  m_enable_3d.resize( NumOfSegTOF_Mtx );
  for( int i=0; i<NumOfSegTOF_Mtx; ++i ){
    m_enable_3d[i].resize( NumOfSegSCH );
    for( int j=0; j<NumOfSegSCH; ++j ){
      m_enable_3d[i][j].resize( NumOfSegBH2 );
    }
  }
}

//______________________________________________________________________________
MatrixParamMan::~MatrixParamMan( void )
{
}

//______________________________________________________________________________
Bool_t
MatrixParamMan::Initialize( void )
{
  return Initialize( m_file_name_2d1, m_file_name_2d2, m_file_name_3d );
}

//______________________________________________________________________________
Bool_t
MatrixParamMan::Initialize( const TString& filename_2d1,
                            const TString& filename_2d2,
			    const TString& filename_3d )
{
  { // 2D1
    std::ifstream ifs( filename_2d1.Data() );
    if( !ifs.is_open() ){
      hddaq::cerr << "#E " << FUNC_NAME << " "
                  << "No such 2D1 parameter file : "
                  << filename_2d1 << std::endl;
      return false;
      // std::exit( EXIT_FAILURE );
    }

    TString line;
    Int_t tofseg = 0;
    while( ifs.good() && line.ReadLine( ifs ) ){
      if( line.IsNull() || line[0] == '#' ) continue;
      TString param[2];
      std::istringstream iss( line.Data() );
      iss >> param[0] >> param[1];
      Int_t schseg = param[0].ReplaceAll( "SCH", "" ).Atoi();
      Int_t enable = param[1].Atoi();
      m_enable_2d1[tofseg][schseg] = enable;
      if( schseg == NumOfSegSCH-1 ) ++tofseg;
    }
  }
  { // 2D2
    std::ifstream ifs( filename_2d2.Data() );
    if( !ifs.is_open() ){
      hddaq::cerr << "#E " << FUNC_NAME << " "
                  << "No such 2D2 parameter file : "
                  << filename_2d2 << std::endl;
      return false;
      // std::exit( EXIT_FAILURE );
    }

    TString line;
    Int_t tofseg = 0;
    while( ifs.good() && line.ReadLine( ifs ) ){
      if( line.IsNull() || line[0] == '#' ) continue;
      TString param[2];
      std::istringstream iss( line.Data() );
      iss >> param[0] >> param[1];
      Int_t schseg = param[0].ReplaceAll( "SCH", "" ).Atoi();
      Int_t enable = param[1].Atoi();
      m_enable_2d2[tofseg][schseg] = enable;
      if( schseg == NumOfSegSCH-1 ) ++tofseg;
    }
  }
  { // 3D
    std::ifstream ifs( filename_3d.Data() );
    if( !ifs.is_open() ){
      hddaq::cerr << "#E " << FUNC_NAME << " "
		  << "No such 3D parameter file : "
                  << filename_3d << std::endl;
      return false;
      // std::exit( EXIT_FAILURE );
    }

    TString line;
    Int_t tofseg = 0;
    while( ifs.good() && line.ReadLine( ifs ) ){
      if( line.IsNull() || line[0] == '#' ) continue;
      TString param[2];
      std::istringstream iss( line.Data() );
      iss >> param[0] >> param[1];
      Int_t schseg = param[0].ReplaceAll( "SCH", "" ).Atoi();
      for( Int_t i=0; i<NumOfSegBH2; ++i ){
        Int_t enable = TString( param[1][i] ).Atoi();
        m_enable_3d[tofseg][schseg][i] = enable;
      }
      if( schseg == NumOfSegSCH-1 ) ++tofseg;
    }
  }
  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
MatrixParamMan::IsAccept2D1( UInt_t detA, UInt_t detB ) const
{
  if( m_enable_2d1.size() <= detA ){
    hddaq::cerr << "#W " << FUNC_NAME << " detA is too much : "
		<< detA << "/" << m_enable_2d1.size() << std::endl;
    return false;
  }

  if( m_enable_2d1[detA].size() <= detB ){
    hddaq::cerr << "#W " << FUNC_NAME << " detB is too much : "
		<< detB << "/" << m_enable_2d1[detA].size() << std::endl;
    return false;
  }
  return ( m_enable_2d1[detA][detB] == 1 );
}

//______________________________________________________________________________
Bool_t
MatrixParamMan::IsAccept2D2( UInt_t detA, UInt_t detB ) const
{
  if( m_enable_2d2.size() <= detA ){
    hddaq::cerr << "#W " << FUNC_NAME << " detA is too much : "
		<< detA << "/" << m_enable_2d2.size() << std::endl;
    return false;
  }

  if( m_enable_2d2[detA].size() <= detB ){
    hddaq::cerr << "#W " << FUNC_NAME << " detB is too much : "
		<< detB << "/" << m_enable_2d2[detA].size() << std::endl;
    return false;
  }
  return ( m_enable_2d2[detA][detB] == 1 );
}

//______________________________________________________________________________
Bool_t
MatrixParamMan::IsAccept3D( UInt_t detA, UInt_t detB, UInt_t detC ) const
{
  if( m_enable_3d.size() <= detA ){
    hddaq::cerr << "#W " << FUNC_NAME << " detA is too much : "
		<< detA << "/" << m_enable_3d.size() << std::endl;
    return false;
  }

  if( m_enable_3d[detA].size() <= detB ){
    hddaq::cerr << "#W " << FUNC_NAME << " detB is too much : "
		<< detB << "/" << m_enable_3d[detA].size() << std::endl;
    return false;
  }

  if( m_enable_3d[detA][detB].size() <= detC ){
    hddaq::cerr << "#W " << FUNC_NAME << " detC is too much : "
		<< detC << "/" << m_enable_3d[detA][detB].size() << std::endl;
    return false;
  }
  return ( m_enable_3d[detA][detB][detC] == 1 );
}

//______________________________________________________________________________
void
MatrixParamMan::Print2D1( const TString& arg ) const
{
  hddaq::cout << "#D " << FUNC_NAME << " " << arg << std::endl;
  for( Int_t i=NumOfSegTOF_Mtx-1; i>=0; --i ){
    hddaq::cout << " detA = " << std::setw(2)
		<< std::right << i << " : ";
    for( Int_t j=0; j<NumOfSegSCH; ++j ){
      hddaq::cout << m_enable_2d1[i][j];
    }
    hddaq::cout << std::endl;
  }
}

//______________________________________________________________________________
void
MatrixParamMan::Print2D2( const TString& arg ) const
{
  hddaq::cout << "#D " << FUNC_NAME << " " << arg << std::endl;
  for( Int_t i=NumOfSegTOF_Mtx-1; i>=0; --i ){
    hddaq::cout << " detA = " << std::setw(2)
		<< std::right << i << " : ";
    for( Int_t j=0; j<NumOfSegSCH; ++j ){
      hddaq::cout << m_enable_2d2[i][j];
    }
    hddaq::cout << std::endl;
  }
}

//______________________________________________________________________________
void
MatrixParamMan::Print3D( const TString& arg ) const
{
  hddaq::cout << "#D " << FUNC_NAME << " " << arg << std::endl;
  for( Int_t k=0; k<NumOfSegBH2; ++k ){
    hddaq::cout << " detC = " << std::setw(2)
		<< std::right << k << std::endl;
    for( Int_t i=NumOfSegTOF_Mtx-1; i>=0; --i ){
      hddaq::cout << " detA = " << std::setw(2)
		  << std::right << i << " : ";
      for( Int_t j=0; j<NumOfSegSCH; ++j ){
	hddaq::cout << m_enable_3d[i][j][k];
      }
      hddaq::cout << std::endl;
    }
  }
}

//______________________________________________________________________________
void
MatrixParamMan::SetMatrix2D1( const TString& file_name )
{
  m_file_name_2d1 = file_name;
}

//______________________________________________________________________________
void
MatrixParamMan::SetMatrix2D2( const TString& file_name )
{
  m_file_name_2d2 = file_name;
}

//______________________________________________________________________________
void
MatrixParamMan::SetMatrix3D( const TString& file_name )
{
  m_file_name_3d = file_name;
}
