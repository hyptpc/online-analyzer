// -*- C++ -*-

#ifndef MATRIX_PARAM_MAN_HH
#define MATRIX_PARAM_MAN_HH

#include <vector>
#include <map>

#include <TString.h>

//______________________________________________________________________________
class MatrixParamMan
{
public:
  static MatrixParamMan& GetInstance( void );
  static TString&        ClassName( void );
  ~MatrixParamMan( void );

private:
  MatrixParamMan( void );
  MatrixParamMan( const MatrixParamMan& );
  MatrixParamMan& operator =( const MatrixParamMan& );

private:
  typedef std::vector< std::vector<Int_t> > Matrix2D;
  typedef std::vector< std::vector< std::vector<Int_t> > > Matrix3D;
  Bool_t   m_is_ready;
  TString  m_file_name_2d1;
  TString  m_file_name_2d2;
  TString  m_file_name_3d;
  Matrix2D m_enable_2d1;
  Matrix2D m_enable_2d2;
  Matrix3D m_enable_3d;

public:
  Bool_t Initialize( void );
  Bool_t Initialize( const TString& filename_2d1,
                     const TString& filename_2d2,
                     const TString& filename_3d );
  Bool_t IsAccept2D1( UInt_t detA, UInt_t detB ) const;
  Bool_t IsAccept2D2( UInt_t detA, UInt_t detB ) const;
  Bool_t IsAccept3D( UInt_t detA, UInt_t detB, UInt_t detC ) const;
  Bool_t IsReady( void ) const { return m_is_ready; }
  void   Print2D1( const TString& arg="" ) const;
  void   Print2D2( const TString& arg="" ) const;
  void   Print3D( const TString& arg="" ) const;
  void   SetMatrix2D1( const TString& file_name );
  void   SetMatrix2D2( const TString& file_name );
  void   SetMatrix3D( const TString& file_name );

};

//______________________________________________________________________________
inline MatrixParamMan&
MatrixParamMan::GetInstance( void )
{
  static MatrixParamMan s_instance;
  return s_instance;
}

//______________________________________________________________________________
inline TString&
MatrixParamMan::ClassName( void )
{
  static TString s_name( "MatixParamMan" );
  return s_name;
}

#endif
