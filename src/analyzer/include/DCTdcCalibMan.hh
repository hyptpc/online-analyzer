// -*- C++ -*-

#ifndef DC_TDC_CALIB_MAN_HH
#define DC_TDC_CALIB_MAN_HH

#include <map>
#include <string>
#include <vector>

#include <TObject.h>
#include <TString.h>

struct DCTdcCalMap;

//______________________________________________________________________________
class DCTdcCalibMan : public TObject
{
public:
  static DCTdcCalibMan& GetInstance( void );
  ~DCTdcCalibMan( void );

private:
  DCTdcCalibMan( void );
  DCTdcCalibMan( const DCTdcCalibMan& );
  DCTdcCalibMan& operator =( const DCTdcCalibMan& );

private:
  Bool_t                         m_is_ready;
  TString                        m_file_name;
  std::map<UInt_t, DCTdcCalMap*> m_map;

public:
  Bool_t Initialize( void );
  Bool_t Initialize( const TString& file_name ) { SetFileName(file_name); return Initialize(); }
  Bool_t IsReady( void ) const { return m_is_ready; }
  void   SetFileName( const TString& file_name ){ m_file_name = file_name; }
  Bool_t GetTime( Int_t DetId, Int_t PlaneId, Double_t WireId, Int_t tdc, Double_t & time ) const;
  Bool_t GetTdc( Int_t DetId, Int_t PlaneId, Double_t WireId, Double_t time, Int_t & tdc ) const;

private:
  DCTdcCalMap* GetMap( Int_t DetId, Int_t PlaneId, Double_t WireId ) const;
  void         ClearElements( void );

  //  ClassDef(DCTdcCalibMan,0);
};

//______________________________________________________________________________
inline DCTdcCalibMan&
DCTdcCalibMan::GetInstance( void )
{
  static DCTdcCalibMan g_instance;
  return g_instance;
}

#endif
