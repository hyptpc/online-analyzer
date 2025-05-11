// -*- C++ -*-

#ifndef PS_MAKER_HH
#define PS_MAKER_HH

#include <vector>
#include <TString.h>

class TCanvas;
class TPostScript;

//_____________________________________________________________________________
class PsMaker
{
public:
  static PsMaker& getInstance( void );
  virtual ~PsMaker( void );

private:
  PsMaker( void );
  PsMaker( const PsMaker& object );
  PsMaker& operator=( const PsMaker& object );

private:
  enum EPSType {
    kPortrait  = 111,
    kLandscape = 112,
    NPSType
  };
  enum ParameterList {
    kXdiv,
    kYdiv,
    kXrange_min,
    kXrange_max,
    sizeParameterList
  };

  std::vector<TString> m_name_option;
  TCanvas*             m_canvas;
  TPostScript*         m_ps;

public:
  enum OptionList {
    kExpDataSheet,
    kFixXaxis,
    kLogyADC,
    kLogyTDC,
    // kLogyHitPat,
    // kLogyMulti,
    // kAutoSaveAtRunChange,
    sizeOptionList
  };
  void getListOfOption( std::vector<TString>& vec );
  void makePs( void );

private:
  void drawRunNumber( void );
  void create( TString& name );
  void drawOneCanvas( std::vector<Int_t>& id_list,
                      std::vector<Int_t>& par_list,
                      Bool_t flag_xaxis, Bool_t flag_log,
                      const Option_t* optDraw = "" );
  void drawDCEff( void );
  void clearOneCanvas( Int_t npad );
};

//_____________________________________________________________________________
inline PsMaker&
PsMaker::getInstance( void )
{
  static PsMaker s_instance;
  return s_instance;
}

#endif
