// -*- C++ -*-

#ifndef GUI_PS_H
#define GUI_PS_H

#include <vector>

#include <TGFrame.h>
#include <TObject.h>

class TGCompositeFrame;
class TGCheckButton;
class TGTextButton;
class TGTextEntry;

namespace hddaq
{

namespace gui
{

//_____________________________________________________________________________
class GuiPs
{

public:
  enum e_command
    {
      k_Save,
      k_Print,
      k_nCommand
    };

private:
  TGCompositeFrame*           m_frame;
  std::vector<TGTextButton*>  m_command;
  std::vector<TGCheckButton*> m_optButton;
  std::vector<TString>        m_optList;
  std::vector<TGCheckButton*> m_devButton;
  std::vector<TString>        m_devList;
  TGTextEntry*                m_textFilename;
  TString                     m_filename;

public:
  static GuiPs& getInstance( void );
  virtual ~GuiPs();

  static TString getFilename( void );
  void initialize( const std::vector<TString>& optList,
                   const std::vector<TString>& devList );
  static bool isDevOn(int i);
  static bool isOptOn(int i);
  // Bool_t isAutoMode( void );
  void print( void );
  void save( void );
  void setFilename( void );
  void setFilename( const char* filename );
  void toggleAllDevice( void );
  void toggleAllOption( void );

private:
  GuiPs( void );
  GuiPs( const GuiPs& );
  GuiPs& operator=( const GuiPs& );

  ClassDef( hddaq::gui::GuiPs, 0 )
};

}

}
#endif
