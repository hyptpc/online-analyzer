// -*- C++ -*-

#include <string>

#include <TString.h>

#include <escape_sequence.hh>
#include <std_ostream.hh>

#include "Exception.hh"
#include "FuncName.hh"

ClassImp(Exception);

//______________________________________________________________________________
Exception::Exception( const TString& msg )
  : TObject(),
    m_msg()
{
  m_msg = hddaq::unpacker::esc::k_yellow
    + FUNC_NAME + " " + msg
    + hddaq::unpacker::esc::k_default_color;
}

//______________________________________________________________________________
Exception::~Exception( void ) throw()
{
}

//______________________________________________________________________________
void
Exception::hoge( const TString& arg ) const
{
  hddaq::cout << "#D " << FUNC_NAME << " " << arg << std::endl;
}

//______________________________________________________________________________
const Char_t*
Exception::what( void ) const throw()
{
  return m_msg.Data();
}
