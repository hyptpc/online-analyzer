// -*- C++ -*-

// Author: Shuhei Hayakawa

#include <TSystem.h>

//______________________________________________________________________________
class ProcInfo
{
public:
  ProcInfo( void );
  ~ProcInfo( void );

private:
  ProcInfo( const ProcInfo& );
  ProcInfo& operator =( const ProcInfo& );

private:

public:
  static TString& ClassName( void );
  void            CheckMemory( void );

};

//______________________________________________________________________________
inline TString&
ProcInfo::ClassName( void )
{
  static TString g_name("ProcInfo");
  return g_name;
}
