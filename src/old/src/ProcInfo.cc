// -*- C++ -*-

// Author: Shuhei Hayakawa

#include <iostream>

#include "FuncName.hh"
#include "ProcInfo.hh"

//______________________________________________________________________________
ProcInfo::ProcInfo( void )
{
}

//______________________________________________________________________________
ProcInfo::~ProcInfo( void )
{
}

//______________________________________________________________________________
void
ProcInfo::CheckMemory( void )
{
  ProcInfo_t p;
  gSystem->GetProcInfo( &p );

  std::cout << "#D " << FUNC_NAME << std::endl;
  Printf("   Rmem = %8.3f MB, Vmem = %8.3f MB \n",
	 p.fMemResident * 1.e-3,
	 p.fMemVirtual  * 1.e-3 );
}
