// -*- C++ -*-

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "ConfMan.hh"
#include "FuncName.hh"
#include "GeAdcCalibMan.hh"

ClassImp(GeAdcCalibMan);

#define frand() (-0.5+(Double_t) rand() / RAND_MAX)

namespace
{

const Int_t BufSize = 200;

//______________________________________________________________________________
inline UInt_t MakeKey( Int_t seg )
{
  return seg;
}

}

//Adc1,Adc2,Adc3,Adc4,Adc5
//______________________________________________________________________________
struct GeAdcCalibMap
{
  GeAdcCalibMap( Double_t q0, Double_t q1, Double_t q2, Double_t q3, Double_t q4 )
    : p0(q0), p1(q1), p2(q2), p3(q3), p4(q4)
  {}
  Double_t p0, p1, p2, p3, p4;
};

//______________________________________________________________________________
GeAdcCalibMan::GeAdcCalibMan( void )
  : TObject()
{
}

//______________________________________________________________________________
GeAdcCalibMan::~GeAdcCalibMan( void )
{
  clearElements();
}

//______________________________________________________________________________
void
GeAdcCalibMan::clearElements( void )
{
  std::map <UInt_t, GeAdcCalibMap *>::iterator itr;
  for( itr=Cont_.begin(); itr!=Cont_.end(); ++itr ){
    delete itr->second;
    itr->second = 0;
  }
  Cont_.clear();
}

//______________________________________________________________________________
Bool_t
GeAdcCalibMan::Initialize( const TString& filename )
{
  ParamFileName_ = filename;

  Int_t seg;
  Double_t p0, p1, p2 ,p3, p4;
  FILE *fp;
  char buf[BufSize];

  if((fp=fopen(ParamFileName_.Data(),"r"))==0){
    std::cerr << FUNC_NAME << ": file open fail" << std::endl;
    return false;
  }

  while( fgets( buf, BufSize, fp ) ){
    if( buf[0]!='#' ){
      if( sscanf( buf, "%d %lf %lf %lf %lf %lf",
		  &seg, &p0, &p1, &p2, &p3, &p4 )==6 ){
	UInt_t key = MakeKey( seg );
	GeAdcCalibMap *p=new GeAdcCalibMap( p0, p1, p2, p3, p4 );
	if(p){
	  if( Cont_[key] ) delete Cont_[key];
	  Cont_[key]=p;
	}
	else{
	  std::cerr << FUNC_NAME << ": new fail. "
		    << " Seg=" << std::setw(3) << std::dec << seg
		    << std::endl;
	}
      }
      else{
	std::cerr << FUNC_NAME << ": Bad format => "
		  << TString(buf) << std::endl;
      }
    }
  }
  fclose(fp);

  std::cout << FUNC_NAME << ": Initialization finished" << std::endl;
  return true;
}

//______________________________________________________________________________
GeAdcCalibMap*
GeAdcCalibMan::getMap( Int_t seg ) const
{
  UInt_t key = MakeKey( seg );
  return Cont_[key];
}

//______________________________________________________________________________
Double_t
GeAdcCalibMan::CalibAdc( Int_t seg, Int_t adc, Double_t &a_ )
{
  GeAdcCalibMap *p=getMap( seg );
  if(p){
    Double_t adcCH=adc+frand();
    if( adc>0 ){
      a_=(p->p0)
	+(p->p1)*adcCH
	+(p->p2)*adcCH*adcCH
	+(p->p3)*adcCH*adcCH*adcCH
	+(p->p4)*adcCH*adcCH*adcCH*adcCH;
    }
    return true;
  }
  else{
    std::cerr << FUNC_NAME << ": No record. "
	      << " Seg=" << std::setw(3) << std::dec << seg
    	      << std::endl;
    return -1;
  }
}
