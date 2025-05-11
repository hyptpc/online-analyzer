// BeamLineTrackMan.cpp
#include <map>

#include "LocalTrack.hh"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
static const double SpatialResolutionOfBLDC=0.02; // [cm]

// ----------------------------- //
// class LocalTrack              //
// ------------------------------//

LocalTrack::LocalTrack()
{
  m_hit[0].clear(); 
  m_hit[1].clear(); 
  A=B=C=D=E=F=-999.;
  GA=GB=GC=GD=GE=GF=-999.;
  xzDof=yzDof=0;
  xzChi=yzChi=-999.;
  CID=-1;
}

LocalTrack::~LocalTrack()
{
}
void LocalTrack::Clear()
{
  m_hit[0].clear();
  m_hit[1].clear();
}

bool LocalTrack::XYLocalPosatZ( const double &z, double &x, double &y)
{
  if( fabs(A) < 1.0e-9 ) return false;
  if( fabs(D) < 1.0e-9 ) return false;

  double tmpB=B,tmpE=E;
  x = -1.*( tmpB*z + C )/A;
  y = -1.*( tmpE*z + F )/D;
  return true;
}

bool LocalTrack::XYPosatZ( const double &z, double &x, double &y )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  x = -1.*( GB*(z-GZ) + GC )/GA;
  y = -1.*( GE*(z-GZ) + GF )/GD;   
  return true;
}

bool LocalTrack::ZYPosatX( const double &x, double &z, double &y )
{
  if( fabs(GB) < 1.0e-9 ) return false;
  if( fabs(GD) < 1.0e-9 ) return false;

  z = -1.*( GA*x +GC ) / GB + GZ;
  y = -1.*( GE*(z-GZ) + GF )/GD;     
  return true;
}

bool LocalTrack::ZXPosatY( const double &y, double &z, double &x )
{
  if( fabs(GA) < 1.0e-9 ) return false;
  if( fabs(GE) < 1.0e-9 ) return false;
  
  z = -1.*( GD*y +GF ) / GE + GZ;
  x = -1.*( GB*(z-GZ) + GC )/GA;
  return true;
}

void LocalTrack::ConvLocalToGlobal(double rotation)  
{
  //  std::cout<<"convL2G nhit "<<nhit(0)<<std::endl;
  TVector3 pos(-C/A,-F/D,0);
  pos.RotateZ(rotation*TMath::DegToRad());
  TVector3 dir(-B/A,-E/D,1);
  dir.RotateZ(rotation*TMath::DegToRad());
  GA = 1.;
  GB = -dir.X()/dir.Z();
  GC = -(pos.X()-dir.X()/dir.Z()*pos.Z());
  GD = 1.; 
  GE = -dir.Y()/dir.Z();
  GF = -(pos.Y()-dir.Y()/dir.Z()*pos.Z());
  GZ = 0;
}

TVector3 LocalTrack::GetLocalPosatZ(const double &z)
{
  double tmpx,tmpy;
  XYLocalPosatZ(z,tmpx,tmpy);
  return TVector3(tmpx,tmpy,z);
}
TVector3 LocalTrack::GetPosatZ(const double &z)
{
  double tmpx,tmpy;
  XYPosatZ(z,tmpx,tmpy);
  return TVector3(tmpx,tmpy,z);
}
TVector3 LocalTrack::GetPosatY(const double &y)
{
  double tmpx,tmpz;
  ZXPosatY(y,tmpz,tmpx);
  return TVector3(tmpx,y,tmpz);
}

TVector3 LocalTrack::GetPosatX(const double &x)
{
  double tmpy,tmpz;
  ZXPosatY(x,tmpz,tmpy);
  return TVector3(x,tmpy,tmpz);
}
TVector3 LocalTrack::GetMomDir()
{
  TVector3 tmpdir(-GB/GA,-GE/GD,1);
  return tmpdir.Unit();
}
TVector3 LocalTrack::GetLocalMomDir()
{
  TVector3 tmpdir(-B/A,-E/D,1);
  return tmpdir.Unit();
}


static const unsigned int LRMASK=0x01;
inline static unsigned int LR( int hid, unsigned int key )
{
  key = key >> hid;
  key = key & LRMASK;
  return key;
}
//####################################################
// bool LocalTrack::dRLineFit( BeamLineHitMan *blMan,ConfMan *conf, TString option )
// {
//   BLDCFittingParamMan *BLDCParam=conf->GetBLDCFittingParamManager();
//   int MaxNumOfHitsInTrack=BLDCParam->GetMaxHitInTrack();
//   if( nhit() > MaxNumOfHitsInTrack )
//     return false;

//   TVector3 wirepos[20];
//   TVector3 wiredir[20];
//   double weight[MaxNumOfHitsInTrack];
//   double drift[MaxNumOfHitsInTrack];
  
//   int allhit=0;
//   double sigma=SpatialResolutionOfBLDC;
//   if(option.Contains("mwpc")) sigma=hit(blMan,0)->dxy()/sqrt(12);
//   for(int n=0;n<(int)nhit();n++)
//     {
//       ChamberLikeHit *tmphit=hit(blMan,n);
//       tmphit->wposdir(wirepos[allhit],wiredir[allhit]);
//       weight[allhit]=sigma;
//       drift[allhit]=tmphit->dl();
//       allhit++;
//     }      

//   double tmppar[6]={x(),y(),0,
// 		    dx(),dy(),1};
//   LineFit *fit=new LineFit(tmppar,wirepos,wiredir,weight,drift,allhit);
//   fit->GetParameters(tmppar); 
//   xzChi=fit->chisquare();
//   xzDof=fit->dof();
//   delete fit;

//   TVector3 pos(tmppar[0],tmppar[1],tmppar[2]);
//   TVector3 dir(tmppar[3],tmppar[4],tmppar[5]);
//   dir=1./dir.Z()*dir;
//   pos=pos-dir*pos.Z();
//   A=1;
//   B=-dir.X();
//   C=-pos.X();
//   D=1;
//   E=-dir.Y();
//   F=-pos.Y();
//   return true;
// }

bool LocalTrack::DoFit()
{
  return LeastSquareFit(0) &&  LeastSquareFit(1);
}
bool LocalTrack::LeastSquareFit( const int &xy, TString option )
{
#if DEBUG
  std::cout << "!!! LocalTrack::LeastSquareFit()  "<<option << std::endl;
#endif
  int np = nhit(xy);
  if( np < 2 ){
#if DEBUG 
    std::cout<<"LocalTrack::LeastSquareFit():  too few hits "<<std::endl;
#endif 
    return false;
  }
  double sigma=SpatialResolutionOfBLDC;
  double minchi = 1.0e+9;
  double canda=0., candb=0., candc=0.;
  unsigned int candkey=0x0;
  unsigned int key = 0x0;
  int nconb = (int)pow(2,np);
  if(option.Contains("mwpc")) nconb=1;  

  double x[16],y[16],w[16];
  double tmp,tmp2;
  for( int ic=0; ic<nconb; ic++ ){      
    int j=0;
    for( int i=0; i<np;i++){
      x[j] = m_hit[xy][i].zpos;
      y[j] = m_hit[xy][i].wpos;
      if(!option.Contains("mwpc")){
	w[j]=sigma;
	double dl = m_hit[xy][i].dl;
	unsigned int lr = LR(j,key);
	if( lr==0 ) y[j] += dl;
	else y[j] -= dl;
      }else{
	w[j]=sigma;
      }
      j++;
    }
    double Sx, Sy, Sxy, Sxx;
    Sx=Sy=Sxy=Sxx=0;
    for( int i=0; i<np; i++ ){
      Sx  += x[i];
      Sy  += y[i];
      Sxy += x[i]*y[i];
      Sxx += x[i]*x[i];
    }
    double D = (double)np*Sxx - Sx*Sx;
    if( D==0 ){
#if DEBUG 
      std::cout<<"LocalTrack::LeastSquareFit():  invalid determinal "<<std::endl;
#endif 
      return false;
    }
    double aa = (Sxx*Sy - Sx*Sxy)/D;
    double bb = ((double)np*Sxy - Sx*Sy)/D;
    double chi = 0;
    for( int i=0; i<np; i++ ) chi += (y[i]-aa-bb*x[i])*(y[i]-aa-bb*x[i])/(w[i]*w[i]);    
    if(np<3) chi=-1;
    else chi /= (double)(np-2);
    //    std::cout<<ic<<" / "<< nconb<<" key:"<<key<<" candkey:"<<candkey<<"  chi:"<<chi<< "  minchi:"<<minchi<<std::endl;
    if( (chi < minchi) || (chi==minchi && candb*candb<bb*bb) ){
      canda = 1.;
      candb = -1.*bb;
      candc = -1.*aa;
      candkey=key;
      minchi=chi;
    }
    key++;
  }
  if(minchi==1.0e+9){
#if DEBUG 
    std::cout<<"LocalTrack::LeastSquareFit():  no track candidate "<<std::endl;
#endif 
    return false;
  }
  if(xy)
    SetDEF(canda,candb,candc);
  else
    SetABC(canda,candb,candc);
  SetDof(xy,np-2);
  SetChisqr(xy,minchi);
  if(option.Contains("mwpc")) return true;
  int j=0;
  //  std::cout<<"-----------------------------__"<<std::endl;
  for( int i=0; i<np; i++ ){
    double dl = m_hit[xy][i].dl;
    double zpos= m_hit[xy][i].zpos; 
    double wpos = m_hit[xy][i].wpos;
    double tmpx,tmpy;
    unsigned int lr = LR(i,candkey);
    if( lr==0 ) m_hit[xy][i].hitpos= wpos + dl;
    else        m_hit[xy][i].hitpos= wpos - dl;
    double tmppos[2];
    XYLocalPosatZ(m_hit[xy][i].zpos,tmppos[0],tmppos[1]);
    m_hit[xy][i].trackpos=tmppos[xy];
    m_hit[xy][i].residual=m_hit[xy][i].hitpos - tmppos[xy];
    //    std::cout<<xy<<"  "<<i<<" / "<<np<<"  z:"<<zpos<<"  w:"<<wpos<<" lr:"<<lr<<" dl: "<<dl<<" trapos: "<<tmppos[xy]<<" hitpos:"<<m_hit[xy][i].hitpos<<"  resid:"<<m_hit[xy][i].residual<<std::endl;
  }
  return true;
}
void LocalTrack::Print(){
  // std::cout<< "DC_ID:" << this->cid() <<std::endl;
  // std::cout<<"cid,layer,wire,dt,cdt,dl,lr,wx,wy,x,y "<<std::endl;
  // for(int i=0;i<this->nhit();i++){
  //   ChamberLikeHit *hit=this->hit(blMan,i);
  //   std::cout<< hit->cid() << "\t"
  // 	     << hit->layer() << "\t"
  // 	     << hit->wire() << "\t"
  // 	     << hit->dt() << "\t"
  // 	     << hit->cdt() << "\t"
  // 	     << hit->dl() << "\t"
  // 	     << hit->leftright() << "\t"
  // 	     << hit->wx() << "\t"
  // 	     << hit->wy() << "\t"
  // 	     << hit->x() << "\t"
  // 	     << hit->y() << "\t"
  // 	     <<std::endl;
  //   TVector3 pos,dir;
  //   hit->wposdir(pos,dir);
  //   pos.Print();
  //   dir.Print();
  // }

  // std::cout<< "x: Chi2,Dof,a,b,c :" << this->chi2xz() << "\t" << this->dofxz() 
  // 	   << "\t" << this->a() << "\t" << this->b() << "\t" << this->c() << std::endl;
  // std::cout<< "           ga,b,c :" << this->chi2xz() << "\t" << this->dofxz() 
  // 	   << "\t" << this->ga() << "\t" << this->gb() << "\t" << this->gc() << std::endl;
  // std::cout<< "y: Chi2,Dof,a,b,c :" << this->chi2yz() << "\t" << this->dofyz() 
  // 	   << "\t" << this->d() << "\t" << this->e() << "\t" << this->f() << std::endl;
  // std::cout<< "           ga,b,c :" << this->chi2yz() << "\t" << this->dofyz() 
  // 	   << "\t" << this->gd() << "\t" << this->ge() << "\t" << this->gf() << std::endl;
  // std::cout<<"------------------------------------------------------"<<std::endl;
}

