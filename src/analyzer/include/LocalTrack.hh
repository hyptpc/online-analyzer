/* LocalTrack.h */
#ifndef LocalTrack_h
#define LocalTrack_h 1

#include <iomanip>
#include <iostream>
#include <vector>
#include <map>

#include "DCHit.hh"
#include "TVector3.h"
#include "TString.h"

class LocalTrack
{
  
public:
  struct TrackHit
  {

    int cid;
    int layer;
    int wire;
    
    double wpos;
    double zpos;
    double dt;
    double dl;
    double angle;
    double hitpos;
    double trackpos;
    double residual;
    
    TrackHit( DCHit* hit )
    {
      cid=hit->GetDetId();
      zpos=hit->GetZ();
      wpos=hit->GetWirePosition();
      layer=hit->GetLayer();
      wire=hit->GetWire();
      dt=hit->GetDriftTime();
      dl=hit->GetDriftLength();
      angle=hit->GetTiltAngle();
    }
    inline void Print( void ) const
    {
      std::cout << std::setw(10)  << "TrackHit"
       		<< std::setw(5)  << cid
       		<< std::setw(5)  << layer
       		<< std::setw(5)  << wire
       		<< std::setw(10)  << zpos
       		<< std::setw(10)  << wpos
		<< std::setw(10) << dt << std::endl;
    }
  };

private:
  std::vector<TrackHit> m_hit[2];
  
public:
  LocalTrack();
  ~LocalTrack();
  
  int AddHit(DCHit *hit, int xy){
    TrackHit newhit(hit);
    CID=hit->GetDetId();
    m_hit[xy].push_back(newhit);
    //    newhit.Print();
    return nhit(xy);
  }
  int nhit( const int &xy )  { return m_hit[xy].size(); }
  double resid( const int &xy, int i )  { return m_hit[xy][i].residual; }
  int layer( const int &xy, int i )  { return m_hit[xy][i].layer; }
  
private:
  int CID;
  double A, B, C;		/* a track in XZ plane : Ax + Bz + C = 0 */
  double D, E, F;		/* a track in YZ plane : Dy + Ez + F = 0 */
  double GA, GB, GC;		/* Parameters for a global track */
  double GD, GE, GF;		/* Parameters for a global track */
  double GZ;

  int xzDof, yzDof;
  double xzChi, yzChi;

 public:  
  int cid()  const { return CID; }
  double a() const { return A; }
  double b() const { return B; }
  double c() const { return C; }
  double d() const { return D; }
  double e() const { return E; }
  double f() const { return F; }
  double x() const  { return -C/A; }
  double dx() const { return -B/A; }
  double y() const  { return -F/D; }
  double dy() const { return -E/D; }
  double ga() const { return GA; }
  double gb() const { return GB; }
  double gc() const { return GC; }
  double gd() const { return GD; }
  double ge() const { return GE; }
  double gf() const { return GF; }
double gx() const  { return -GC/GA; }
double gdx() const { return -GB/GA; }
  double gy() const  { return -GF/GD; }
  double gdy() const { return -GE/GD; }
  double gz() const { return GZ; }
  void labc( double &a, double &b, double &c ) { a=A; b=B; c=C; }
  void ldef( double &d, double &e, double &f ) { d=D; e=E; f=F; }
  void abc( double &a, double &b, double &c ) const { a=A; b=B; c=C; }
  void def( double &d, double &e, double &f ) const { d=D; e=E; f=F; }
  void gabc( double &a, double &b, double &c ) const { a=GA; b=GB; c=GC; }
  void gdef( double &d, double &e, double &f ) const { d=GD; e=GE; f=GF; }
  double chi2xz() const { return xzChi; }
  double chi2yz() const { return yzChi; }
  double chi2all() const { return (xzChi*xzDof + yzChi*yzDof)/(double)(xzDof+yzDof); }
  int dofxz() const { return xzDof; }
  int dofyz() const { return yzDof; }

void SetABC( const double &a, const double &b, const double &c ){ A=a; B=b; C=c; }
void SetDEF( const double &d, const double &e, const double &f ){ D=d; E=e; F=f; }
void SetChisqr( const int &xy, const double &val ) { 
    if(xy==0) xzChi=val; if(xy==1) yzChi=val; }
  void SetDof( const int &xy, const int &val ) {
    if(xy==0) xzDof=val; if(xy==1) yzDof=val; }
  double chi2(const int &xy) const { return xy ? yzChi : xzChi; }
  int dof(const int &xy) const { return xy ? yzDof : xzDof ; }

 public:
  bool XYLocalPosatZ( const double &z, double &x, double &y);
  bool XYPosatZ( const double &z, double &x, double &y );
  bool ZXPosatY( const double &y, double &z, double &x );
  bool ZYPosatX( const double &x, double &z, double &y );
  TVector3 GetLocalPosatZ(const double &z);
  TVector3 GetPosatZ(const double &z);
  TVector3 GetPosatX(const double &z);
  TVector3 GetPosatY(const double &z);
  TVector3 GetMomDir();
  TVector3 GetLocalMomDir();
  bool DoFit();
  void ConvLocalToGlobal(double rotation);
  bool LeastSquareFit( const int &xy, TString option="" );
  // bool LinearFit( BeamLineHitMan *blMan,ConfMan *conf, const bool &CHECKLR=false, TString option="" );
  // bool dRLineFit( BeamLineHitMan *blMan,ConfMan *conf, TString option="" );
  
public:
  // void Calculate();
  void Clear();
  void Print();

};

#endif
