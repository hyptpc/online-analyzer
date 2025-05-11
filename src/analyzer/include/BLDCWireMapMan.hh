// BLDCWireMapMan.h

#ifndef BLDCWireMapMan_h
#define BLDCWireMapMan_h 1

#include <map>
#include <string>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

#include "TVector3.h"
#include "TMath.h"
#include "TString.h"
//#include "GeomMapMan.h"

class GeomMap
{
 public:
  GeomMap(){};
  ~GeomMap() {};
 private:
  static const int npar=11;
  double R,Theta;
  double Pos[3];
  double Rot[3];
  double Size[4];
  double LightVelocity;
  
public:
  inline void SetParam(const double *par);
  inline void GetParam(double *par);
  void SetPos( const TVector3 &pos ) { Pos[0]=pos.X(),Pos[1]=pos.Y(),Pos[2]=pos.Z(); }          
  void SetX( const double &x ) { Pos[0]=x; }
  void SetY( const double &y ) { Pos[1]=y; }
  void SetZ( const double &z ) { Pos[2]=z; }
  void SetRTheta( const double &r, const double &theta)
  { R=r; Theta=theta; Pos[0]=r*TMath::Cos(theta*TMath::DegToRad());  Pos[1]=r*TMath::Sin(theta*TMath::DegToRad()); }
  void SetRot( const TVector3 &rot ) { Rot[0]=rot.X(),Rot[1]=rot.Y(),Rot[2]=rot.Z(); }          
  void SetRotX( const double &x ) { Rot[0]=x; }
  void SetRotY( const double &y ) { Rot[1]=y; }
  void SetRotZ( const double &z ) { Rot[2]=z; }
  void SetSize( const double *size ) { for( int i =0; i<4;i++ ) Size[i]=size[i]; }
  void SetLength( const double &len ) { Size[1] = len; }
  void SetWidth( const double &w ) { Size[0] = w; }
  void SetThick( const double &t ) { Size[2] = t; }
  void SetLightVelocity( const double &lv ) { LightVelocity = lv; }

  TVector3 GetPos() { return TVector3(Pos); }
  inline TVector3 GetPos(const double &ctsub);
  TVector3 GetRot() { return TVector3(Rot); }
  double* const GetSize() { return Size; }

  double GetX() { return Pos[0]; }
  double GetY() { return Pos[1]; }
  double GetZ() { return Pos[2]; }
  double GetR() { return R; }
  double GetTheta() { return Theta; }

  double GetRotX() { return Rot[0]; }
  double GetRotY() { return Rot[1]; }
  double GetRotZ() { return Rot[2]; }

  double GetLength() { return IsBox() ? Size[1] : Size[3]; }
  double GetWidth() { return IsBox() ? Size[0] : Size[0]*Size[2]*TMath::DegToRad(); }
  double GetThick() { return IsBox() ? Size[2] : Size[1]-Size[0]; }

  double GetRmin() { return Size[0]; }
  double GetRmax() { return Size[1]; }
  double GetRmean() { return (Size[0]+Size[0])/2.; }
  double GetdPhi() { return Size[2]; }

  bool IsBox() { return Size[3]>0 ? false : true; }
  bool IsTube() { return Size[3]>0 ? true : false; }
  bool IsCartesian() { return Size[3]<0 ? false : true; }
  bool IsCyl() { return Size[3]<0 ? true : false; }

  double GetLightVelocity() { return LightVelocity; }
  void PrintMap( std::ostream &p_out = std::cout );
};

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
inline void GeomMap::PrintMap( std::ostream &p_out ){
  double tmppar[npar];
  GetParam(tmppar);
  for(int i=0;i<npar;i++)  
    p_out<<std::setw(10)<<"   "<<tmppar[i];
  p_out<<std::endl;
}

inline void GeomMap::SetParam( const double *param ){
  for(int i=0;i<3;i++)
    Pos[i]=param[i];
  for(int i=0;i<3;i++)
    Rot[i]=param[i+3];
  for(int i=0;i<4;i++)
    Size[i]=param[i+6];
  if(IsCyl()) SetRTheta(param[0],param[1]);
  if(IsTube()&&Size[2]<90){
    double tmpr=(GetRmin()+GetRmax())/2.; double tmptheta=GetRotZ();
    SetRTheta(tmpr,tmptheta);
    Pos[0]+=param[0];
    Pos[1]+=param[1];
    R=param[0];
    Theta=param[1];
  }
  LightVelocity=param[10];
}

inline void GeomMap::GetParam( double *param ){
  for(int i=0;i<3;i++)
    param[i]=Pos[i];
  for(int i=0;i<3;i++)
    param[i+3]=Rot[i];
  for(int i=0;i<4;i++)
    param[i+6]=Size[i];
  if(IsCyl()||(IsTube()&&Size[2]<90)){ param[0]=R; param[1]=Theta; }
  param[10]=LightVelocity;
}

inline TVector3 GeomMap::GetPos(const double &ctsub){
  TVector3 dis(0,ctsub*LightVelocity,0);
  if(IsTube()) dis.SetXYZ(0,0,ctsub*LightVelocity);
  dis.RotateX(Rot[0]);
  dis.RotateY(Rot[1]);
  dis.RotateZ(Rot[2]);
  return TVector3(Pos)+dis;
}

class BLDCWireMap
{
public:
  BLDCWireMap();
  ~BLDCWireMap() {};

 private:
  int nWire;
  double Z;
  int XY;
  double XY0, dXY, WireLength, TiltAngle, RotationAngle;
  double WirePhi;

 public:
  void SetParam( const int &nw, const double &z, const int &xy, const double &xy0, const double &dxy,
		 const double &wl, const double &tilt, const double &ra );
  void SetParam(double *param);

  void SetNWire( const int &nw ) { nWire = nw; }
  void SetZ(     const double &z ) { Z = z; }
  void SetXY( const int &xy ) { XY = xy; }
  void SetXY0( const double &xy0 ) { XY0 = xy0; }
  void SetdXY( const double &dxy ) { dXY = dxy; }
  void SetWireLength( const double &wl ) { WireLength = wl; }
  void SetTiltAngle( const double &tilt ) { TiltAngle = tilt; }
  void SetRotationAngle( const double &ra ) { RotationAngle = ra; }
  void SetWirePhi( const double &phi ) { WirePhi = phi; }

  int GetNWire()const  { return nWire; }
  double GetZ() const { return Z; }
  int GetXY()   const  { return XY; }
  double GetXY0() const { return XY0; } 
  double GetdXY() const { return dXY; }
  double GetWireLength() const { return WireLength; }
  double GetTiltAngle()const  { return TiltAngle; }
  double GetRotationAngle()const { return RotationAngle; }
  double GetWirePhi() const{ return WirePhi; }
};

class BLDCWireMapMan
{
public:
  static BLDCWireMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~BLDCWireMapMan(void) {}
  void SetFileName( const TString & filename ) { FileName = filename; }
  bool Initialize();
  bool Initialize( const char *file_name );
  bool Initialize( const std::string& file_name );
  
private:
  BLDCWireMapMan( void );

 private:
  TString FileName;

  typedef std::map < unsigned int, BLDCWireMap > BLDCWireMapContainer;
  BLDCWireMapContainer bldcContainer;
  typedef std::map < unsigned int, GeomMap > GeomMapContainer;
  GeomMapContainer geomContainer;

  static const unsigned int KEYMASK = 0x000F;
  static const unsigned int CMASK   = 0x00FF;
  static const unsigned int LMASK   = 0x00FF;
  static const int          CSHIFT  = 4;
  static const int          LSHIFT  = 16;
  static const unsigned int KEYFLAG = 0x0003; 
  inline int KEY(const int &cid, const int &layer) const
  {  return ((((cid)&CMASK)<<CSHIFT) | (((layer)&LMASK)<<LSHIFT) | KEYFLAG ); }
  static const int MAXCHAR = 256;
  bool m_isready;

  
 public:
  bool IsReady( void ) const { return m_isready; }
  GeomMap * GetGMap( const int &cid, const int &layer);
  const BLDCWireMap * GetWireMap( const int &cid, const int &layer) const;

  double CalcWirePosition( const int &cid, const int &layer, const int &wire) const;
  double GetTiltAngle( const int &cid, const int & layer) const ;
  double GetLocalZ( const int &cid, const int &layer) const ;
  bool GetParam( const int &cid, const int &layer,
		 int &nw, double &z, int &xy, double &xy0, double &dxy, 
		 double &wl, double &tilt, double &ra ) const;
  bool GetXY0( const int &cid, const int &layer, double &xy0) const;
  // bool SetXY0( const int &cid, const int &layer,const double &xy0);
  // bool SetParam( const int &cid, const int &layer,
  // 		 const int &nw,const double &z,const int &xy,const double &xy0,const double &dxy, 
  // 		 const double &wl,const double &tilt,const double &ra );
  int  GetNWire( const int &cid, const int &layer ) const;
  bool GetGParam( const int &cid, TVector3 &pos, TVector3 &rot );

  bool GetGParam( const int &cid, const int &seg, double *par);
  bool SetGParam( const int &cid, const int &seg, const double *par);

  TString GetFileName() { return FileName; }

  void PrintMap( const int &id,std::ostream &p_out = std::cout );
  void PrintMapBL(std::ostream &p_out = std::cout );
};
//______________________________________________________________________________
inline BLDCWireMapMan&
BLDCWireMapMan::GetInstance( void )
{
  static BLDCWireMapMan g_instance;
  return g_instance;
}
inline const std::string&
BLDCWireMapMan::ClassName( void )
{
  static std::string g_name("BLDCWireMan");
  return g_name;
}
#endif
