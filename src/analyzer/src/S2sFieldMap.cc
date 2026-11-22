// -*- C++ -*-

#include "S2sFieldMap.hh"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <std_ostream.hh>

//#include "ConfMan.hh"
#include "FuncName.hh"

#define DebugDisp 0

//namespace
//{
//const auto& gConf = ConfMan::GetInstance();
//const auto& valueNMR  = ConfMan::Get<Double_t>("FLDNMR");
//const auto& valueCalc = ConfMan::Get<Double_t>("FLDCALC");
//}

//_____________________________________________________________________________
S2sFieldMap::S2sFieldMap(const TString& file_name)
  : m_is_ready(false),
    m_file_name(file_name),
    B(),
    Nx(0),
    Ny(0),
    Nz(0),
    X0(),
    Y0(),
    Z0(),
    dX(),
    dY(),
    dZ()
{
}

//_____________________________________________________________________________
S2sFieldMap::S2sFieldMap(const TString& file_name, const Double_t measure, const Double_t calc)
  : m_is_ready(false),
    m_file_name(file_name),
    B(),
    Nx(0),
    Ny(0),
    Nz(0),
    X0(),
    Y0(),
    Z0(),
    dX(),
    dY(),
    dZ(),
    valueMeasure(measure),
    valueCalc(calc)
{
}


//_____________________________________________________________________________
S2sFieldMap::~S2sFieldMap()
{
  ClearField();
}

//_____________________________________________________________________________
Bool_t
S2sFieldMap::Initialize()
{

  const auto t_start = std::chrono::steady_clock::now();

  if(m_is_ready){
    hddaq::cerr << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }

  std::vector<char> buf;
  {
    std::ifstream ifs(m_file_name, std::ios::binary);
    if(!ifs.is_open()){
      hddaq::cerr << FUNC_NAME << " file open fail : "
                  << m_file_name << std::endl;
      return false;
    }
    ifs.seekg(0, std::ios::end);
    const std::streamoff size = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    buf.resize(size + 1);
    ifs.read(buf.data(), size);
    buf[ifs.gcount()] = '\0';
  }

  ClearField();

  char* p = buf.data();
  char* e = nullptr;
  Bool_t ok = true;
  auto next_int = [&](Int_t& v){
    v = std::strtol(p, &e, 10); if(e == p) ok = false; p = e; };
  auto next_double = [&](Double_t& v){
    v = std::strtod(p, &e); if(e == p) ok = false; p = e; };
  next_int(Nx); next_int(Ny); next_int(Nz);
  next_double(X0); next_double(Y0); next_double(Z0);
  next_double(dX); next_double(dY); next_double(dZ);
  if(!ok){
    hddaq::cerr << FUNC_NAME << " invalid format" << std::endl;
    return false;
  }

  if(Nx<0 || Ny<0 || Nz<0){
    hddaq::cerr << FUNC_NAME << " Nx, Ny, Nz must be positive" << std::endl;
    return false;
  }

  B.resize(Nx);
  for(Int_t ix=0; ix<Nx; ++ix){
    B[ix].resize(Ny);
    for(Int_t iy=0; iy<Ny; ++iy){
      B[ix][iy].resize(Nz);
    }
  }

  if(valueCalc==0. || !std::isfinite(valueCalc) ||
     valueMeasure==0.  || !std::isfinite(valueMeasure) ){
    hddaq::cout << FUNC_NAME << " S2sField is zero : "
                << " Calc = " << valueCalc
                << " Measure = " << valueMeasure << std::endl
                << " -> skip reading fieldmap" << std::endl;
    return true;
  }
  const Double_t factor = valueMeasure/valueCalc;

  hddaq::cout << " reading fieldmap " << std::flush;

  Int_t line = 0;
  Double_t v[6];
  while(true){
    Int_t k = 0;
    for(; k<6; ++k){
      v[k] = std::strtod(p, &e);
      if(e == p) break;
      p = e;
    }
    if(k < 6) break; // end of file (or a truncated last line)
    if(line++%1000000==0)
      hddaq::cout << "." << std::flush;
    const Double_t x = v[0], y = v[1], z = v[2];
    Int_t ix = Int_t((x-X0+0.1*dX)/dX);
    Int_t iy = Int_t((y-Y0+0.1*dY)/dY);
    Int_t iz = Int_t((z-Z0+0.1*dZ)/dZ);
    if(ix>=0 && ix<Nx && iy>=0 && iy<Ny && iz>=0 && iz<Nz){
      B[ix][iy][iz].x = v[3]*factor;
      B[ix][iy][iz].y = v[4]*factor;
      B[ix][iy][iz].z = v[5]*factor;
    }
  }

  const Double_t dt = std::chrono::duration<Double_t>(
    std::chrono::steady_clock::now() - t_start).count();
  hddaq::cout << " done (" << line << " points, " << dt << " s)" << std::endl;
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
S2sFieldMap::GetFieldValue(const Double_t pointCM[3],
                              Double_t* BfieldTesla) const
{
  Double_t xt = pointCM[0];
  Double_t yt = pointCM[1];
  Double_t zt = pointCM[2];

  Int_t ix1, ix2, iy1, iy2, iz1, iz2;
  ix1 = Int_t((xt-X0)/dX);
  iy1 = Int_t((yt-Y0)/dY);
  iz1 = Int_t((zt-Z0)/dZ);

  Double_t wx1, wx2, wy1, wy2, wz1, wz2;
  if(ix1<0) { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if(ix1>=Nx-1) { ix1=ix2=Nx-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(X0+dX*ix2-xt)/dX; wx2=1.-wx1; }

  if(iy1<0) { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if(iy1>=Ny-1) { iy1=iy2=Ny-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(Y0+dY*iy2-yt)/dY; wy2=1.-wy1; }

  if(iz1<0) { iz1=iz2=0; wz1=1.; wz2=0.; }
  else if(iz1>=Nz-1) { iz1=iz2=Nz-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(Z0+dZ*iz2-zt)/dZ; wz2=1.-wz1; }

  Double_t bx1 = wx1*wy1*B[ix1][iy1][iz1].x + wx1*wy2*B[ix1][iy2][iz1].x
    + wx2*wy1*B[ix2][iy1][iz1].x + wx2*wy2*B[ix2][iy2][iz1].x;
  Double_t bx2 = wx1*wy1*B[ix1][iy1][iz2].x + wx1*wy2*B[ix1][iy2][iz2].x
    + wx2*wy1*B[ix2][iy1][iz2].x + wx2*wy2*B[ix2][iy2][iz2].x;
  Double_t bx  = wz1*bx1 + wz2*bx2;

  Double_t by1 = wx1*wy1*B[ix1][iy1][iz1].y + wx1*wy2*B[ix1][iy2][iz1].y
    + wx2*wy1*B[ix2][iy1][iz1].y + wx2*wy2*B[ix2][iy2][iz1].y;
  Double_t by2 = wx1*wy1*B[ix1][iy1][iz2].y + wx1*wy2*B[ix1][iy2][iz2].y
    + wx2*wy1*B[ix2][iy1][iz2].y + wx2*wy2*B[ix2][iy2][iz2].y;
  Double_t by  = wz1*by1 + wz2*by2;

  Double_t bz1 = wx1*wy1*B[ix1][iy1][iz1].z + wx1*wy2*B[ix1][iy2][iz1].z
    + wx2*wy1*B[ix2][iy1][iz1].z + wx2*wy2*B[ix2][iy2][iz1].z;
  Double_t bz2 = wx1*wy1*B[ix1][iy1][iz2].z + wx1*wy2*B[ix1][iy2][iz2].z
    + wx2*wy1*B[ix2][iy1][iz2].z + wx2*wy2*B[ix2][iy2][iz2].z;
  Double_t bz  = wz1*bz1 + wz2*bz2;

  //Default
  BfieldTesla[0] = bx;
  BfieldTesla[1] = by;
  BfieldTesla[2] = bz;

  return true;
}

//_____________________________________________________________________________
void
S2sFieldMap::ClearField()
{
  for(Int_t ix=0; ix<Nx; ++ix){
    for(Int_t iy=0; iy<Ny; ++iy){
      B[ix][iy].clear();
    }
    B[ix].clear();
  }
  B.clear();
}
