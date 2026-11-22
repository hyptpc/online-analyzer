// -*- C++ -*-

#ifndef FIELD_MAN_HH
#define FIELD_MAN_HH

#include <TString.h>

#include "ThreeVector.hh"

class S2sFieldMap;

//_____________________________________________________________________________
// HS (SHS) field map for the beam RK propagation, ported from the offline
// analyzer (FieldMan, SHS part only). Positions are in the HypTPC local
// frame [mm], fields in [T]. The map is scaled by Hall/Calc as offline:
//   B = B_map * (kHall*SHSMAPSCALE) / kCalc
class FieldMan
{
public:
  static FieldMan& GetInstance();
  ~FieldMan();

  // E72 HS field was stable through beamtime (offline hs_field constants)
  static constexpr Double_t kHall = 0.7872;
  static constexpr Double_t kCalc = 0.780;

private:
  FieldMan();
  FieldMan(const FieldMan&);
  FieldMan& operator =(const FieldMan&);

  Bool_t       m_is_ready;
  TString      m_file_name;
  S2sFieldMap* m_shs_map;

public:
  Bool_t      Initialize(const TString& file_name, Double_t map_scale=1.);
  Bool_t      IsReady() const { return m_is_ready; }
  ThreeVector GetField(const ThreeVector& position) const;
};

//_____________________________________________________________________________
inline FieldMan&
FieldMan::GetInstance()
{
  static FieldMan s_instance;
  return s_instance;
}

#endif
