// -*- C++ -*-

#include "FieldMan.hh"

#include <iostream>

#include "ConfMan.hh"
#include "S2sFieldMap.hh"

//_____________________________________________________________________________
FieldMan::FieldMan()
  : m_is_ready(false),
    m_file_name(),
    m_shs_map(nullptr)
{
}

//_____________________________________________________________________________
FieldMan::~FieldMan()
{
  delete m_shs_map;
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize(const TString& file_name, Double_t map_scale)
{
  if(m_is_ready){
    std::cerr << "#W FieldMan already initialized" << std::endl;
    return false;
  }
  m_file_name = file_name;
  delete m_shs_map;
  // S2sFieldMap applies factor = measure/calc to the map B
  m_shs_map = new S2sFieldMap(m_file_name, kHall*map_scale, kCalc);
  if(!m_shs_map->Initialize()){
    delete m_shs_map;
    m_shs_map = nullptr;
    return false;
  }
  std::cout << "#D FieldMan SHS map " << m_file_name
            << " factor = " << kHall*map_scale/kCalc << std::endl;
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
ThreeVector
FieldMan::GetField(const ThreeVector& position) const
{
  ThreeVector field(0., 0., 0.);
  if(m_shs_map){
    // the map is in [cm]
    const Double_t p[3] = { position.x()*0.1, position.y()*0.1,
                            position.z()*0.1 };
    Double_t b[3];
    if(m_shs_map->GetFieldValue(p, b))
      field.SetXYZ(b[0], b[1], b[2]);
  }
  return field;
}

//______________________________________________________________________________
// Optional: without a "FLDMAP:" key the beam RK histograms stay empty.
void
ConfMan::initializeFieldMan()
{
  if(name_file_["FLDMAP:"] != ""){
    if(!FieldMan::GetInstance().Initialize(name_file_["FLDMAP:"])){
      std::cout << "#W ConfMan:: FieldMan initialization failed,"
                << " beam RK histograms will stay empty" << std::endl;
    }
  }else{
    std::cout << "#W ConfMan:: no FLDMAP: key in the conf file,"
              << " beam RK histograms will stay empty" << std::endl;
  }
}
