// -*- C++ -*-

#ifndef SSD_ANALYZER_HH
#define SSD_ANALYZER_HH

#include <vector>
#include <Rtypes.h>

//______________________________________________________________________________
class SsdAnalyzer
{
public:
  SsdAnalyzer( void );
  ~SsdAnalyzer( void );

private:
  SsdAnalyzer( const SsdAnalyzer& );
  SsdAnalyzer& operator =( const SsdAnalyzer& );

private:
  // sampling trigger timing
  std::vector<Double_t> m_ssdt;
  std::vector<Double_t> m_ssdct;
  // [layer][segment]
  std::vector< std::vector<Double_t> > m_adc;
  std::vector< std::vector<Double_t> > m_tdc;
  std::vector< std::vector<Double_t> > m_de;
  std::vector< std::vector<Double_t> > m_time;
  std::vector< std::vector<Double_t> > m_chisqr;

public:
  Bool_t   Calculate( void );
  Double_t GetSSDT( Int_t seg ) const { return m_ssdt[seg]; }
  Double_t GetSSDCT( Int_t seg ) const { return m_ssdct[seg]; }
  Double_t GetAdc( Int_t layer, Int_t seg ) const { return m_adc[layer][seg]; }
  Double_t GetTdc( Int_t layer, Int_t seg ) const { return m_tdc[layer][seg]; }
  Double_t GetDe( Int_t layer, Int_t seg ) const { return m_de[layer][seg]; }
  Double_t GetTime( Int_t layer, Int_t seg ) const { return m_time[layer][seg]; }
  Double_t GetChisqr( Int_t layer, Int_t seg ) const { return m_chisqr[layer][seg]; }

};

#endif
