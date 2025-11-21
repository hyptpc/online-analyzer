// -*- C++ -*-

#ifndef MACRO_BUILDER_HH
#define MACRO_BUILDER_HH

#include "HistMaker.hh"

class TCanvas;
class TObject;
class TString;

namespace analyzer
{

//______________________________________________________________________________
namespace macro
{
TObject* Get( TString name );
}

//______________________________________________________________________________
// for HttpServer
namespace http
{
const std::vector<Double_t>& GetEnvelopeXMean();
const std::vector<Double_t>& GetEnvelopeXfitMean();
const std::vector<Double_t>& GetEnvelopeXRMS();
const std::vector<Double_t>& GetEnvelopeXfitSigma();
const std::vector<Double_t>& GetEnvelopeYMean();
const std::vector<Double_t>& GetEnvelopeYfitMean();
const std::vector<Double_t>& GetEnvelopeYRMS();
const std::vector<Double_t>& GetEnvelopeYfitSigma();

TCanvas* Correlation( void );
TCanvas* CorrelationFHT( void );
TCanvas* HitPattern( void );
TCanvas* HitPatternBeam( void );
TCanvas* HitPatternScat( void );
TCanvas* TriggerFlagU( void );
TCanvas* TriggerFlagD( void );
TCanvas* TriggerFlagHitPat( void );
TCanvas* BH1ADC( void );
TCanvas* BH1TDC( void );
TCanvas* BFT( void );
TCanvas* BH2ADC( void );
TCanvas* BH2TDC( void );
TCanvas* HTOFADCU( void );
TCanvas* HTOFADCD( void );
TCanvas* HTOFTDCU( void );
TCanvas* HTOFTDCD( void );
TCanvas* BAC( void );
TCanvas* ACs( void );
TCanvas* BeamHitMulti(void);
TCanvas* SCHTDCU( void );
TCanvas* SCHTDCD( void );
TCanvas* SCHTOTU( void );
TCanvas* SCHTOTD( void );
TCanvas* SCHHitMulti( void );
TCanvas* TOFADCU( void );
TCanvas* TOFADCD( void );
TCanvas* TOFTDCU( void );
TCanvas* TOFTDCD( void );
TCanvas* BVHTDC( void );
TCanvas* LACTDC( void );
TCanvas* WCADCU( void );
TCanvas* WCADCD( void );
TCanvas* WCADCSUM( void );
TCanvas* WCTDCU( void );
TCanvas* WCTDCD( void );
TCanvas* WCTDCSUM( void );
TCanvas* ScatHitMulti(void);
TCanvas* MsTTDC( void );
TCanvas* MsTHitPat( void );
TCanvas* BC3TDCTOT( void );
TCanvas* BC3HitMulti( void );
TCanvas* BC4TDCTOT( void );
TCanvas* BC4HitMulti( void );
TCanvas* SDC1TDCTOT( void );
TCanvas* SDC1HitMulti( void );
TCanvas* SDC2TDCTOT( void );
TCanvas* SDC2HitMulti( void );
TCanvas* SDC3TDCTOT( void );
TCanvas* SDC3HitMulti( void );
TCanvas* SDC4TDCTOT( void );
TCanvas* SDC4HitMulti( void );
TCanvas* T0( void );
TCanvas* BcOutEfficiency( void );
TCanvas* SdcInOutEfficiency( void );
TCanvas* BFTSCHTOT( void );
TCanvas* BFTSFTSCHTOT( void );
TCanvas* DAQ( void );
void     UpdateTOTPeakFitting( void );

// BeamProfile
TCanvas* BeamProfileX( ParticleType p );
TCanvas* BeamProfileY( ParticleType p );
TCanvas* BeamProfileXY( ParticleType p );
TCanvas* BeamProfileFF( void );
TCanvas* BeamEnvelope( ParticleType p );

// E42
TCanvas* TPC( void );
TCanvas* TPC2D( void );
TCanvas* TPC3D( void );
TCanvas* TPCADCPAD( void );
TCanvas* TPCTDCPAD( void );
TCanvas* TPCHTOFPAD( void );
TCanvas* BH2HTOFPAD( void );
TCanvas* TPCFADC( void );

// E03
TCanvas* GeADC( void );
TCanvas* GeCRM( void );
TCanvas* GeTFA( void );
TCanvas* GeRST( void );
TCanvas* GeMultiCRM( void );
TCanvas* GeMultiTFA( void );
TCanvas* GeTFACRM( void );
TCanvas* GeTFAADC( void );
TCanvas* GeRSTADC( void );
TCanvas* BGOTDC( void );
TCanvas* HitPatternE03( void );
void     UpdateTOTPeakFittingE03( void );

// E07
TCanvas* HitPatternE07( void );
TCanvas* SSD1ADCTDC( void );
TCanvas* SSD2ADCTDC( void );
TCanvas* SSD1HitMulti( void );
TCanvas* SSD2HitMulti( void );
TCanvas* SSD1CHitMulti( void );
TCanvas* SSD2CHitMulti( void );
TCanvas* SSD1dETime( void );
TCanvas* SSD2dETime( void );
TCanvas* SSD12Chisqr( void );
TCanvas* SSDEfficiency( void );
TCanvas* EMC( void );
TCanvas* FBHTDC( void );
TCanvas* FBHTOT( void );
TCanvas* FBHHitMulti( void );

// E40
TCanvas* FHTTOT( void );
TCanvas* SAC( void );
TCanvas* LC( void );
TCanvas* TOF_HT( void );
TCanvas* FHT1TDC( void );
TCanvas* FHT1TOT( void );
TCanvas* FHT1HitMulti( void );
TCanvas* FHT2TDC( void );
TCanvas* FHT2TOT( void );
TCanvas* FHT2HitMulti( void );
TCanvas* SFTTDCTOT( void );
TCanvas* SFTHitMulti( void );
TCanvas* SFT2D( void );
TCanvas* HitPatternE40( void );
TCanvas* SFT( void );
TCanvas* FHT1TDC( void );
TCanvas* FHT1TOT( void );
TCanvas* FHT1HitMulti( void );
TCanvas* FHT2TDC( void );
TCanvas* FHT2TOT( void );
TCanvas* FHT2HitMulti( void );
TCanvas* CFTTDC( void );
TCanvas* CFTTDC2D( void );
TCanvas* CFTTOT( void );
TCanvas* CFTTOT2D( void );
TCanvas* CFTHighGain( void );
TCanvas* CFTHighGain2D( void );
TCanvas* CFTLowGain( void );
TCanvas* CFTLowGain2D( void );
TCanvas* CFTPedestal( void );
TCanvas* CFTPedestal2D( void );
TCanvas* CFTHitPat( void );
TCanvas* CFTMulti( void );
TCanvas* CFTEfficiency( void );
TCanvas* CFTClusterHighGain( void );
TCanvas* CFTClusterHighGain2D( void );
TCanvas* CFTClusterLowGain( void );
TCanvas* CFTClusterLowGain2D( void );
TCanvas* CFTClusterTDC( void );
TCanvas* CFTClusterTDC2D( void );
// TCanvas* BGOFADC( void );
// TCanvas* BGOADC( void );
// TCanvas* BGOTDC( void );
// TCanvas* BGOADCTDC2D( void );
// TCanvas* BGOHitMulti( void );
TCanvas* PiIDTDC( void );
TCanvas* PiIDHighGain( void );
TCanvas* PiIDLowGain( void );
TCanvas* PiIDHitMulti( void );
void     UpdateBcOutEfficiency( void );
void     UpdateSdcInOutEfficiency( void );
void     UpdateCFTEfficiency( void );
void     UpdateSSDEfficiency( void );
void     UpdateT0PeakFitting( void );
void     UpdateTOTPeakFittingE40( void );
void     UpdateBeamProfile( ParticleType p );
}

};

#endif
