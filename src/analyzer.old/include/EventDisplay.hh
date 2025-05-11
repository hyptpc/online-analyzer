// -*- C++ -*-

#ifndef EVENT_DISPLAY_HH
#define EVENT_DISPLAY_HH

#include "DetectorID.hh"
#include "LorentzVector.hh"
#include "ThreeVector.hh"

class DCLocalTrack;
class K18TrackD2U;
class KuramaTrack;

class TApplication;
class TBox;
class TBRIK;
class TCanvas;
class TFile;
class TGeometry;
class TGraph;
class TH1;
class TLatex;
class TMarker;
class TMarker3DBox;
class TMixture;
class TNode;
class TPad;
class TPave;
class TPaveLabel;
class TPaveText;
class TPolyLine;
class TPolyLine3D;
class TPolyMarker;
class TPolyMarker3D;
class TROOT;
class TRint;
class TRotMatrix;
class TStyle;
class TTRD1;
class TTRD2;
class TTUBS;
class TView;

//______________________________________________________________________________
class EventDisplay// public TObject
{
public:
  static EventDisplay& GetInstance( void );
  static const std::string& ClassName( void );
  ~EventDisplay( void );

private:
  EventDisplay( void );
  EventDisplay( const EventDisplay& );
  EventDisplay& operator =( const EventDisplay& );

private:
  Bool_t                     m_is_ready;
  TLatex                    *m_text;
  TApplication              *m_theApp;
  TGeometry                 *m_geometry;
  TNode                     *m_node;
  TCanvas                   *m_canvas;
  TCanvas                   *m_canvas_vertex;
  TCanvas                   *m_canvas_hist;
  TH1                       *m_hist_vertex_x;
  TH1                       *m_hist_vertex_y;
  TH1                       *m_hist_p;
  TH1                       *m_hist_m2;
  TNode                     *m_target_node;
  TNode                     *m_kurama_inner_node;
  TNode                     *m_kurama_outer_node;
  TNode                     *m_BH2wall_node;
  std::vector<TNode*>        m_BH2seg_node;
  std::vector<TNode*>        m_BC3x1_node;
  std::vector<TNode*>        m_BC3x2_node;
  std::vector<TNode*>        m_BC3v1_node;
  std::vector<TNode*>        m_BC3v2_node;
  std::vector<TNode*>        m_BC3u1_node;
  std::vector<TNode*>        m_BC3u2_node;
  std::vector<TNode*>        m_BC4u1_node;
  std::vector<TNode*>        m_BC4u2_node;
  std::vector<TNode*>        m_BC4v1_node;
  std::vector<TNode*>        m_BC4v2_node;
  std::vector<TNode*>        m_BC4x1_node;
  std::vector<TNode*>        m_BC4x2_node;
  std::vector<TNode*>        m_SFTu_node;
  std::vector<TNode*>        m_SFTv_node;
  std::vector<TNode*>        m_SFTx_node;
  std::vector<TNode*>        m_SDC1v1_node;
  std::vector<TNode*>        m_SDC1v2_node;
  std::vector<TNode*>        m_SDC1x1_node;
  std::vector<TNode*>        m_SDC1x2_node;
  std::vector<TNode*>        m_SDC1u1_node;
  std::vector<TNode*>        m_SDC1u2_node;
  std::vector<TNode*>        m_SDC2x1_node;
  std::vector<TNode*>        m_SDC2x2_node;
  std::vector<TNode*>        m_SDC2y1_node;
  std::vector<TNode*>        m_SDC2y2_node;
  std::vector<TNode*>        m_SDC3x1_node;
  std::vector<TNode*>        m_SDC3x2_node;
  std::vector<TNode*>        m_SDC3y1_node;
  std::vector<TNode*>        m_SDC3y2_node;
  std::vector<TNode*>        m_SDC4y1_node;
  std::vector<TNode*>        m_SDC4y2_node;
  std::vector<TNode*>        m_SDC4x1_node;
  std::vector<TNode*>        m_SDC4x2_node;
  std::vector<TNode*>        m_SSD1y1_node;
  TNode                     *m_FBHwall_node;
  std::vector<TNode*>        m_FBHseg_node;
  TNode                     *m_SCHwall_node;
  std::vector<TNode*>        m_SCHseg_node;
  TNode                     *m_TOFwall_node;
  std::vector<TNode*>        m_TOFseg_node;
  TPolyMarker3D             *m_init_step_mark;
  std::vector<TPolyLine3D*>  m_BcOutTrack;
  std::vector<TPolyLine3D*>  m_SdcInTrack;
  std::vector<TPolyLine3D*>  m_SdcOutTrack;
  TPolyMarker3D             *m_kurama_step_mark;
  // vertex
  TBox                      *m_TargetXZ_box;
  TBox                      *m_TargetYZ_box;
  TBox                      *m_EmulsionXZ_box;
  TBox                      *m_EmulsionYZ_box;
  std::vector<TBox*>         m_FBHseg_box;
  TH1                       *m_SSD_x_hist;
  TH1                       *m_SSD_y_hist;
  TMarker                   *m_VertexPointXZ;
  TMarker                   *m_VertexPointYZ;
  std::vector<TPolyLine*>    m_BcOutXZ_line;
  std::vector<TPolyLine*>    m_BcOutYZ_line;
  std::vector<TPolyLine*>    m_SdcInXZ_line;
  std::vector<TPolyLine*>    m_SdcInYZ_line;
  TPolyMarker               *m_KuramaMarkVertexX;
  TPolyMarker               *m_KuramaMarkVertexY;
  TPolyLine                 *m_MissMomXZ_line;
  TPolyLine                 *m_MissMomYZ_line;

private:
  static const Int_t NumOfMaxrixElements = 9;
  // for Print method
  enum eK18Tracking { kP3rd, kDelta, kDelta3rd, kK18XY, kK18UV,
		      nK18Tracking };
  enum eKuramaTracking { kNHit, kChiSqr, kP, kPath, kPolarity,
			 kMassSqr, kKuramaXY, kKuramaUV,
			 nKuramaTracking };
  enum eReaction { kTheta, kVertex, kClose, kMissMass, kMissMom,
		   nReaction };

public:
  Bool_t Initialize( void );
  Bool_t IsReady( void ) const { return m_is_ready; }
  Bool_t ConstructEmulsion( void );
  Bool_t ConstructTarget( void );
  Bool_t ConstructBH2( void );
  Bool_t ConstructKURAMA( void );
  Bool_t ConstructCollimator( void );
  Bool_t ConstructBcOut( void );
  Bool_t ConstructSdcIn( void );
  Bool_t ConstructSdcOut( void );
  Bool_t ConstructSSD( void );
  Bool_t ConstructFBH( void );
  Bool_t ConstructSCH( void );
  Bool_t ConstructTOF( void );
  void   DrawInitTrack( Int_t nStep, ThreeVector *StepPoint );
  void   DrawInitTrack( void );
  void   DrawHitWire( Int_t lid, Int_t hit_wire,
		      Bool_t range_check=true, Bool_t tdc_check=true );
  void   DrawText( void );
  void   DrawTrigger( const std::vector<Int_t>& flag );
  void   DrawHitHodoscope( Int_t lid, Int_t seg, Int_t Tu=1, Int_t Td=1 );
  void   DrawBcOutLocalTrack( DCLocalTrack *tp );
  void   DrawBcOutLocalTrack( Double_t x0, Double_t y0, Double_t u0, Double_t v0 );
  void   DrawSdcInLocalTrack( DCLocalTrack *tp );
  void   DrawSdcOutLocalTrack( DCLocalTrack *tp );
  void   DrawSsdHit( Int_t lid, Int_t seg, Double_t de );
  void   DrawVertex( const ThreeVector& vertex );
  void   DrawKuramaTrack( Int_t nStep, ThreeVector *StepPoint, Int_t Polarity );
  void   DrawTarget( void );
  void   DrawMissingMomentum( const ThreeVector& mom,
			      const ThreeVector& pos );
  void   DrawMomentum( Double_t momentum );
  void   DrawMassSquare( Double_t mass_square );
  void   EndOfEvent( void );
  void   PrintHit( TString name, Double_t y, Bool_t hit_flag );
  void   PrintK18( K18TrackD2U *tp );
  void   PrintKurama( KuramaTrack *tp, Double_t m2 );
  void   PrintReaction( Double_t theta, const ThreeVector& vertex, Double_t close, const LorentzVector& miss );
  void   ResetVisibility( void );
  void   CalcRotMatrix( Double_t TA, Double_t RA1, Double_t RA2, Double_t *rotMat );
  void   GetCommand( void ) const;
  void   Run( Bool_t flag=kTRUE );
  void   SetStyle( void );
  void   Update( void );

private:
  void   ResetVisibility( TNode *& node, Color_t c=kBlack );
  void   ResetVisibility( std::vector<TNode*>& node, Color_t c=kBlack );

  //  ClassDef(EventDisplay,0);
};

//______________________________________________________________________________
inline EventDisplay&
EventDisplay::GetInstance( void )
{
  static EventDisplay g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
EventDisplay::ClassName( void )
{
  static std::string g_name("EventDisplay");
  return g_name;
}

#endif
