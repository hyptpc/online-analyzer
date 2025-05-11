/**
 *  file: DCHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <std_ostream.hh>

#include "DebugCounter.hh"

//typedef std::vector<bool>   BoolVec;
typedef std::deque<bool>    BoolVec;
typedef std::vector<int>    IntVec;
typedef std::vector<double> DoubleVec;

class DCLTrackHit;

//______________________________________________________________________________
class DCHit
{
public:
  DCHit( void );
  DCHit( int layer );
  DCHit( int layer, double wire );
  DCHit( int cid, int layer, double wire );
  ~DCHit( void );

private:
  DCHit( const DCHit& );
  DCHit& operator =( const DCHit& );

private:
  int       m_cid;
  int       m_layer;
  double    m_wire;
  IntVec    m_tdc;
  IntVec    m_adc;
  IntVec    m_trailing;
  DoubleVec m_dt;
  DoubleVec m_dl;
  DoubleVec m_trailing_time;

  double  m_wpos;
  double  m_angle;
  BoolVec m_belong_track;
  BoolVec m_dl_range;

  ///// for MWPC
  int    m_cluster_size;
  bool   m_mwpc_flag;
  double m_mwpc_wire;
  double m_mwpc_wpos;

  ///// for SSD
  bool      m_is_ssd;
  bool      m_zero_suppressed;
  bool      m_time_corrected;
  bool      m_good_waveform;
  int       m_pedestal;
  int       m_peak_height;
  int       m_peak_position;
  double    m_deviation;
  double    m_amplitude;
  double    m_peak_time;
  double    m_adc_sum;
  double    m_de;
  double    m_rms;
  double    m_chisqr;
  DoubleVec m_time;
  DoubleVec m_waveform;
  bool      m_belong_kaon;
  ///// for TOF
  double    m_z;

  mutable std::vector <DCLTrackHit *> m_register_container;

public:
  bool CalcDCObservables( void );

  void SetDetId( int detid )              { m_cid = detid; }
  void SetLayer( int layer )              { m_layer = layer; }
  void SetWire( double wire )             { m_wire  = wire; }
  void SetTdcVal( int tdc );
  void SetAdcVal( int adc );
  void SetTdcTrailing(int tdc)            { m_trailing.push_back(tdc); }
  void SetDriftTime( double dt )          { m_dt.push_back(dt); }
  void ClearDriftTime( void )             { m_dt.clear(); }
  void SetDriftLength( double dl )        { m_dl.push_back(dl); }
  void SetDriftLength( int ith, double dl ) { m_dl[ith] = dl; }
  void ClearDriftLength( void )           { m_dl.clear();}
  void SetTiltAngle( double angleDegree ) { m_angle = angleDegree; }
  void SetTrailingTime( double t )        { m_trailing_time.push_back(t); }

  void SetClusterSize( int size )          { m_cluster_size   = size; }
  void SetMWPCFlag( bool flag )            { m_mwpc_flag = flag; }
  void SetMeanWire( double mwire )         { m_mwpc_wire    = mwire; }
  void SetMeanWirePosition( double mwpos ) { m_mwpc_wpos    = mwpos; }
  void SetWirePosition( double wpos )      { m_wpos     = wpos; }

  ///// for SSD
  void SetSsdFlag( bool flag=true )          { m_is_ssd        = flag;          }
  void SetGoodWaveForm( bool good=true )     { m_good_waveform = good;          }
  void SetPedestal( int pedestal )           { m_pedestal      = pedestal;      }
  void SetRms( double rms )                  { m_rms           = rms;           }
  void SetAdcSum( double sum )               { m_adc_sum       = sum;           }
  void SetDe( double de )                    { m_de            = de;            }
  void SetDeviation( double deviation )      { m_deviation     = deviation;     }
  void SetTime( DoubleVec time )             { m_time          = time;          }
  void SetWaveform( DoubleVec waveform )     { m_waveform      = waveform;      }
  void SetAmplitude( double amplitude )      { m_amplitude     = amplitude;     }
  void SetPeakTime( double peaktime )        { m_peak_time     = peaktime;      }
  void SetPeakHeight( int height )           { m_peak_height   = height;        }
  void SetPeakPosition( int position )       { m_peak_position = position;      }
  void SetChisquare( double chisqr )         { m_chisqr        = chisqr;        }

  ///// for TOF
  void SetZ( double z ) { m_z = z; }

  int GetDetId( void ) const { return m_cid; }
  int GetLayer( void ) const { return m_layer; }
  double GetWire( void )  const {
    return int(m_wire);
  }

  int GetTdcSize( void ) const { return m_tdc.size(); }
  int GetAdcSize( void ) const { return m_adc.size(); }
  int GetDriftTimeSize( void ) const { return m_dt.size(); }
  int GetDriftLengthSize( void ) const { return m_dl.size(); }
  int GetTdcVal( int nh=0 ) const { return m_tdc[nh]; }
  int GetAdcVal( int nh=0 ) const { return m_adc[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return m_trailing[nh]; }
  int GetTdcTrailingSize( void ) const { return m_trailing.size(); }

  double GetResolution( void ) const;

  double GetDriftTime( int nh=0 ) const { return m_dt[nh]; }
  double GetDriftLength( int nh=0 ) const { return m_dl[nh]; }
  double GetTrailingTime( int nh=0 ) const { return m_trailing_time[nh]; }

  double GetTiltAngle( void ) const { return m_angle; }
  double GetWirePosition( void ) const {
    return m_wpos;
    // if( m_mwpc_flag ) return m_mwpc_wpos;
    // else return m_wpos;
  }

  int GetClusterSize( void ) const { return m_cluster_size; }
  double GetMeamWire( void ) const { return m_mwpc_wire; }
  double GetMeamWirePosition( void ) const { return m_mwpc_wpos; }

  ///// for SSD
  bool      IsSsd( void )              const { return m_is_ssd;          }
  bool      IsTimeCorrected( void )    const { return m_time_corrected;  }
  bool      IsGoodWaveForm( void )     const { return m_good_waveform;   }
  bool      IsZeroSuppressed( void )   const { return m_zero_suppressed; }
  int       GetPedestal( void )        const { return m_pedestal;        }
  DoubleVec GetTime( void )            const { return m_time;            }
  DoubleVec GetWaveform( void )        const { return m_waveform;        }
  double    GetAmplitude( void )       const { return m_amplitude;       }
  double    GetDeviation( void )       const { return m_deviation;       }
  double    GetAdcSum( void )          const { return m_adc_sum;         }
  double    GetDe( void )              const { return m_de;              }
  double    GetPeakTime( void )        const { return m_peak_time;       }
  double    GetRms( void )             const { return m_rms;             }
  double    GetAdcPeakHeight( void )   const { return m_peak_height;     }
  double    GetAdcPeakPosition( void ) const { return m_peak_position;   }
  double    GetChisquare( void )       const { return m_chisqr;          }
  bool      DoTimeCorrection( double offset );
  void      JoinKaonTrack( void ) { m_belong_kaon = true; }
  void      QuitKaonTrack( void ) { m_belong_kaon = false; }
  bool      BelongToKaonTrack( void ) const { return m_belong_kaon; }

  ///// for TOF
  double GetZ( void ) const { return m_z; }

  void JoinTrack( int nh=0 ) { m_belong_track[nh] = true; }
  void QuitTrack( int nh=0 ) { m_belong_track[nh] = false; }
  bool BelongToTrack( int nh=0 ) const { return m_belong_track[nh]; }
  bool IsWithinRange( int nh=0 ) const { return m_dl_range[nh]; }

  void RegisterHits( DCLTrackHit *hit ) const
  { m_register_container.push_back(hit); }

  bool ReCalcDC( bool applyRecursively=false ) { return CalcDCObservables(); }

  void Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;

protected:
  void ClearRegisteredHits( void );
};

//_____________________________________________________________________
inline std::ostream&
operator <<( std::ostream& ost, const DCHit& hit )
{
  hit.Print( "", ost );
  return ost;
}

#endif
