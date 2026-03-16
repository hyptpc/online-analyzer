/**
 *  file: UserParamMan.hh
 *  date: 2017.04.10
 *
 */

#ifndef USER_PARAM_MAN_HH
#define USER_PARAM_MAN_HH

#include <TString.h>
#include <string>
#include <map>
#include <vector>

//______________________________________________________________________________
class UserParamMan
{
public:
  static UserParamMan&      GetInstance( void );
  static const std::string& ClassName( void );
  ~UserParamMan( void );

private:
  UserParamMan( void );
  UserParamMan( const UserParamMan&  );
  UserParamMan& operator =( const UserParamMan& );

private:
  typedef std::vector<double>               ParamArray;
  typedef std::map<std::string, ParamArray> ParamMap;
  typedef ParamMap::const_iterator          PIterator;
  bool        m_is_ready;
  std::string m_file_name;
  ParamMap    m_param_map;

public:
  bool    Initialize( void );
  bool    Initialize( const std::string& filename );
  bool    IsReady( void ) const { return m_is_ready; }
  Bool_t   IsInRange(const std::string& key, Double_t val) const;
  int     GetSize( const std::string& key ) const;
  double  GetParameter( const std::string& key, std::size_t i=0 ) const;
  double  Get(const TString& key, int i=0)const
  {return GetParameter(std::string(key.Data()),i);}
  void    SetFileName( std::string& name ) { m_file_name=name; }
  void    Print( const std::string& arg="" ) const;
};

//______________________________________________________________________________
inline UserParamMan&
UserParamMan::GetInstance( void )
{
  static UserParamMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
UserParamMan::ClassName( void )
{
  static std::string g_name("UserParamMan");
  return g_name;
}

#endif
