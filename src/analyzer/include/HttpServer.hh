// -*- C++ -*-

#ifndef HTTP_SERVER_HH
#define HTTP_SERVER_HH

class TCanvas;
class THttpServer;
class TMacro;
class TObject;

//______________________________________________________________________________
class HttpServer
{
public:
  static HttpServer& GetInstance( void );
  ~HttpServer( void );

private:
  HttpServer( void );
  HttpServer( const HttpServer& );
  HttpServer& operator =( const HttpServer& );

private:
  THttpServer           *m_server;
  Int_t                  m_port;
  std::vector<TH1*>      m_th1_list;

public:
  void Open( void );
  void Begin( void );
  void Register( TObject *obj, TString subdir="/");
  void Register( TList *list, TList *parent=0 );
  void Register( TMacro *macro );
  void Register( TString dir, TString command );
  void ResetAll( void );
  void SetPort( Int_t port ){ m_port = port; }

};

//______________________________________________________________________________
inline HttpServer&
HttpServer::GetInstance( void )
{
  static HttpServer g_instance;
  return g_instance;
}

#endif
