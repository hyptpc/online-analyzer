#include "HistMaker.hh"

// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

//_____________________________________________________________________
void
dispCFTTDC( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  std::vector<Int_t> id = { kTDC };
  for( Int_t i=0, n=id.size(); i<n; ++i ){
    TCanvas *c = (TCanvas*)gROOT->FindObject( Form("c%d", i+1) );
    c->Clear();
    c->Divide( 4, 2 );
    for( Int_t l=0; l<NumOfLayersCFT; ++l ){
      c->cd(l+1);
      TH1 *h = GHist::get(HistMaker::getUniqueID(kCFT, kCluster, id[i], l+1));
      if(!h) continue;
      h->Draw("colz");
    }
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
