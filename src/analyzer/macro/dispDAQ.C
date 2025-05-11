// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispDAQ( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  //gStyle->SetOptStat(1111110);

  const Int_t n_hist = 6;
  Int_t hist_id[n_hist] = {
    HistMaker::getUniqueID(kDAQ, kEB,       kHitPat,   1),
    HistMaker::getUniqueID(kDAQ, kVME,      kHitPat2D, 1),
    HistMaker::getUniqueID(kDAQ, kCLite,    kHitPat2D, 1),
    HistMaker::getUniqueID(kDAQ, kEASIROC,  kHitPat2D, 1),
    HistMaker::getUniqueID(kDAQ, kCAMAC,    kHitPat2D, 1),
    HistMaker::getUniqueID(kDAQ, kMiscNode, kHitPat2D, 1)
  };

  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->Clear();
    c->Divide(3,2);
    for( Int_t i=0; i<n_hist; ++i ){
      c->cd(i+1);
      gPad->SetGrid();
      TH1 *h = GHist::get( hist_id[i] );
      if( h ) h->Draw("colz");
    }
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
