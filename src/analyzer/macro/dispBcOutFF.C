// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispBcOutFF( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const int n_hist = 5;
    
  // XY position
  for( int i=0; i<n_hist; ++i ){
    TCanvas *c = (TCanvas*)gROOT->FindObject( Form("c%d", i+1 ) );
    c->Clear();
    c->Divide(2,2);

    int hx_id = HistMaker::getUniqueID(kMisc, 0, kHitPat, i+1);
    c->cd(1);
    TH1 *h_x = (TH1*)GHist::get( hx_id );
    if( h_x ){
      // h_x->GetXaxis()->SetRangeUser(-200,200);
      h_x->Draw();
    }

    int hy_id = HistMaker::getUniqueID(kMisc, 0, kHitPat, i+6);
    c->cd(2);
    TH1 *h_y = (TH1*)GHist::get( hy_id );
    if( h_y ){
      // h_y->GetXaxis()->SetRangeUser(-200,200);
      h_y->Draw();
    }

    int hxy_id = HistMaker::getUniqueID(kMisc, 0, kHitPat, i+11);
    c->cd(3);
    TH1 *h_xy = (TH1*)GHist::get( hxy_id );
    if( h_xy ){
      // h_y->GetXaxis()->SetRangeUser(-200,200);
      h_xy->Draw("colz");
    }

    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
