// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispSDC1( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  gStyle->SetOptStat(1111110);
  int n_layer = 6;

  // draw TDC
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->Clear();
    c->Divide(3,2);
    int base_id = HistMaker::getUniqueID(kSDC1, 0, kTDC);
    for( int i=0; i<n_layer; ++i ){
      c->cd(i+1);
      TH1 *h = (TH1*)GHist::get(base_id + i);;
      if( h ) h->Draw();
    }
    c->Update();
  }

  // draw TDC1st
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
    c->Clear();
    c->Divide(3,2);
    int base_id = HistMaker::getUniqueID(kSDC1, 0, kTDC2D);
    for( int i=0; i<n_layer; ++i ){
      c->cd(i+1);
      TH1 *h = (TH1*)GHist::get(base_id + i);;
      if( h ) h->Draw();
    }
    c->Update();
  }

  // draw HitPat
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c3");
    c->Clear();
    c->Divide(3,2);
    int base_id = HistMaker::getUniqueID(kSDC1, 0, kHitPat);
    for( int i=0; i<n_layer; ++i ){
      c->cd(i+1);
      TH1 *h = (TH1*)GHist::get(base_id + i);
      if( h ) h->Draw();
    }
    c->Update();
  }

  // draw Multi
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c4");
    c->Clear();
    c->Divide(3,2);
    int base_id = HistMaker::getUniqueID(kSDC1, 0, kMulti);
    for( int i=0; i<n_layer; ++i ){
      c->cd(i+1);
      TH1 *h_wot = (TH1*)GHist::get(base_id + i);//->Clone();
      TH1 *h_wt  = (TH1*)GHist::get(base_id + i + n_layer);//->Clone();
      if( !h_wot || !h_wt ) continue;
      // h_wot->SetMaximum(h_wt->GetMaximum()*1.1);
      h_wot->Draw();
      h_wt->SetLineColor(2);
      h_wt->Draw("same");
    }
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
