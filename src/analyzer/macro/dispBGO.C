// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispBGO()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  int n_seg = 12;

  // draw ADC  1-12
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->Clear();
    c->Divide(4,3);
    int fadc_id     = HistMaker::getUniqueID(kBGO, 0, kFADC);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      gPad->SetLogy();
      TH1 *h = (TH1*)GHist::get( fadc_id + i );
      if( !h ) continue;
      h->Draw();
    }
    c->Update();
  }

  // draw ADC 13-24
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
    c->Clear();
    c->Divide(4,3);
    int fadc_id     = HistMaker::getUniqueID(kBGO, 0, kFADC, 12);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      gPad->SetLogy();
      TH1 *h = (TH1*)GHist::get( fadc_id + i );
      if( !h ) continue;
      h->Draw();
    }
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
