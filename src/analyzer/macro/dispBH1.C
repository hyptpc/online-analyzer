// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispBH1()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  int n_seg = 11;

  // draw ADC U
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->Clear();
    c->Divide(4,3);
    int adc_id     = HistMaker::getUniqueID(kBH1, 0, kADC, 1);
    int adcwtdc_id = HistMaker::getUniqueID(kBH1, 0, kADCwTDC, 1);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      gPad->SetLogy();
      TH1 *h = (TH1*)GHist::get( adc_id + i );
      if( !h ) continue;
      h->GetXaxis()->SetRangeUser(0,2000);
      h->Draw();
      TH1 *hh = (TH1*)GHist::get( adcwtdc_id + i );
      if( !hh ) continue;
      hh->GetXaxis()->SetRangeUser(0,2000);
      hh->SetLineColor( kRed );
      hh->Draw("same");
    }
    c->Update();
  }

  // draw ADC D
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
    c->Clear();
    c->Divide(4,3);
    int adc_id     = HistMaker::getUniqueID(kBH1, 0, kADC, 1+n_seg);
    int adcwtdc_id = HistMaker::getUniqueID(kBH1, 0, kADCwTDC, 1+n_seg);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      gPad->SetLogy();
      TH1 *h = (TH1*)GHist::get( adc_id + i );
      if( !h ) continue;
      h->GetXaxis()->SetRangeUser(0,2000);
      h->Draw();
      TH1 *hh = (TH1*)GHist::get( adcwtdc_id + i );
      if( !hh ) continue;
      hh->GetXaxis()->SetRangeUser(0,2000);
      hh->SetLineColor( kRed );
      hh->Draw("same");
    }
    c->Update();
  }

  // draw TDC U
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c3");
    c->Clear();
    c->Divide(4,3);
    int tdc_id = HistMaker::getUniqueID(kBH1, 0, kTDC, 1);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      TH1 *h = (TH1*)GHist::get( tdc_id + i );
      if( !h ) continue;
      h->Draw();
    }
    c->Update();
  }

  // draw TDC D
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c4");
    c->Clear();
    c->Divide(4,3);
    int tdc_id = HistMaker::getUniqueID(kBH1, 0, kTDC, 1+n_seg);
    for( int i=0; i<n_seg; ++i ){
      c->cd(i+1);
      TH1 *h = (TH1*)GHist::get( tdc_id + i);
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
