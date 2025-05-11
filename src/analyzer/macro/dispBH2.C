// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispBH2()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    TCanvas *c3 = (TCanvas*)gROOT->FindObject("c3");
    c->Clear();
    c->Divide(4,2);
    c3->Clear();
    c3->Divide(4,2);
    int icanvas = 0;
    int icanvas3 = 0;
    for( int ud=0; ud<8; ++ud ){
      TH1 *h  = NULL;
      TH1 *hh = NULL;
      // draw ADC
      c->cd( ++icanvas );
      gPad->SetLogy();
      h = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kADC, ud+1) );
      if( h ){
	h->GetXaxis()->SetRangeUser(0,2000);
	h->Draw();
      }
      hh = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kADCwTDC, ud+1) );
      if( hh ){
	hh->GetXaxis()->SetRangeUser(0,2000);
	hh->SetLineColor( kRed );
	hh->Draw("same");
      }
  }
    for( int ud=0; ud<8; ++ud ){
      TH1 *h  = NULL;
      // draw TDC
      c3->cd( ++icanvas3 );
      gPad->SetLogy();
      h = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kTDC, ud+1) );
      if( h ){
	h->GetXaxis()->SetRangeUser(0,2000);
	h->Draw();
      }
    }
      c->Update();
      c->cd(0);
      c3->Update();
      c3->cd(0);
  }

  {
    TCanvas *c2 = (TCanvas*)gROOT->FindObject("c2");
    TCanvas *c4 = (TCanvas*)gROOT->FindObject("c4");
    c2->Clear();
    c4->Clear();
    c2->Divide(4,2);
    c4->Divide(4,2);
    int icanvas2 = 0;
    int icanvas4 = 0;
    for( int ud=8; ud<16; ++ud ){
      TH1 *h  = NULL;
      TH1 *hh = NULL;
      // draw ADC
      c2->cd( ++icanvas2 );
      gPad->SetLogy();
      h = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kADC, ud+1) );
      if( h ){
	h->GetXaxis()->SetRangeUser(0,2000);
	h->Draw();
      }
      hh = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kADCwTDC, ud+1) );
      if( hh ){
	hh->GetXaxis()->SetRangeUser(0,2000);
	hh->SetLineColor( kRed );
	hh->Draw("same");
      }
    }
    for( int ud=8; ud<16; ++ud ){
      TH1 *h  = NULL;
      // draw TDC
      c4->cd( ++icanvas4 );
      gPad->SetLogy();
      h = (TH1*)GHist::get( HistMaker::getUniqueID(kBH2, 0, kTDC, ud+1) );
      if( h ){
	h->GetXaxis()->SetRangeUser(0,2000);
	h->Draw();
      }
    }
      c2->Update();
      c2->cd(0);
      c4->Update();
      c4->cd(0);
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
