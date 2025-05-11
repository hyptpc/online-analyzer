// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispSCH()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const int n_seg     = 64;

  int sch_id_c1[] = {
    HistMaker::getUniqueID(kSCH, 0, kTDC,    1),
    HistMaker::getUniqueID(kSCH, 0, kADC,    1),
    HistMaker::getUniqueID(kSCH, 0, kTDC,    n_seg +1),
    HistMaker::getUniqueID(kSCH, 0, kADC,    n_seg +1),
    HistMaker::getUniqueID(kSCH, 0, kHitPat, 1),
    HistMaker::getUniqueID(kSCH, 0, kTDC2D,  1),
    HistMaker::getUniqueID(kSCH, 0, kADC2D,  1),
    HistMaker::getUniqueID(kSCH, 0, kMulti,  1)
  };

  // // Draw TDC 01-32
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
  //   c->Clear();
  //   c->Divide(8,4);
  //   for( int i=0; i<n_seg/2; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( sch_id_c1[0]+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TDC 33-64
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
  //   c->Clear();
  //   c->Divide(8,4);
  //   for( int i=0; i<n_seg/2; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( sch_id_c1[0]+i+n_seg/2 );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TOT 01-32
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c3");
  //   c->Clear();
  //   c->Divide(8,4);
  //   for( int i=0; i<n_seg/2; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( sch_id_c1[1]+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TOT 33-64
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c4");
  //   c->Clear();
  //   c->Divide(8,4);
  //   for( int i=0; i<n_seg/2; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( sch_id_c1[1]+i+n_seg/2 );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // draw TDC/TOT
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c5");
    c->Clear();
    c->Divide(3,2);
    for(int i=0; i<6; i++){
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( sch_id_c1[2+i] );
      h->Draw("colz");
    }
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
