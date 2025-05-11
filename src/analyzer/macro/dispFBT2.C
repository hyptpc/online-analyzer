// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispFBT2( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const int n_seg1     = 48;
  const int n_seg2     = 64;
  const int FBTOffset = 200;
    
  int fbt2_1_U_id[] = {
//    HistMaker::getUniqueID(kFBT2, 0, kTDC,    1),
//    HistMaker::getUniqueID(kFBT2, 0, kADC,    1),
    HistMaker::getUniqueID(kFBT2, 0, kTDC,    n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 0, kADC,    n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 0, kHitPat, 1),
    HistMaker::getUniqueID(kFBT2, 0, kTDC2D,  1),
    HistMaker::getUniqueID(kFBT2, 0, kADC2D,  1),
    HistMaker::getUniqueID(kFBT2, 0, kMulti,  1)
  };

  int fbt2_1_D_id[] = {
//    HistMaker::getUniqueID(kFBT2, 0, kTDC,    1),
//    HistMaker::getUniqueID(kFBT2, 0, kADC,    1),
    HistMaker::getUniqueID(kFBT2, 0, kTDC,    FBTOffset+n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 0, kADC,    FBTOffset+n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 0, kHitPat, FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 0, kTDC2D,  FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 0, kADC2D,  FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 0, kMulti,  FBTOffset+1)
  };


  int fbt2_2_U_id[] = {
//    HistMaker::getUniqueID(kFBT2, 1, kTDC,    1),
//    HistMaker::getUniqueID(kFBT2, 1, kADC,    1),
    HistMaker::getUniqueID(kFBT2, 1, kTDC,    n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 1, kADC,    n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 1, kHitPat, 1),
    HistMaker::getUniqueID(kFBT2, 1, kTDC2D,  1),
    HistMaker::getUniqueID(kFBT2, 1, kADC2D,  1),
    HistMaker::getUniqueID(kFBT2, 1, kMulti,  1)
  };

  int fbt2_2_D_id[] = {
//    HistMaker::getUniqueID(kFBT2, 1, kTDC,    1),
//    HistMaker::getUniqueID(kFBT2, 1, kADC,    1),
    HistMaker::getUniqueID(kFBT2, 1, kTDC,    FBTOffset+n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 1, kADC,    FBTOffset+n_seg2+1),
    HistMaker::getUniqueID(kFBT2, 1, kHitPat, FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 1, kTDC2D,  FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 1, kADC2D,  FBTOffset+1),
    HistMaker::getUniqueID(kFBT2, 1, kMulti,  FBTOffset+1)
  };


  // // Draw TDC UP
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
  //   c->Clear();
  //   c->Divide(4,4);
  //   for( int i=0; i<n_seg; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( fbt2_id_c1[0]+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TDC DOWN
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
  //   c->Clear();
  //   c->Divide(4,4);
  //   for( int i=0; i<n_seg; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( fbt2_id_c1[0]+n_seg+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TOT UP
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c3");
  //   c->Clear();
  //   c->Divide(4,4);
  //   for( int i=0; i<n_seg; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( fbt2_id_c1[1]+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  // // Draw TOT DOWN
  // {
  //   TCanvas *c = (TCanvas*)gROOT->FindObject("c4");
  //   c->Clear();
  //   c->Divide(4,4);
  //   for( int i=0; i<n_seg; ++i ){
  //     c->cd(i+1)->SetGrid();
  //     TH1 *h = (TH1*)GHist::get( fbt2_id_c1[1]+n_seg+i );
  //     h->Draw();
  //   }
  //   c->Update();
  // }

  //
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->Clear();
    c->Divide(3,2);
    for( int i=0; i<6; ++i ){
     if(i<3){
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_1_U_id[i] );
      h->Draw("colz");
    }else{
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_1_D_id[i-3] );
      h->Draw("colz");
    }
    c->Update();
  }
  }
    
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c2");
    c->Clear();
    c->Divide(3,2);
    for( int i=0; i<6; ++i ){
     if(i<3){
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_1_U_id[i+3] );
      h->Draw("colz");
    }else{
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_1_D_id[i] );
      h->Draw("colz");
    }
    c->Update();
  }
  }

  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c3");
    c->Clear();
    c->Divide(3,2);
    for( int i=0; i<6; ++i ){
     if(i<3){
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_2_U_id[i] );
      h->Draw("colz");
    }else{
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_2_D_id[i-3] );
      h->Draw("colz");
    }
    c->Update();
  }
  }
    
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c4");
    c->Clear();
    c->Divide(3,2);
    for( int i=0; i<6; ++i ){
     if(i<3){
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_2_U_id[i+3] );
      h->Draw("colz");
    }else{
      c->cd(i+1);//->SetGrid();
      TH1 *h = (TH1*)GHist::get( fbt2_2_D_id[i] );
      h->Draw("colz");
    }
    c->Update();
  }
  }


//  {
//    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
//    c->Clear();
//    c->Divide(3,2);
//    for( int i=0; i<6; ++i ){
//      c->cd(i+1)->SetGrid();
//      TH1 *h = (TH1*)GHist::get( fbt2_1_id[i] );
//      h->Draw("colz");
//    }
//    c->Update();
//  }
//    
//  {
//    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
//    c->Clear();
//    c->Divide(3,2);
//    for( int i=0; i<6; ++i ){
//      c->cd(i+1)->SetGrid();
//      TH1 *h = (TH1*)GHist::get( fbt2_1_id[i] );
//      h->Draw("colz");
//    }
//    c->Update();
//  }
//    
  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
