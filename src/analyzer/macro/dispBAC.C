// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispBAC()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const int n_seg_ac = 4;
  
  // draw BAC ADC & TDC
  TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
  c->Clear();
  c->Divide(2,2);
  int idac[n_seg_ac] = {
    HistMaker::getUniqueID(kBAC,  0, kADC, 1),
    HistMaker::getUniqueID(kBAC,  0, kADC, 2),
    HistMaker::getUniqueID(kBAC,  0, kTDC, 1),
    HistMaker::getUniqueID(kBAC,  0, kTDC, 2)
  };

  for(int i = 0; i<n_seg_ac; ++i){
    c->cd(i+1);
    gPad->SetLogy();
    TH1 *h = (TH1*)GHist::get(idac[i]);
    h->GetXaxis()->SetRangeUser(0,4000);
    h->Draw();
  }

  c->Update();

  c->cd(0);

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
