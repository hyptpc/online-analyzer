// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispGeTdc()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const int NumOfSegGe = 32;
  ////////// Ge TDC
  {
    const int n_type = 4;
    const int tdc_type[n_type] = { kCRM, kTFA, kPUR, kRST };
    // draw TDC
    for(int t=0; t<n_type; t++){
      TCanvas *c = (TCanvas*)gROOT->FindObject(Form("c%d", t+1));
      c->Clear();
      c->Divide(8,4);
      int base_id = HistMaker::getUniqueID(kGe, 0, tdc_type[t], 1);
      for(int i=0; i<NumOfSegGe; i++){
	c->cd(i+1);
	gPad->SetLogy();
	TH1 *h = (TH1*)GHist::get(base_id + i);
	h->Draw();
      }
      c->Update();
    }
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
