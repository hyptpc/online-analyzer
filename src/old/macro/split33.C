// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void split33()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  TCanvas *c = (TCanvas*)gROOT->FindObject(gPad->GetName());
  c->Divide(3,3);
  c->cd(1);

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
