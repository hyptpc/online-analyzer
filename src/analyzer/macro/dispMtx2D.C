// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispMtx2D( void )
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  const Int_t n_seg_tof = 24;
  const Int_t n_seg_sch = 64;

  // Draw Mtx2D pattern and TOFxSCH
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c1");
    c->cd(0);
    Int_t mtx2d_id  = HistMaker::getUniqueID(kMisc,  kHul2D,      kHitPat2D);
    Int_t tofsch_id = HistMaker::getUniqueID(kMtx3D, kHul2DHitPat,kHitPat2D);
    TH2 *hmtx2d  = (TH2*)GHist::get(mtx2d_id);
    TH2 *htofsch = (TH2*)GHist::get(tofsch_id);

    hmtx2d->Reset();
    hmtx2d->SetLineWidth(2);
    hmtx2d->SetLineColor(2);
    MatrixParamMan& gMatrix = MatrixParamMan::GetInstance();
    for( Int_t i_tof=0; i_tof<n_seg_tof; ++i_tof ){
      for( Int_t i_sch=0; i_sch<n_seg_sch; ++i_sch ){
	bool hul2d_flag = gMatrix.IsAccept( i_tof, i_sch );
	if( !hul2d_flag ) hmtx2d->Fill( i_sch, i_tof );
      } // for(i_sch)
    } // for(i_tof)

    hmtx2d->GetXaxis()->SetRangeUser(0, n_seg_sch);
    hmtx2d->GetYaxis()->SetRangeUser(0, n_seg_tof+1);
    hmtx2d->Draw("box");

    htofsch->GetXaxis()->SetRangeUser(0, n_seg_sch);
    htofsch->GetYaxis()->SetRangeUser(0, n_seg_tof);
    htofsch->Draw("same");

    gPad->SetGridx();
    gPad->SetGridy();

    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
