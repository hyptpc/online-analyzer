// Updater belongs to the namespace hddaq::gui
using namespace hddaq::gui;

void dispGeAdc_60Co()
{
  // You must write these lines for the thread safe
  // ----------------------------------
  if(Updater::isUpdating()){return;}
  Updater::setUpdating(true);
  // ----------------------------------

  gStyle->SetTitleW(.5);
  gStyle->SetStatH(.45);
  const int n_calib = 2;
  const int n_seg = 32;
  const double energy[n_calib] = { 1.17, 1.33 };
  const int hmin[n_calib] = { 1200, 1350 };
  const int hmax[n_calib] = { 1400, 1550 };
  const double width = 20.;
  
  double slot[n_calib][n_seg]  = {};
  double Mean[n_calib][n_seg]  = {};
  double eMean[n_calib][n_seg] = {};
  double Sigma[n_calib][n_seg]  = {};
  double eSigma[n_calib][n_seg] = {};
  
  TF1 *f0 = new TF1("f0", "gaus");
  TF1 *f1 = new TF1("f1", "gaus(0) + pol1(3)");
  f0->SetLineColor(kRed); f1->SetLineColor(kGreen);
  
  TLatex t_name, t_val;
  t_name.SetNDC(); t_val.SetNDC();
  t_name.SetTextSize(0.07); t_val.SetTextSize(0.1);

  ////////// Ge ADC at 1.17MeV / 1.33MeV
  for(int i=0; i<n_calib; i++){
    for(int ud=0; ud<2; ud++){
      TCanvas *c = (TCanvas*)gROOT->FindObject(Form("c%d", i*2+ud+1));
      c->Clear();
      c->Divide(4,4);
      int base_id = HistMaker::getUniqueID(kGe, 0, kADC, 1);
      for(int j=0; j<n_seg/2; j++){
	int seg = j + ud*n_seg/2;
	slot[i][seg] = seg + 1;
	c->cd(j+1);
	//gPad->SetLogy();
	TH1 *h = (TH1*)GHist::get(base_id +seg)->Clone();
	h->GetXaxis()->SetRangeUser(hmin[i], hmax[i]);
	if(h->GetEntries()==0){
	  Sigma[i][seg] = 0.; eSigma[i][seg] = 0.;
	  continue;
	}
	double center = h->GetBinCenter(h->GetMaximumBin());
	double param[5] = { 2., center, 2., 0., 0. };
	h->Fit("f0", "Q", "", center-width, center+width);
	f0->GetParameters(param);
	Mean[i][seg]   = f0->GetParameter(1);
	eMean[i][seg]  = f0->GetParError(1);
	Sigma[i][seg]  = f0->GetParameter(2);
	eSigma[i][seg] = f0->GetParError(2);
	double min = param[1] - param[2]*20.;
	double max = param[1] + param[2]*20.;
	f1->SetParameters(param);
	f1->SetParLimits(0, 0., 10000000.);
	f1->SetParLimits(1, hmin[i], hmax[i]);
	f1->SetParLimits(2, 0., 100.);
	for(int it=0; it<5; it++){
	  if(param[2]<20.){
	    h->Fit("f1", "Q", "", min, max);
	    f1->GetParameters(param);
	    Mean[i][seg]  = f1->GetParameter(1);
	    eMean[i][seg] = f1->GetParError(1);
	    Sigma[i][seg]  = f1->GetParameter(2);
	    eSigma[i][seg] = f1->GetParError(2);
	    min = param[1] - param[2]*20.;
	    max = param[1] + param[2]*20.;
	  }
	}
	h->GetXaxis()->SetRangeUser(param[1]-width*1.4, param[1]+width*2.6);
	h->SetTitle(Form("Ge ADC %d @%.2fMeV", seg+1, energy[i]));
	t_name.DrawLatex(0.55, 0.56, "FWHM [ch]");
	t_val.DrawLatex(0.54, 0.44, Form("%.2lf#pm %.2lf",
					 Sigma[i][seg]*2.35482, eSigma[i][seg]*2.35482));
      }
      c->Update();
    }
  }
  ////////// Calc.
  double Gain[n_seg];  // [keV/ch]
  double eGain[n_seg];
  double FWHM[n_seg];  // [keV] @1.33MeV
  double eFWHM[n_seg];
  std::cout<<std::endl<<"Slot\tGain[keV/ch]\tSigma[ch]\tFWHM[keV]"<<std::endl;
  for(int i=0; i<n_seg; i++){
    Gain[i]  = ( 1332.5 - 1173.2 ) / ( Mean[1][i] - Mean[0][i] );
    eGain[i] = ( 1332.5 - 1173.2 ) / pow( Mean[1][i] - Mean[0][i], 2. )
      * sqrt( pow(eMean[0][i],2.) + pow(eMean[1][i],2.) );
    FWHM[i]  = Gain[i] * Sigma[1][i] * 2.35482;
    eFWHM[i] = sqrt( pow(Gain[i]*eSigma[1][i],2.) + pow(eGain[i]*Sigma[1][i],2.) ) * 2.35482;
    if(Sigma[1][i]<=0){
      Gain[i]  = 0.; eGain[i] = 0.;
      FWHM[i]  = 0.; eFWHM[i] = 0.;
    }
    std::cout<<Form("%02d\t%.6f\t%.6f\t%.6f +/- %.6f",
		    i+1, Gain[i], Sigma[1][i], FWHM[i], eFWHM[i])<<std::endl;
  }
  std::cout<<std::endl;
  ////////// Ge ADC FWHM % Slot
  {
    TCanvas *c = (TCanvas*)gROOT->FindObject("c5");
    c->Clear();
    c->SetGrid();
    gStyle->SetEndErrorSize(0);
    TGraph *gr = new TGraphErrors(n_seg, slot[1], FWHM, 0, eFWHM );
    gr->SetTitle("Ge ADC @ 1330keV (60Co)");
    gr->GetXaxis()->SetTitle("slot number");
    gr->GetXaxis()->SetRangeUser(0, n_seg+1);
    gr->GetYaxis()->SetTitle("FWHM [keV]  ");
    gr->GetYaxis()->SetTitleOffset(1.2);
    gr->GetYaxis()->SetRangeUser(0, 15);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(3);
    gr->SetMarkerColor(kGreen);
    gr->SetLineWidth(4);
    gr->Draw("AP");
    c->Update();
  }

  // You must write these lines for the thread safe
  // ----------------------------------
  Updater::setUpdating(false);
  // ----------------------------------
}
