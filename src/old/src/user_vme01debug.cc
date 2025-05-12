// Author: Tomonori Takahashi

#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include <TGFileBrowser.h>
#include <TH1.h>
#include <TString.h>

#include "user_analyzer.hh"
#include "Controller.hh"
#include "ConfMan.hh"
#include "MacroBuilder.hh"

#include <DAQNode.hh>
#include <filesystem_util.hh>
#include <Unpacker.hh>
#include <UnpackerManager.hh>

namespace analyzer
{
  using namespace hddaq::unpacker;
  using namespace hddaq;

  namespace{
    const int n_v792    = 5;
    const int n_ch_v792 = 32;
    TH1F *h792[n_v792][n_ch_v792];
  }

//____________________________________________________________________________
int
process_begin(const std::vector<std::string>& argv)
{
  ConfMan& gConfMan = ConfMan::GetInstance();
  gConfMan.Initialize(argv);
  //  GUnpacker::get_instance().set_decode_mode(false);

  hddaq::gui::Controller& gCon = hddaq::gui::Controller::getInstance();
  TGFileBrowser *tab_hist  = gCon.makeFileBrowser("Hist");
  TGFileBrowser *tab_macro = gCon.makeFileBrowser("Macro");

  // Add macros to the Macro tab
  tab_macro->Add(macro::Get("clear_all_canvas"));
  tab_macro->Add(macro::Get("clear_canvas"));
  tab_macro->Add(macro::Get("dispVme01"));

  // v792 ________________________________________
  for(int i = 0; i<n_v792; ++i){
    for(int j = 0; j<n_ch_v792; ++j){

      h792[i][j] = new TH1F(Form("h%d-%d",i,j), Form("ADC-%d-%d", i, j),
			    4096, 0, 4096);

      tab_hist->Add(h792[i][j]);
    }// for(j)
  }// for(i)

  return 0;
}

//____________________________________________________________________________
int
process_end()
{
  return 0;
}

//____________________________________________________________________________
int
process_event()
{
  static UnpackerManager& gUnpacker = GUnpacker::get_instance();

  {
    static const int k_device = gUnpacker.get_device_id("v792");
    static const int k_adc    = gUnpacker.get_data_id("v792", "adc");

    for(int mod = 0; mod<n_v792; ++mod){
      for(int ch = 0; ch<n_ch_v792; ++ch){
	int nhit = gUnpacker.get_entries( k_device, mod, ch, 0, k_adc );
	if( nhit>0){
	  int adc = gUnpacker.get( k_device, mod, ch, 0, k_adc );
	  h792[mod][ch]->Fill(adc);
	}
      }// for(v792)
    }// for(ch)
  }

  return 0;
}

}
