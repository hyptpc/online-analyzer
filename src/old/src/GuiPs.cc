// -*- C++ -*-

#include "GuiPs.hh"

#include <iostream>

#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGTextEntry.h>
#include <TRootBrowser.h>

#include "Controller.hh"
#include "PsMaker.hh"

ClassImp( hddaq::gui::GuiPs )

namespace hddaq
{
namespace gui
{

//_____________________________________________________________________________
GuiPs::GuiPs( void )
  : m_frame(0),
    m_command(k_nCommand),
    m_optButton(),
    m_optList(),
    m_devButton(),
    m_devList(),
    m_textFilename(0),
    m_filename("")
{
}

//_____________________________________________________________________________
GuiPs::~GuiPs( void )
{
}

//_____________________________________________________________________________
TString
GuiPs::getFilename( void )
{
  GuiPs& g = GuiPs::getInstance();
  return g.m_filename;
}

//_____________________________________________________________________________
GuiPs&
GuiPs::getInstance( void )
{
  static GuiPs g_instance;
  return g_instance;
}

//_____________________________________________________________________________
void
GuiPs::initialize( const std::vector<TString>& optList,
                   const std::vector<TString>& devList )
{
  m_optList = optList;
  m_optButton.resize(m_optList.size());

  m_devList = devList;
  m_devButton.resize(m_devList.size());

  hddaq::gui::Controller& gcon = hddaq::gui::Controller::getInstance();
  TRootBrowser* browser = gcon.getBrowser();

  browser->StartEmbedding(0);
  m_frame = new TGMainFrame(gClient->GetRoot());
  browser->StopEmbedding("Ps");

  // check button of "Option"
  TGGroupFrame* optFrame = new TGGroupFrame(m_frame, "Option");
  optFrame->SetTitlePos(TGGroupFrame::kCenter);
  for( Int_t i=0, n=m_optList.size(); i<n; ++i ){
    TGCheckButton* c = new TGCheckButton( optFrame, m_optList[i] );
    optFrame->AddFrame(c, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    m_optButton[i] = c;
    // if( i==k_OptAll )
    //   c->Connect("Toggled(Bool_t)", "hddaq::gui::GuiPs", this,
    //              "toggleAllOption()");
  }
  m_frame->AddFrame(optFrame, new TGLayoutHints(kLHintsExpandX));

  // check button of "Device"
  TGGroupFrame* devFrame = new TGGroupFrame(m_frame, "Device");
  devFrame->SetTitlePos(TGGroupFrame::kCenter);
  for( Int_t i=0, n=m_devList.size(); i<n; ++i ){
    TGCheckButton* c
      = new TGCheckButton( devFrame, m_devList[i] );
    devFrame->AddFrame(c, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    m_devButton[i] = c;
    //       if (i==k_DevAll)
    //	 c->Connect("Toggled(Bool_t)", "hddaq::gui::GuiPs", this,
    //		    "toggleAllDevice()");
  }
  m_frame->AddFrame(devFrame, new TGLayoutHints(kLHintsExpandX));

  // text entry of " filename "
  TGHorizontalFrame* hf = new TGHorizontalFrame(m_frame);
  m_frame->AddFrame(hf, new TGLayoutHints(kLHintsLeft));

  m_textFilename = new TGTextEntry(hf);
  hf->AddFrame(m_textFilename,
	       new TGLayoutHints(kLHintsLeft));//, 2, 2, 2 ,2));
  m_textFilename->SetToolTipText("please input filename and press return");
  m_textFilename->Resize(100, m_textFilename->GetDefaultHeight());
  m_textFilename->Connect("ReturnPressed()", "hddaq::gui::GuiPs", this,
			  "setFilename()");
  m_textFilename->SetText( m_filename );

  // text button of "Save"
  TGHorizontalFrame* hframe = new TGHorizontalFrame(m_frame);

  m_command[k_Save] = new TGTextButton(hframe, "Save");
  m_command[k_Save]->Connect("Clicked()", "hddaq::gui::GuiPs", this, "save()");
  hframe->AddFrame(m_command[k_Save], new TGLayoutHints(kLHintsLeft));

  m_frame->AddFrame(hframe, new TGLayoutHints(kLHintsLeft));

  // text button of "Print"
  m_command[k_Print] = new TGTextButton(hframe, "Print");
  m_command[k_Print]->Connect("Clicked()", "hddaq::gui::GuiPs", this, "print()");
  hframe->AddFrame(m_command[k_Print], new TGLayoutHints(kLHintsLeft));

  m_frame->AddFrame(hframe, new TGLayoutHints(kLHintsLeft));

  m_frame->SetCleanup(kDeepCleanup);
  m_frame->MapSubwindows();
  m_frame->Resize(m_frame->GetDefaultSize());
  m_frame->MapWindow();

  return;
}

//_____________________________________________________________________________
// Bool_t
// GuiPs::isAutoMode( void )
// {
//   return isOptOn( PsMaker::kAutoSaveAtRunChange );
// }

//_____________________________________________________________________________
bool
GuiPs::isDevOn(int i)
{
  GuiPs& g = GuiPs::getInstance();
  return g.m_devButton[i]->IsOn();
}

//_____________________________________________________________________________
bool
GuiPs::isOptOn(int i)
{
  GuiPs& g = GuiPs::getInstance();
  return g.m_optButton[i]->IsOn();
}

//_____________________________________________________________________________
void
GuiPs::print()
{
  std::cout << "#D GuiPs::print() this function is now no longer supported."
            << "   Nothing was done."
            << std::endl;
  // std::cout << "\n#D GuiPs::print()" << std::endl;
  // TString command = "lpr " + m_filename;
  // std::cout << "#D " << command << std::endl;
  return;
}

//_____________________________________________________________________________
void
GuiPs::save()
{
  std::cout << "\n#D GuiPs::save()" << std::endl;
  const std::vector<TGCheckButton*>& oB = m_optButton;
  const std::vector<TGCheckButton*>& dB = m_devButton;

  for (int i = 0; i<(int)m_optList.size(); ++i)
    {
      std::cout << "option = " << m_optList[i] << " : "
		<< oB[i]->IsOn()
		<< "\n";
    }
  std::cout << "\n";
  for (int i = 0; i<(int)m_devList.size(); ++i)
    {
      std::cout << "device = " << m_devList[i] << " : "
		<< dB[i]->IsOn()
		<< "\n";
    }
  std::cout.flush();

  PsMaker::getInstance().makePs();

  return;
}

//_____________________________________________________________________________
void
GuiPs::setFilename()
{
  m_filename = m_textFilename->GetText();
  std::cout << "#D GuiPs::setFilename()"
	    << " file = " << m_filename << std::endl;
  return;
}

//_____________________________________________________________________________
void
GuiPs::setFilename(const char* filename)
{
  m_filename = filename;
  if( m_textFilename )
    m_textFilename->SetText( m_filename );
  return;
}

//_____________________________________________________________________________
void
GuiPs::toggleAllDevice()
{
  // if (m_devButton[k_DevAll]->IsOn())
  //   for (int i=0; i<k_nDevice; ++i)
  //     m_devButton[i]->SetOn();
  // else
  //   for (int i=0; i<k_nDevice; ++i)
  //     m_devButton[i]->SetOn(false);

  return;
}

//_____________________________________________________________________________
void
GuiPs::toggleAllOption()
{
//   if (m_optButton[k_OptAll]->IsOn())
//     {
//       for (int i=0; i<k_nOption; ++i)
// 	{
// 	  if (i!=k_DEBUG)
// 	    m_optButton[i]->SetOn();
// 	}
//     }
//   else
//     for (int i=0; i<k_nOption; ++i)
//       m_optButton[i]->SetOn(false);

  return;
}

}
}
