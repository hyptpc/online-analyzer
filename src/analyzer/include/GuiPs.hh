#ifndef GUIPS_H
#define GUIPS_H

#include <vector>

#include <TGFrame.h>

class TGCompositeFrame;
class TGCheckButton;
class TGTextButton;
class TGTextEntry;

namespace hddaq
{
  namespace gui
  {

class GuiPs
{

public:

  enum e_command
    {
      k_Save,
      k_Print,
      k_nCommand
    };

private:
  TGCompositeFrame*           m_frame;
  std::vector<TGTextButton*>  m_command;
  std::vector<TGCheckButton*> m_optButton;
  std::vector<std::string>    m_optList;
  std::vector<TGCheckButton*> m_devButton;
  std::vector<std::string>    m_devList;
  TGTextEntry*                m_textFilename;
  std::string                 m_filename;

public:
  static GuiPs& getInstance();
  virtual ~GuiPs();

  static std::string getFilename();
  void initialize(const std::vector<std::string>& optList,
		  const std::vector<std::string>& devList);
  static bool isDevOn(int i);
  static bool isOptOn(int i);
  void print();
  void save();
  void setFilename();
  void setFilename(const char* filename);
  void toggleAllDevice();
  void toggleAllOption();

private:
  GuiPs();
  GuiPs(const GuiPs&);
  GuiPs& operator=(const GuiPs&);

  ClassDef(hddaq::gui::GuiPs,0)

    };

  }
}
#endif


