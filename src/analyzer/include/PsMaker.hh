#ifndef H_PSMAKERIMPL
#define H_PSMAKERIMPL

#include<vector>
#include<string>

class TCanvas;
class TPostScript;

class PsMaker{
  // Option list which are displayed in the Ps tab
  enum OptionList{
    kExpDataSheet, 
    kFixXaxis,
    kLogyADC, kLogyTDC,//, kLogyHitPat, kLogyMulti,
    sizeOptionList
  };
  std::vector<std::string> name_option_;

  // Parameter list controlling the histogram drawing
  enum ParameterList{
    kXdiv, kYdiv, kXrange_min, kXrange_max,
    sizeParameterList
  };

  // The instance to make the ps file
  TCanvas     *cps_;
  TPostScript *ps_;

public:
  virtual ~PsMaker();
  static PsMaker& getInstance();
  
  void getListOfOption(std::vector<std::string>& vec);
  void makePs();

private:
  void drawRunNumber();
  void create(std::string& name);
  void drawOneCanvas(std::vector<int>& id_list, std::vector<int>& par_list,
		     bool flag_xaxis, bool flag_log,
		     const char* optDraw = ""
		     );
  void drawDCEff( void );
  void clearOneCanvas(int npad);  

  PsMaker();
  PsMaker(const PsMaker& object);
  PsMaker& operator=(const PsMaker& object);
};

inline PsMaker& PsMaker::getInstance()
{
  static PsMaker object;
  return object;
}

#endif
