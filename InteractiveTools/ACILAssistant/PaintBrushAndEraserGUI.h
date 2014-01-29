// generated by Fast Light User Interface Designer (fluid) version 1.0107

#ifndef PaintBrushAndEraserGUI_h
#define PaintBrushAndEraserGUI_h
#include <FL/Fl.H>
#include "itkImage.h"
#include <vector>
using namespace std;
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Round_Button.H>
#include <FL/Fl_Spinner.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Choice.H>

class PaintBrushAndEraserGUI {
typedef itk::Image< unsigned short, 3 > LabelMapType;
public:
  PaintBrushAndEraserGUI();
  Fl_Double_Window *paintBrushAndEraserWindow;
  Fl_Round_Button *paintBrushRadio;
  Fl_Round_Button *eraserRadio;
  Fl_Spinner *toolSizeSpinner;
  Fl_Group *eraserGroup;
  Fl_Button *eraseButton;
  Fl_Check_Button *eraseSelectedCheck;
  Fl_Group *toolThresholdsGroup;
  Fl_Check_Button *tool3DCheck;
  Fl_Value_Input *minValueInput;
  Fl_Value_Input *maxValueInput;
  Fl_Choice *chestRegionChoice;
  Fl_Choice *chestTypeChoice;
private:
  static void paintBrushRadio_CB( Fl_Widget* o, void* v );
  void paintBrushRadio_CB_i();
  static void eraserRadio_CB( Fl_Widget* o, void* v );
  void eraserRadio_CB_i();
  static void eraseButton_CB( Fl_Widget* o, void* v );
  void eraseButton_CB_i();
public:
  unsigned char GetPaletteSelection();
  unsigned char GetChestRegion();
  unsigned char GetChestType();
  unsigned int GetToolRadius();
  bool GetPaintBrushSelected();
  bool GetEraserSelected();
  bool GetEraseSelectedSelected();
  void SetLabelMapImage( LabelMapType::Pointer labelMapImage );
  void SetUpdateViewerFunction( void (*f)() );
  void SetPaintedIndices( std::vector< LabelMapType::IndexType >* indices );
private:
  LabelMapType::Pointer m_LabelMapImage;
  void (*m_UpdateViewerFunction)();
  unsigned char m_PaletteSelection;
  unsigned char m_ChestRegion;
  unsigned char m_ChestType;
  std::vector< LabelMapType::IndexType >* m_PaintedIndices;
public:
  short GetToolLowerThreshold();
  short GetToolUpperThreshold();
private:
  static void undefinedRegionMenuItem_CB( Fl_Widget* o, void* v );
  void undefinedRegionMenuItem_CB_i();

  static void ambiguousBronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v );
  void ambiguousBronchiectaticAirwayMenuItem_CB_i();

  static void nonBronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v );
  void nonBronchiectaticAirwayMenuItem_CB_i();

  static void bronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v );
  void bronchiectaticAirwayMenuItem_CB_i();

  static void rightLungMenuItem_CB( Fl_Widget* o, void* v );
  void rightLungMenuItem_CB_i();

  static void leftLungMenuItem_CB( Fl_Widget* o, void* v );
  void leftLungMenuItem_CB_i();

  static void rightUpperLobeMenuItem_CB( Fl_Widget* o, void* v );
  void rightUpperLobeMenuItem_CB_i();

  static void rightMiddleLobeMenuItem_CB( Fl_Widget* o, void* v );
  void rightMiddleLobeMenuItem_CB_i();

  static void rightLowerLobeMenuItem_CB( Fl_Widget* o, void* v );
  void rightLowerLobeMenuItem_CB_i();

  static void leftUpperLobeMenuItem_CB( Fl_Widget* o, void* v );
  void leftUpperLobeMenuItem_CB_i();

  static void leftLowerLobeMenuItem_CB( Fl_Widget* o, void* v );
  void leftLowerLobeMenuItem_CB_i();

  static void leftMenuItem_CB( Fl_Widget* o, void* v );
  void leftMenuItem_CB_i();

  static void rightMenuItem_CB( Fl_Widget* o, void* v );
  void rightMenuItem_CB_i();

  static void abdomen_CB( Fl_Widget* o, void* v );
  void abdomen_CB_i();

  static void aorta_CB( Fl_Widget* o, void* v );
  void aorta_CB_i();

  static void paravertebral_CB( Fl_Widget* o, void* v );
  void paravertebral_CB_i();

  static void muscle_CB( Fl_Widget* o, void* v );
  void muscle_CB_i();

  static void liverMenuItem_CB( Fl_Widget* o, void* v );
  void liverMenuItem_CB_i();

  static void spleenMenuItem_CB( Fl_Widget* o, void* v );
  void spleenMenuItem_CB_i();

  static void undefinedTypeMenuItem_CB( Fl_Widget* o, void* v );
  void undefinedTypeMenuItem_CB_i();

  static void normalParenchymaMenuItem_CB( Fl_Widget* o, void* v );
  void normalParenchymaMenuItem_CB_i();

  static void airwayGeneration3MenuItem_CB( Fl_Widget* o, void* v );
  void airwayGeneration3MenuItem_CB_i();

  static void airwayGeneration4MenuItem_CB( Fl_Widget* o, void* v );
  void airwayGeneration4MenuItem_CB_i();

  static void airwayGeneration5MenuItem_CB( Fl_Widget* o, void* v );
  void airwayGeneration5MenuItem_CB_i();

  static void pectoralisMinorMenuItem_CB( Fl_Widget* o, void* v );
  void pectoralisMinorMenuItem_CB_i();

  static void pectoralisMajorMenuItem_CB( Fl_Widget* o, void* v );
  void pectoralisMajorMenuItem_CB_i();

  static void subcutaneousFatMenuItem_CB( Fl_Widget* o, void* v );
  void subcutaneousFatMenuItem_CB_i();

  static void visceralFatMenuItem_CB( Fl_Widget* o, void* v );
  void visceralFatMenuItem_CB_i();

  static void obliqueFissureMenuItem_CB( Fl_Widget* o, void* v );
  void obliqueFissureMenuItem_CB_i();

  static void horizontalFissureMenuItem_CB( Fl_Widget* o, void* v );
  void horizontalFissureMenuItem_CB_i();

  static void anteriorScaleneMenuItem_CB( Fl_Widget* o, void* v );
  void anteriorScaleneMenuItem_CB_i();

  static void mildCentrilobularMenuItem_CB( Fl_Widget* o, void* v );
  void mildCentrilobularMenuItem_CB_i();

  static void moderateCentrilobularMenuItem_CB( Fl_Widget* o, void* v );
  void moderateCentrilobularMenuItem_CB_i();

  static void severeCentrilobularMenuItem_CB( Fl_Widget* o, void* v );
  void severeCentrilobularMenuItem_CB_i();

  static void panlobularMenuItem_CB( Fl_Widget* o, void* v );
  void panlobularMenuItem_CB_i();

  static void paraseptalMenuItem_CB( Fl_Widget* o, void* v );
  void paraseptalMenuItem_CB_i();
};
#endif
