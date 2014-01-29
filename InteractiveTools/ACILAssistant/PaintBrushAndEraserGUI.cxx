// generated by Fast Light User Interface Designer (fluid) version 1.0107

#include "PaintBrushAndEraserGUI.h"
#include "cipConventions.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <iostream>

PaintBrushAndEraserGUI::PaintBrushAndEraserGUI() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = paintBrushAndEraserWindow = new Fl_Double_Window(380, 255, "Paint Brush & Eraser");
    w = o;
    o->box(FL_UP_BOX);
    o->user_data((void*)(this));
    { Fl_Group* o = new Fl_Group(13, 28, 357, 47, "Tool Type");
      o->box(FL_ENGRAVED_FRAME);
      o->align(FL_ALIGN_TOP_LEFT);
      { Fl_Round_Button* o = paintBrushRadio = new Fl_Round_Button(19, 28, 115, 42, "Paint Brush");
        o->type(102);
        o->down_box(FL_ROUND_DOWN_BOX);
        o->callback((Fl_Callback*)paintBrushRadio_CB, (void*)(this));
        paintBrushRadio->value(1);
      }
      { Fl_Round_Button* o = eraserRadio = new Fl_Round_Button(156, 33, 68, 32, "Eraser");
        o->type(102);
        o->down_box(FL_ROUND_DOWN_BOX);
        o->callback((Fl_Callback*)eraserRadio_CB, (void*)(this));
      }
      { Fl_Spinner* o = toolSizeSpinner = new Fl_Spinner(315, 38, 40, 25, "Tool Size");
        o->maximum( 99 );
      }
      o->end();
    }
    { Fl_Group* o = eraserGroup = new Fl_Group(13, 133, 357, 37);
      o->box(FL_ENGRAVED_BOX);
      { Fl_Button* o = eraseButton = new Fl_Button(246, 139, 110, 23, "Erase");
        o->callback((Fl_Callback*)eraseButton_CB, (void*)(this));
        eraseButton->deactivate();
      }
      { Fl_Check_Button* o = eraseSelectedCheck = new Fl_Check_Button(19, 142, 129, 15, "Erase Selected");
        o->down_box(FL_DOWN_BOX);
        eraseSelectedCheck->deactivate();
      }
      { Fl_Group* o = new Fl_Group(30, 150, 340, 20);
        o->end();
      }
      o->end();
    }
    { Fl_Group* o = toolThresholdsGroup = new Fl_Group(15, 87, 355, 38);
      o->box(FL_ENGRAVED_FRAME);
      o->labeltype(FL_ENGRAVED_LABEL);
      { Fl_Check_Button* o = tool3DCheck = new Fl_Check_Button(291, 96, 16, 20, "3D Tool");
        o->down_box(FL_DOWN_BOX);
      }
      { Fl_Value_Input* o = minValueInput = new Fl_Value_Input(91, 96, 46, 20, "Min Value");
        o->value(-1024);
      }
      { Fl_Value_Input* o = maxValueInput = new Fl_Value_Input(227, 96, 46, 20, "Max Value");
        o->value(1024);
      }
      o->end();
    }
    { Fl_Choice* o = chestRegionChoice = new Fl_Choice(103, 177, 269, 28, "Chest Region");
      o->down_box(FL_BORDER_BOX);
      //o->add("Undefined Region", 0, (Fl_Callback*)undefinedRegionMenuItem_CB, (void*)(this));
      o->add("Left Lung", 0, (Fl_Callback*)leftLungMenuItem_CB, (void*)(this));
      o->add("Right Lung", 0, (Fl_Callback*)rightLungMenuItem_CB, (void*)(this));
      o->add("Right Upper Lobe", 0, (Fl_Callback*)rightUpperLobeMenuItem_CB, (void*)(this));
      o->add("Right Middle Lobe", 0, (Fl_Callback*)rightMiddleLobeMenuItem_CB, (void*)(this));
      o->add("Right Lower Lobe", 0, (Fl_Callback*)rightLowerLobeMenuItem_CB, (void*)(this));
      o->add("Left Upper Lobe", 0, (Fl_Callback*)leftUpperLobeMenuItem_CB, (void*)(this));
      o->add("Left Lower Lobe", 0, (Fl_Callback*)leftLowerLobeMenuItem_CB, (void*)(this));
      o->add("Left", 0, (Fl_Callback*)leftMenuItem_CB, (void*)(this));
      o->add("Right", 0, (Fl_Callback*)rightMenuItem_CB, (void*)(this));
      o->add("Liver", 0, (Fl_Callback*)liverMenuItem_CB, (void*)(this));
      o->add("Spleen", 0, (Fl_Callback*)spleenMenuItem_CB, (void*)(this));
      o->add("Abdomen", 0, (Fl_Callback*)abdomen_CB, (void*)(this));
      o->add("Aorta", 0, (Fl_Callback*)aorta_CB, (void*)(this));
      o->add("Paravertebral", 0, (Fl_Callback*)paravertebral_CB, (void*)(this));
    }
    { Fl_Choice* o = chestTypeChoice = new Fl_Choice(103, 214, 269, 28, "Chest Type");
      o->down_box(FL_BORDER_BOX);
      o->add("Undefined Type", 0, (Fl_Callback*)undefinedTypeMenuItem_CB, (void*)(this));
      //o->add("Normal Parenchyma", 0, (Fl_Callback*)normalParenchymaMenuItem_CB, (void*)(this));
      o->add("Airway Generation 3", 0, (Fl_Callback*)airwayGeneration3MenuItem_CB, (void*)(this));
      o->add("Airway Generation 4", 0, (Fl_Callback*)airwayGeneration4MenuItem_CB, (void*)(this));
      o->add("Pec Minor", 0, (Fl_Callback*)pectoralisMinorMenuItem_CB, (void*)(this));
      o->add("Pec Major", 0, (Fl_Callback*)pectoralisMajorMenuItem_CB, (void*)(this));
      o->add("Subcutaneous Fat", 0, (Fl_Callback*)subcutaneousFatMenuItem_CB, (void*)(this));
      o->add("Visceral Fat", 0, (Fl_Callback*)visceralFatMenuItem_CB, (void*)(this));
      o->add("Oblique Fissure", 0, (Fl_Callback*)obliqueFissureMenuItem_CB, (void*)(this));
      o->add("Horizontal Fissure", 0, (Fl_Callback*)horizontalFissureMenuItem_CB, (void*)(this));
      o->add("Bronchiectatic Airway", 0, (Fl_Callback*)bronchiectaticAirwayMenuItem_CB, (void*)(this));
      o->add("Non-Bronchiectatic Airway", 0, (Fl_Callback*)nonBronchiectaticAirwayMenuItem_CB, (void*)(this));
      o->add("Ambiguous-Bronchiectatic Airway", 0, (Fl_Callback*)ambiguousBronchiectaticAirwayMenuItem_CB, (void*)(this));
      o->add("Muscle", 0, (Fl_Callback*)muscle_CB, (void*)(this));
      //o->add("Airway Generation 5", 0, (Fl_Callback*)airwayGeneration5MenuItem_CB, (void*)(this));
      //o->add("Mild Centrilobular Emphysema", 0, (Fl_Callback*)mildCentrilobularMenuItem_CB, (void*)(this));
      //o->add("Moderate Centrilobular Emphysema", 0, (Fl_Callback*)moderateCentrilobularMenuItem_CB, (void*)(this));
      //o->add("Severe Centrilobular Emphysema", 0, (Fl_Callback*)severeCentrilobularMenuItem_CB, (void*)(this));
      //o->add("Paraseptal Emphysema", 0, (Fl_Callback*)paraseptalMenuItem_CB, (void*)(this));
      //o->add("Panlobular Emphysema", 0, (Fl_Callback*)panlobularMenuItem_CB, (void*)(this));
    }
    o->end();
  }
  this->m_PaletteSelection =  static_cast< unsigned char >( cip::OBLIQUEFISSURE );
  this->m_LabelMapImage = LabelMapType::New();
  this->m_ChestRegion = static_cast< unsigned char >( cip::UNDEFINEDREGION );
  this->m_ChestType = static_cast< unsigned char >( cip::UNDEFINEDTYPE );
}

void PaintBrushAndEraserGUI::paintBrushRadio_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->paintBrushRadio_CB_i();
}

void PaintBrushAndEraserGUI::paintBrushRadio_CB_i() {
  this->eraseButton->deactivate();
this->eraseSelectedCheck->deactivate();
}

void PaintBrushAndEraserGUI::eraserRadio_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->eraserRadio_CB_i();
}

void PaintBrushAndEraserGUI::eraserRadio_CB_i() {
  this->eraseSelectedCheck->activate();
this->eraseButton->activate();
}

void PaintBrushAndEraserGUI::eraseButton_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->eraseButton_CB_i();
}

void PaintBrushAndEraserGUI::eraseButton_CB_i() 
{
  unsigned char cipRegionSelection = this->GetChestRegion();
  unsigned char cipTypeSelection   = this->GetChestType();

  cip::ChestConventions conventions;

  typedef itk::ImageRegionIteratorWithIndex< LabelMapType > LabelMapIteratorType;

  //
  // First delete painted indices as appropriate
  //
  bool stillDeleting = true;
  while ( stillDeleting )
    {
    stillDeleting = false;
    for ( unsigned int i=0; i < this->m_PaintedIndices->size(); i++ )
      {
      unsigned short label = this->m_LabelMapImage->GetPixel( (*this->m_PaintedIndices)[i] );
      
      unsigned char currentRegion = conventions.GetChestRegionFromValue( label );
      unsigned char currentType   = conventions.GetChestTypeFromValue( label );
      
      if ( this->eraseSelectedCheck->value() )
        {
        if ( (cipRegionSelection == static_cast< unsigned char >( cip::UNDEFINEDREGION ) && currentType == cipTypeSelection) ||
             (currentRegion == cipRegionSelection && cipTypeSelection == static_cast< unsigned char >( cip::UNDEFINEDTYPE )) )
          {
          this->m_PaintedIndices->erase( (*this->m_PaintedIndices).begin()+i );
          stillDeleting = true;
          break;
          }
        else if ( currentRegion == cipRegionSelection && currentType == cipTypeSelection )
          {
          this->m_PaintedIndices->erase( (*this->m_PaintedIndices).begin()+i );
          stillDeleting = true;
          break;
          }
        } 
      }  
    }
  if ( !this->eraseSelectedCheck->value() && this->m_PaintedIndices->size() > 0 )
    {
    this->m_PaintedIndices->clear();
    }

  //
  // Now update label map appropriately
  //
  LabelMapIteratorType lIt( this->m_LabelMapImage, this->m_LabelMapImage->GetBufferedRegion() );

  lIt.GoToBegin();
  while ( !lIt.IsAtEnd() )
    {
    if ( lIt.Get() != 0 )
      {
      unsigned short currentLabel = lIt.Get();
      unsigned char  currentRegion = 
        conventions.GetChestRegionFromValue( currentLabel );
      unsigned char  currentType = 
        conventions.GetChestTypeFromValue( currentLabel );
      
      if ( this->eraseSelectedCheck->value() )
        {
        if ( cipRegionSelection == static_cast< unsigned char >( cip::UNDEFINEDREGION ) && currentType == cipTypeSelection )
          {
          unsigned short newLabel = 
            conventions.GetValueFromChestRegionAndType( currentRegion, static_cast< unsigned char >( cip::UNDEFINEDTYPE ) );
          lIt.Set( newLabel );
          }
        else if ( currentRegion == cipRegionSelection && cipTypeSelection == static_cast< unsigned char >( cip::UNDEFINEDTYPE ) )
          {
          unsigned short newLabel = 
            conventions.GetValueFromChestRegionAndType( static_cast< unsigned char >( cip::UNDEFINEDREGION ), currentType );
          lIt.Set( newLabel );
          }
        else if ( currentRegion == cipRegionSelection && currentType == cipTypeSelection )
          {
          lIt.Set( 0 );
          }
        } 
      else
        {
        lIt.Set( 0 );
        }
      }
    
    ++lIt;	
    }  

//   if ( this->eraseSelectedCheck->value() )
//     {
//     (*this->m_PaintedIndices)[paletteSelection].clear();
//     }
//   else
//     {
//     for ( int i=0; i<256; i++ )
//       { 
//       (*this->m_PaintedIndices)[i].clear();
//       }
//     }

  this->m_UpdateViewerFunction();
}

unsigned char PaintBrushAndEraserGUI::GetPaletteSelection() {
  return this->m_PaletteSelection;
}

unsigned char PaintBrushAndEraserGUI::GetChestRegion() {
  return this->m_ChestRegion;
}

unsigned char PaintBrushAndEraserGUI::GetChestType() {
  return this->m_ChestType;
}

unsigned int PaintBrushAndEraserGUI::GetToolRadius() {
  unsigned int radius = this->toolSizeSpinner->value() - 1;

return radius;
}

bool PaintBrushAndEraserGUI::GetPaintBrushSelected() {
  return static_cast< bool >( this->paintBrushRadio->value() );
}

bool PaintBrushAndEraserGUI::GetEraserSelected() {
  return static_cast< bool >( this->eraserRadio->value() );
}

bool PaintBrushAndEraserGUI::GetEraseSelectedSelected() {
  return static_cast< bool >( this->eraseSelectedCheck->value() );
}

void PaintBrushAndEraserGUI::SetLabelMapImage( LabelMapType::Pointer labelMapImage ) {
  this->m_LabelMapImage = labelMapImage;
}

void PaintBrushAndEraserGUI::SetUpdateViewerFunction( void (*f)() ) {
  this->m_UpdateViewerFunction = f;
}

void PaintBrushAndEraserGUI::SetPaintedIndices( std::vector< LabelMapType::IndexType >* indices ) {
  this->m_PaintedIndices = indices;
}

short PaintBrushAndEraserGUI::GetToolLowerThreshold() {
  return static_cast< short >( this->minValueInput->value() );
}

short PaintBrushAndEraserGUI::GetToolUpperThreshold() {
  return static_cast< short >( this->maxValueInput->value() );
}

void PaintBrushAndEraserGUI::rightLungMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->rightLungMenuItem_CB_i();
}

void PaintBrushAndEraserGUI::rightLungMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::RIGHTLUNG );
}


void PaintBrushAndEraserGUI::undefinedRegionMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->undefinedRegionMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::undefinedRegionMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::UNDEFINEDREGION );
}


void PaintBrushAndEraserGUI::leftLungMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->leftLungMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::leftLungMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::LEFTLUNG );
}


void PaintBrushAndEraserGUI::rightUpperLobeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->rightUpperLobeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::rightUpperLobeMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::RIGHTSUPERIORLOBE );
}


void PaintBrushAndEraserGUI::rightMiddleLobeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->rightMiddleLobeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::rightMiddleLobeMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::RIGHTMIDDLELOBE );
}


void PaintBrushAndEraserGUI::rightLowerLobeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->rightLowerLobeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::rightLowerLobeMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::RIGHTINFERIORLOBE );
}


void PaintBrushAndEraserGUI::leftLowerLobeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->leftLowerLobeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::leftLowerLobeMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::LEFTINFERIORLOBE );
}


void PaintBrushAndEraserGUI::leftUpperLobeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->leftUpperLobeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::leftUpperLobeMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::LEFTSUPERIORLOBE );
}


void PaintBrushAndEraserGUI::leftMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->leftMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::leftMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::LEFT );
}

void PaintBrushAndEraserGUI::rightMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->rightMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::rightMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::RIGHT );
}

void PaintBrushAndEraserGUI::abdomen_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->abdomen_CB_i();
}
void PaintBrushAndEraserGUI::abdomen_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::ABDOMEN );
}

void PaintBrushAndEraserGUI::aorta_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->aorta_CB_i();
}
void PaintBrushAndEraserGUI::aorta_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::AORTA );
}

void PaintBrushAndEraserGUI::paravertebral_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->paravertebral_CB_i();
}
void PaintBrushAndEraserGUI::paravertebral_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::PARAVERTEBRAL );
}

void PaintBrushAndEraserGUI::muscle_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->muscle_CB_i();
}
void PaintBrushAndEraserGUI::muscle_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::MUSCLE );
  this->minValueInput->value(-50);
  this->maxValueInput->value(90);
}

void PaintBrushAndEraserGUI::liverMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->liverMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::liverMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::LIVER );
}

void PaintBrushAndEraserGUI::spleenMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->spleenMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::spleenMenuItem_CB_i() {
  this->m_ChestRegion = static_cast< unsigned char >( cip::SPLEEN );
}

void PaintBrushAndEraserGUI::undefinedTypeMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->undefinedTypeMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::undefinedTypeMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::UNDEFINEDTYPE );
}


void PaintBrushAndEraserGUI::normalParenchymaMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->normalParenchymaMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::normalParenchymaMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::NORMALPARENCHYMA );
}


void PaintBrushAndEraserGUI::airwayGeneration3MenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->airwayGeneration3MenuItem_CB_i();
}
void PaintBrushAndEraserGUI::airwayGeneration3MenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::AIRWAYGENERATION3 );
}


void PaintBrushAndEraserGUI::airwayGeneration4MenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->airwayGeneration4MenuItem_CB_i();
}
void PaintBrushAndEraserGUI::airwayGeneration4MenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::AIRWAYGENERATION4 );
}


void PaintBrushAndEraserGUI::airwayGeneration5MenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->airwayGeneration5MenuItem_CB_i();
}
void PaintBrushAndEraserGUI::airwayGeneration5MenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::AIRWAYGENERATION5 );
}


void PaintBrushAndEraserGUI::pectoralisMinorMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->pectoralisMinorMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::pectoralisMinorMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::PECTORALISMINOR );
  this->minValueInput->value(-50);
  this->maxValueInput->value(90);
}

void PaintBrushAndEraserGUI::pectoralisMajorMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->pectoralisMajorMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::pectoralisMajorMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::PECTORALISMAJOR );
  this->minValueInput->value(-50);
  this->maxValueInput->value(90);
}

void PaintBrushAndEraserGUI::subcutaneousFatMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->subcutaneousFatMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::subcutaneousFatMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::SUBCUTANEOUSFAT );
  this->minValueInput->value(-200);
  this->maxValueInput->value(0); 
}

void PaintBrushAndEraserGUI::visceralFatMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->visceralFatMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::visceralFatMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::VISCERALFAT );
  this->minValueInput->value(-250);
  this->maxValueInput->value(-50);
}

void PaintBrushAndEraserGUI::obliqueFissureMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->obliqueFissureMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::obliqueFissureMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::OBLIQUEFISSURE );
}

void PaintBrushAndEraserGUI::ambiguousBronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->ambiguousBronchiectaticAirwayMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::ambiguousBronchiectaticAirwayMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::AMBIGUOUSBRONCHIECTATICAIRWAY );
}

void PaintBrushAndEraserGUI::nonBronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->nonBronchiectaticAirwayMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::nonBronchiectaticAirwayMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::NONBRONCHIECTATICAIRWAY );
}

void PaintBrushAndEraserGUI::bronchiectaticAirwayMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->bronchiectaticAirwayMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::bronchiectaticAirwayMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::BRONCHIECTATICAIRWAY );
}

void PaintBrushAndEraserGUI::horizontalFissureMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->horizontalFissureMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::horizontalFissureMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::HORIZONTALFISSURE );
}

void PaintBrushAndEraserGUI::anteriorScaleneMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->anteriorScaleneMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::anteriorScaleneMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::ANTERIORSCALENE );
}


void PaintBrushAndEraserGUI::paraseptalMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->paraseptalMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::paraseptalMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::PARASEPTALEMPHYSEMA );
}


void PaintBrushAndEraserGUI::panlobularMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->panlobularMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::panlobularMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::PANLOBULAREMPHYSEMA );
}


void PaintBrushAndEraserGUI::mildCentrilobularMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->mildCentrilobularMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::mildCentrilobularMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::MILDCENTRILOBULAREMPHYSEMA );
}


void PaintBrushAndEraserGUI::moderateCentrilobularMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->moderateCentrilobularMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::moderateCentrilobularMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::MODERATECENTRILOBULAREMPHYSEMA );
}


void PaintBrushAndEraserGUI::severeCentrilobularMenuItem_CB( Fl_Widget* o, void* v ) {
  ((PaintBrushAndEraserGUI*)v)->severeCentrilobularMenuItem_CB_i();
}
void PaintBrushAndEraserGUI::severeCentrilobularMenuItem_CB_i() {
  this->m_ChestType = static_cast< unsigned char >( cip::SEVERECENTRILOBULAREMPHYSEMA );
}

