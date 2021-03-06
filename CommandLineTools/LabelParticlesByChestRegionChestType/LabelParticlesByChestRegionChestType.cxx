/** \file
 *  \ingroup commandLineTools 
 *  \details This program is used to label particles datasets by chest
 *  region and chest type. The user must specify the type of the input
 *  particles, but the chest region can either be determined by an
 *  input label map or be specified at the command line. 
 * 
 *  $Date: 2012-07-18 13:23:15 -0700 (Wed, 18 Jul 2012) $
 *  $Revision: 196 $
 *  $Author: rjosest $
 *
 *  USAGE: 
 *
 *   LabelParticlesByChestRegionChestType  [-t \<unsigned char\>] [-r
 *                                         \<unsigned char\>] -o \<string\> -i
 *                                         \<string\> [-l \<string\>] [--]
 *                                         [--version] [-h]
 *
 *  Where: 
 *   -t \<unsigned char\>,  --type \<unsigned char\>
 *     Chest type for particles labeling. UNDEFINEDTYPE by default
 *
 *   -r \<unsigned char\>,  --region \<unsigned char\>
 *     Chest region for particles labeling. UNDEFINEDREGION by default
 *
 *   -o \<string\>,  --outParticles \<string\>
 *     (required)  Output particles file name
 *
 *   -i \<string\>,  --inParticles \<string\>
 *     (required)  Input particles file name
 *
 *   -l \<string\>,  --labelMap \<string\>
 *     Input label map file name. If specified the 'ChestRegion'value will be
 *     determined from the label map
 *
 *   --,  --ignore_rest
 *     Ignores the rest of the labeled arguments following this flag.
 *
 *   --version
 *     Displays version information and exits.
 *
 *   -h,  --help
 *     Displays usage information and exits.
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <tclap/CmdLine.h>
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "cipConventions.h"


typedef itk::Image< unsigned short, 3 >      ImageType;
typedef itk::ImageFileReader< ImageType >    ReaderType;


void InitializeParticleChestRegionChestTypeArrays( vtkSmartPointer< vtkPolyData > );


int main( int argc, char *argv[] )
{
  //
  // Begin by defining the arguments to be passed
  //
  std::string labelMapFileName     = "NA";
  std::string inParticlesFileName  = "NA";
  std::string outParticlesFileName = "NA";
  float       cipRegion            = static_cast< float >( cip::UNDEFINEDREGION );
  float       cipType              = static_cast< float >( cip::UNDEFINEDTYPE );

  //
  // Program and argument descriptions for user help
  //
  std::string programDescription = "This program is used to label particles datasets by chest region and \
chest type. The user must specify the type of the input particles, but the chest region \
can either be determined by an input label map or be specified at the command line. ";

  std::string labelMapFileNameDescription     = "Input label map file name. If specified the 'ChestRegion'\
value will be determined from the label map";
  std::string inParticlesFileNameDescription  = "Input particles file name";
  std::string outParticlesFileNameDescription = "Output particles file name";
  std::string cipRegionDescription            = "Chest region for particles labeling. UNDEFINEDREGION by default";
  std::string cipTypeDescription              = "Chest type for particles labeling. UNDEFINEDTYPE by default";

  //
  // Parse the input arguments
  //
  try
    {
    TCLAP::CmdLine cl( programDescription, ' ', "$Revision: 196 $" );

    TCLAP::ValueArg<std::string> labelMapFileNameArg( "l", "labelMap", labelMapFileNameDescription, false, labelMapFileName, "string", cl );
    TCLAP::ValueArg<std::string> inParticlesFileNameArg( "i", "inParticles", inParticlesFileNameDescription, true, inParticlesFileName, "string", cl );
    TCLAP::ValueArg<std::string> outParticlesFileNameArg( "o", "outParticles", outParticlesFileNameDescription, true, outParticlesFileName, "string", cl );
    TCLAP::ValueArg<float> cipRegionArg( "r", "region", cipRegionDescription, false, cipRegion, "unsigned char", cl );
    TCLAP::ValueArg<float> cipTypeArg( "t", "type", cipTypeDescription, false, cipType, "unsigned char", cl );

    cl.parse( argc, argv );

    labelMapFileName     = labelMapFileNameArg.getValue();
    inParticlesFileName  = inParticlesFileNameArg.getValue();
    outParticlesFileName = outParticlesFileNameArg.getValue();
    cipRegion            = cipRegionArg.getValue();
    cipType              = cipTypeArg.getValue();

    }
  catch ( TCLAP::ArgException excp )
    {
    std::cerr << "Error: " << excp.error() << " for argument " << excp.argId() << std::endl;
    return cip::ARGUMENTPARSINGERROR;
    }

  //
  // Instantiate ChestConventions for later use
  //
  cip::ChestConventions conventions;

  //
  // Read the particles
  //
  std::cout << "Reading polydata..." << std::endl;
  vtkSmartPointer< vtkPolyDataReader > particlesReader = vtkSmartPointer< vtkPolyDataReader >::New();
    particlesReader->SetFileName( inParticlesFileName.c_str() );
    particlesReader->Update();    

  //
  // Initialize chest region and chest type field data arrays. We
  // don't assume that the incoming particles have these arrays. If
  // they don't we add them. If they do, we initialize them with
  // 'UNDEFINEDREGION' and 'UNDEFINEDTYPE'.
  //
  std::cout << "Initializing chest-region chest-type arrays..." << std::endl;
  InitializeParticleChestRegionChestTypeArrays( particlesReader->GetOutput() );

  //
  // If specified, read the input label map
  //
  if ( labelMapFileName.compare( "NA" ) != 0 )
    {
    std::cout << "Reading label map..." << std::endl;
    ReaderType::Pointer labelMapReader = ReaderType::New();
      labelMapReader->SetFileName( labelMapFileName );
    try
      {
      labelMapReader->Update();
      }
    catch ( itk::ExceptionObject &excp )
      {
      std::cerr << "Exception caught reading label map:";
      std::cerr << excp << std::endl;     

      return cip::LABELMAPREADFAILURE;
      }

    ImageType::PointType point;
    ImageType::IndexType index;

    //
    // Loop through the particles to label them
    //
    for ( unsigned int i=0; i<particlesReader->GetOutput()->GetNumberOfPoints(); i++ )
      {
      point[0] = particlesReader->GetOutput()->GetPoint(i)[0];
      point[1] = particlesReader->GetOutput()->GetPoint(i)[1];
      point[2] = particlesReader->GetOutput()->GetPoint(i)[2];
      
      labelMapReader->GetOutput()->TransformPhysicalPointToIndex( point, index );    

      unsigned short labelValue = labelMapReader->GetOutput()->GetPixel( index );

      cipRegion  = static_cast< float >( conventions.GetChestRegionFromValue( labelValue ) );

      particlesReader->GetOutput()->GetPointData()->GetArray( "ChestRegion" )->SetTuple( i, &cipRegion );
      particlesReader->GetOutput()->GetPointData()->GetArray( "ChestType" )->SetTuple( i, &cipType );
      }
    }
  else
    {
    //
    // If here, no label map was specified, and we must assign region
    // and types based on user specification. Loop through the
    // particles to label them 
    //
    for ( unsigned int i=0; i<particlesReader->GetOutput()->GetNumberOfPoints(); i++ )
      {
      particlesReader->GetOutput()->GetPointData()->GetArray( "ChestRegion" )->SetTuple( i, &cipRegion );
      particlesReader->GetOutput()->GetPointData()->GetArray( "ChestType" )->SetTuple( i, &cipType );
      }
    }

  //
  // Write the labeled particles
  //
  std::cout << "Writing labeled particles..." << std::endl;
  vtkSmartPointer< vtkPolyDataWriter > particlesWriter = vtkSmartPointer< vtkPolyDataWriter >::New();
    particlesWriter->SetFileName( outParticlesFileName.c_str() );
    particlesWriter->SetInput( particlesReader->GetOutput() );
    particlesWriter->SetFileTypeToBinary();
    particlesWriter->Update();    

  std::cout << "DONE." << std::endl;

  return 0;
}


void InitializeParticleChestRegionChestTypeArrays( vtkSmartPointer< vtkPolyData > particles )
{
  unsigned int numberPointDataArrays = particles->GetPointData()->GetNumberOfArrays();
  unsigned int numberParticles       = particles->GetNumberOfPoints();

  //
  // The input particles may or may not have 'ChestType' and
  // 'ChestRegion' data arrays. As we loop through the input, we will
  // check their existence
  //
  bool foundChestRegionArray = false;
  bool foundChestTypeArray   = false;

  for ( unsigned int i=0; i<numberPointDataArrays; i++ )
    {
    std::string name( particles->GetPointData()->GetArray(i)->GetName() );
    if ( name.compare( "ChestType" ) == 0 )
      {
      foundChestTypeArray = true;
      }
    if ( name.compare( "ChestRegion" ) == 0 )
      {
      foundChestRegionArray = true;
      }
    }

  //
  // The 'chestRegionArray' and 'chestTypeArray' defined here will
  // only be used if they have not already been found in 'particles' 
  //
  vtkSmartPointer< vtkFloatArray > chestRegionArray = vtkSmartPointer< vtkFloatArray >::New();
    chestRegionArray->SetNumberOfComponents( 1 );
    chestRegionArray->SetName( "ChestRegion" );

  vtkSmartPointer< vtkFloatArray > chestTypeArray = vtkSmartPointer< vtkFloatArray >::New();
    chestTypeArray->SetNumberOfComponents( 1 );
    chestTypeArray->SetName( "ChestType" );
  
  //
  // Now loop through the particles to initialize the arrays
  //
  float cipRegion = static_cast< float >( cip::UNDEFINEDREGION );
  float cipType   = static_cast< float >( cip::UNDEFINEDTYPE );

  for ( unsigned int i=0; i<numberParticles; i++ )
    {
    if ( foundChestRegionArray )
      {
      particles->GetPointData()->GetArray( "ChestRegion" )->SetTuple( i, &cipRegion );
      }
    else
      {
      chestRegionArray->InsertTuple( i, &cipRegion );
      }

    if ( foundChestTypeArray )
      {
      particles->GetPointData()->GetArray( "ChestType" )->SetTuple( i, &cipType );
      }
    else
      {
      chestTypeArray->InsertTuple( i, &cipType );
      }
    }

  if ( !foundChestRegionArray )
    {
    particles->GetPointData()->AddArray( chestRegionArray );
    }
  if ( !foundChestTypeArray )
    {
    particles->GetPointData()->AddArray( chestTypeArray );
    }
}


#endif
