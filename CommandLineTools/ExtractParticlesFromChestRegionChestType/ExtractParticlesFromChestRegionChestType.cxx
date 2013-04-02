/** \file
 *  \ingroup commandLineTools 
 *  \details This program allows you to extract particles from an
 *  input particles data set using either an input label map or the
 *  particles themselves. Many particles datasets contain a
 *  'ChestType' and 'ChestRegion' data field. The values of these
 *  fields can be used to isolate particles of
 *  interest. Alternatively, a label map can be specified, and only
 *  the particles falling in the specified region of interest will be
 *  written. Additionally, even if the input particles do not have the
 *  'ChestType' or 'ChestRegion' data arrays, the output particles
 *  data set will have these with region and type values specified at
 *  the command line.
 * 
 *  $Date: 2013-03-25 13:23:52 -0400 (Mon, 25 Mar 2013) $
 *  $Revision: 383 $
 *  $Author: jross $
 *
 *  USAGE: 
 *
 *  ExtractParticlesFromChestRegionChestType  [-t \<int\>] -r \<int\> 
 *                                            -o \<string\> -i \<string\>
 *                                            [-l \<string\>]  
 *                                            [--] [--version] [-h]
 *
 *  Where: 
 *
 *   -t \<int\>,  --type \<int\>
 *     Chest type for which to extract particles. If specifying a label map
 *     this flag shouldbe used to indicate the type of particles in the input
 *     file for output array labeling purposes (if no value is
 *     specifiedUNDEFINEDTYPE will be set as the particle ChestType field
 *     value
 *
 *   -r \<int\>,  --region \<int\>
 *     (required)  Chest region from which to extract particles
 *
 *   -o \<string\>,  --outParticles \<string\>
 *     (required)  Output particles file name
 *
 *   -i \<string\>,  --inParticles \<string\>
 *     (required)  Input particles file name
 *
 *   -l \<string\>,  --labelMap \<string\>
 *     Input label map file name. This is an optional input. If no label map
 *     is specified,the 'ChestRegion' and 'ChestType' arrays in the input
 *     will be used to extract theregion or type specified with the '-r' and
 *     '-t' flags, respectively
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


void GetOutputParticlesUsingLabelMap( std::string, unsigned char, unsigned char, vtkSmartPointer< vtkPolyData >, vtkSmartPointer< vtkPolyData > );
void GetOutputParticlesUsingChestRegionChestTypeArrays( std::vector< unsigned char >, std::vector< unsigned char >, 
							vtkSmartPointer< vtkPolyData >, vtkSmartPointer< vtkPolyData > );

int main( int argc, char *argv[] )
{
  //
  // Begin by defining the arguments to be passed
  //
  std::string labelMapFileName     = "NA";
  std::string inParticlesFileName  = "NA";
  std::string outParticlesFileName = "NA";
  std::vector< unsigned char > cipRegions;
  std::vector< unsigned char > cipTypes;

  //
  // Program and argument descriptions for user help
  //
  std::string programDescription = "This program allows you to extract particles from an input particles data set using either\
an input label map or the particles themselves. Many particles datasets contain a 'ChestType' and 'ChestRegion' data field.\
The values of these fields can be used to isolate particles of interest. Alternatively, a label map can be specified, and only\
the particles falling in the specified region of interest will be written. Additionally, even if the input particles do not have\
the ChestType' or 'ChestRegion' data arrays, the output particles data set will have these with region and type values specified\
at the command line.";

  std::string labelMapFileNameDescription     = "Input label map file name. This is an optional input. If no label map is specified,\
the 'ChestRegion' and 'ChestType' arrays in the input will be used to extract the\
region or type specified with the '-r' and '-t' flags, respectively";
  std::string inParticlesFileNameDescription  = "Input particles file name";
  std::string outParticlesFileNameDescription = "Output particles file name";
  std::string cipRegionsDescription            = "Chest regions from which to extract particles";
  std::string cipTypesDescription              = "Chest types for which to extract particles. If specifying a label map this flag should\
be used to indicate the type of particles in the input file for output array labeling purposes (if no value is specified\
UNDEFINEDTYPE will be set as the particle ChestType field value";

  //
  // Parse the input arguments
  //
  try
    {
    TCLAP::CmdLine cl( programDescription, ' ', "$Revision: 383 $" );

    TCLAP::ValueArg<std::string> labelMapFileNameArg( "l", "labelMap", labelMapFileNameDescription, false, labelMapFileName, "string", cl );
    TCLAP::ValueArg<std::string> inParticlesFileNameArg( "i", "inParticles", inParticlesFileNameDescription, true, inParticlesFileName, "string", cl );
    TCLAP::ValueArg<std::string> outParticlesFileNameArg( "o", "outParticles", outParticlesFileNameDescription, true, outParticlesFileName, "string", cl );
    TCLAP::MultiArg<int> cipRegionsArg( "r", "region", cipRegionsDescription, false, "int", cl );
    TCLAP::MultiArg<int> cipTypesArg( "t", "type", cipTypesDescription, false, "int", cl );

    cl.parse( argc, argv );

    labelMapFileName     = labelMapFileNameArg.getValue();
    inParticlesFileName  = inParticlesFileNameArg.getValue();
    outParticlesFileName = outParticlesFileNameArg.getValue();
    for ( unsigned int i=0; i<cipRegionsArg.getValue().size(); i++ )
      {
	cipRegions.push_back( (unsigned char)(cipRegionsArg.getValue()[i]) );
      }
    for ( unsigned int i=0; i<cipTypesArg.getValue().size(); i++ )
      {
	cipTypes.push_back( (unsigned char)(cipTypesArg.getValue()[i]) );
      }
    }
  catch ( TCLAP::ArgException excp )
    {
    std::cerr << "Error: " << excp.error() << " for argument " << excp.argId() << std::endl;
    return cip::ARGUMENTPARSINGERROR;
    }

  vtkSmartPointer< vtkPolyData > outParticles = vtkSmartPointer< vtkPolyData >::New();

  std::cout << "Reading polydata..." << std::endl;
  vtkSmartPointer< vtkPolyDataReader > particlesReader = vtkSmartPointer< vtkPolyDataReader >::New();
    particlesReader->SetFileName( inParticlesFileName.c_str() );
    particlesReader->Update();    

  if ( labelMapFileName.compare( "NA" ) != 0 )
    {
    // unsigned char tmpType;
    // if ( cipType == -1 )
    //   {
    //   tmpType = static_cast< unsigned char >( cip::UNDEFINEDTYPE );
    //   }
    // else 
    //   {
    //   tmpType = static_cast< unsigned char >( cipType );
    //   }
    // GetOutputParticlesUsingLabelMap( labelMapFileName, static_cast< unsigned char >( cipRegion ), tmpType, particlesReader->GetOutput(), outParticles );
    }
  else 
    {
    GetOutputParticlesUsingChestRegionChestTypeArrays( cipRegions, cipTypes, particlesReader->GetOutput(), outParticles );
    }

  std::cout << "Writing extracted particles..." << std::endl;
  vtkSmartPointer< vtkPolyDataWriter > particlesWriter = vtkSmartPointer< vtkPolyDataWriter >::New();
    particlesWriter->SetFileName( outParticlesFileName.c_str() );
    particlesWriter->SetInput( outParticles );
    particlesWriter->Write();

  std::cout << "DONE." << std::endl;

  return 0;
}


//
// Note that only the 'cipRegion' is used to isolate the
// particles. The 'cipType' is used here only to fill in information
// (in the 'ChestType' array) in the output.
//
// void GetOutputParticlesUsingLabelMap( std::string fileName, unsigned char cipRegion, unsigned char cipType, vtkSmartPointer< vtkPolyData > inParticles, 
//                                       vtkSmartPointer< vtkPolyData > outParticles )
// {
//   ChestConventions conventions;

//   std::cout << "Reading label map..." << std::endl;
//   ReaderType::Pointer labelMapReader = ReaderType::New();
//     labelMapReader->SetFileName( fileName );
//   try
//     {
//     labelMapReader->Update();
//     }
//   catch ( itk::ExceptionObject &excp )
//     {
//     std::cerr << "Exception caught reading label map:";
//     std::cerr << excp << std::endl;
//     }

//   unsigned int numberFieldDataArrays = inParticles->GetFieldData()->GetNumberOfArrays();
//   unsigned int numberParticles       = inParticles->GetNumberOfPoints();

//   vtkSmartPointer< vtkPoints > outputPoints  = vtkSmartPointer< vtkPoints >::New();
  
//   std::vector< vtkSmartPointer< vtkFloatArray > > arrayVec;

//   //
//   // The input particles may or may not have 'ChestType' and
//   // 'ChestRegion' data arrays. As we loop through the input, we will
//   // check their existence
//   //
//   bool foundChestRegionArray = false;
//   bool foundChestTypeArray   = false;

//   for ( unsigned int i=0; i<numberFieldDataArrays; i++ )
//     {
//     vtkSmartPointer< vtkFloatArray > array = vtkSmartPointer< vtkFloatArray >::New();
//       array->SetNumberOfComponents( inParticles->GetFieldData()->GetArray(i)->GetNumberOfComponents() );
//       array->SetName( inParticles->GetFieldData()->GetArray(i)->GetName() );

//     //
//     // The input particles may not have the 'ChestType' or
//     // 'ChestRegion' float arrays defined. If not, we want to define
//     // and add them
//     //
//     std::string name( inParticles->GetFieldData()->GetArray(i)->GetName() );
//     if ( name.compare( "ChestType" ) == 0 )
//       {
//       foundChestTypeArray = true;
//       }
//     if ( name.compare( "ChestRegion" ) == 0 )
//       {
//       foundChestRegionArray = true;
//       }

//     arrayVec.push_back( array );
//     }
//   if ( !foundChestRegionArray )
//     {
//     //
//     // The 'ChestRegion' data array was not found in the input, so add
//     // it here. Note that we have to increment the number of field
//     // data arrays for processing below
//     //
//     vtkSmartPointer< vtkFloatArray > array = vtkSmartPointer< vtkFloatArray >::New();
//       array->SetNumberOfComponents( 1 );
//       array->SetName( "ChestRegion" );

//     arrayVec.push_back( array );

//     numberFieldDataArrays++;
//     }
//   if ( !foundChestTypeArray )
//     {
//     //
//     // The 'ChestType' data array was not found in the input, so add
//     // it here. Note that we have to increment the number of field
//     // data arrays for processing below
//     //
//     vtkSmartPointer< vtkFloatArray > array = vtkSmartPointer< vtkFloatArray >::New();
//       array->SetNumberOfComponents( 1 );
//       array->SetName( "ChestType" );

//     arrayVec.push_back( array );

//     numberFieldDataArrays++;
//     }

//   unsigned int inc = 0;
//   ImageType::PointType point;
//   ImageType::IndexType index;

//   for ( unsigned int i=0; i<numberParticles; i++ )
//     {
//     point[0] = inParticles->GetPoint(i)[0];
//     point[1] = inParticles->GetPoint(i)[1];
//     point[2] = inParticles->GetPoint(i)[2];
    
//     labelMapReader->GetOutput()->TransformPhysicalPointToIndex( point, index );    

//     unsigned short labelValue   = labelMapReader->GetOutput()->GetPixel( index );
//     unsigned char  labelRegion  = conventions.GetChestRegionFromValue( labelValue );

//     //
//     // If the label map chest region is a subordinate of the requested
//     // chest region, then add this particle to the output
//     //
//     if ( conventions.CheckSubordinateSuperiorChestRegionRelationship( labelRegion, cipRegion ) )
//       {
//       outputPoints->InsertNextPoint( inParticles->GetPoint(i) );

//       for ( unsigned int j=0; j<numberFieldDataArrays; j++ )
//         {
//         //
//         // We have to perform this check because we don't know if the
//         // input has chest region and/or chest type arrays
//         //
//         std::string arrayName( arrayVec[j]->GetName() );

//         if ( arrayName.compare( "ChestRegion" ) == 0 )
//           {
//           float tmp = static_cast< float >( cipRegion );
          
//           arrayVec[j]->InsertTuple( inc, &tmp );
//           }
//         else if ( arrayName.compare( "ChestType" ) == 0 )
//           {
//           float tmp = static_cast< float >( cipType );

//           arrayVec[j]->InsertTuple( inc, &tmp );
//           }
//         else
//           {
//           arrayVec[j]->InsertTuple( inc, inParticles->GetFieldData()->GetArray(j)->GetTuple(i) );
//           }
//         }

//       inc++;
//       }
//     }
  
//   outParticles->SetPoints( outputPoints );
//   for ( unsigned int j=0; j<numberFieldDataArrays; j++ )
//     {
//     outParticles->GetFieldData()->AddArray( arrayVec[j] );
//     }
// }


void GetOutputParticlesUsingChestRegionChestTypeArrays( std::vector< unsigned char > cipRegions, std::vector< unsigned char > cipTypes, 
							vtkSmartPointer< vtkPolyData > inParticles, vtkSmartPointer< vtkPolyData > outParticles )
{
  unsigned int numberFieldDataArrays = inParticles->GetFieldData()->GetNumberOfArrays();
  unsigned int numberPointDataArrays = inParticles->GetPointData()->GetNumberOfArrays();
  unsigned int numberParticles       = inParticles->GetNumberOfPoints();

  vtkSmartPointer< vtkPoints > outputPoints  = vtkSmartPointer< vtkPoints >::New();
  
  std::vector< vtkSmartPointer< vtkFloatArray > > fieldArrayVec;
  std::vector< vtkSmartPointer< vtkFloatArray > > pointArrayVec;

  for ( unsigned int i=0; i<numberFieldDataArrays; i++ )
    {
    vtkSmartPointer< vtkFloatArray > array = vtkSmartPointer< vtkFloatArray >::New();
      array->SetNumberOfComponents( inParticles->GetFieldData()->GetArray(i)->GetNumberOfComponents() );
      array->SetName( inParticles->GetFieldData()->GetArray(i)->GetName() );

    fieldArrayVec.push_back( array );
    }
  for ( unsigned int i=0; i<numberPointDataArrays; i++ )
    {
    vtkSmartPointer< vtkFloatArray > array = vtkSmartPointer< vtkFloatArray >::New();
      array->SetNumberOfComponents( inParticles->GetPointData()->GetArray(i)->GetNumberOfComponents() );
      array->SetName( inParticles->GetPointData()->GetArray(i)->GetName() );

    pointArrayVec.push_back( array );
    }

  unsigned int inc = 0;
  for ( unsigned int i=0; i<numberParticles; i++ )
    {
      for ( unsigned int j=0; j<cipTypes.size(); j++ )
	{
	  if ( static_cast< int >( inParticles->GetFieldData()->GetArray( "ChestType" )->GetTuple(i)[0] ) == cipTypes[j] )
	    {
	      outputPoints->InsertNextPoint( inParticles->GetPoint(i) );

	      for ( unsigned int j=0; j<numberFieldDataArrays; j++ )
		{
		  fieldArrayVec[j]->InsertTuple( inc, inParticles->GetFieldData()->GetArray(j)->GetTuple(i) );
		}

	      inc++;
	    }
	}
    }

  outParticles->SetPoints( outputPoints );
  for ( unsigned int j=0; j<numberFieldDataArrays; j++ )
    {
    outParticles->GetFieldData()->AddArray( fieldArrayVec[j] );
    }
  for ( unsigned int j=0; j<numberPointDataArrays; j++ )
    {
    outParticles->GetPointData()->AddArray( pointArrayVec[j] );
    }
}


#endif