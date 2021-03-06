/** \file
 *  \ingroup commandLineTools 
 *  \details This program can be used to transfer the contents of a VTK 
 *  polydata's field data to point data and vice-versa. Generally, field 
 *  data applies to a dataset as a whole and need not have a one-to-one 
 *  correspondence with the points. However, this may be the case in some 
 *  instances (esp. with the particles datasets). In those cases it may be 
 *  helpful to have the data contained in field data arrays also stored in 
 *  point data arrays (e.g. for rendering purposes). Field data will only 
 *  be transferred provided that the number of tuples in the field data 
 *  array is the same as the number of points.
 * 
 *  USAGE: 
 *
 *  TransferFieldDataToFromPointData  [--mp <bool>] [--mf <bool>] [--pf
 *                                     <bool>] [--fp <bool>] -o <string> -i
 *                                     <string> [--] [--version] [-h]
 *
 *  Where: 
 *
 *  --mp <bool>
 *   Setting this to true will maintain the field data. Setting it to false
 *   will eliminate the field data from the output. Only relevant if
 *   requesting to transfer field data to point data
 *
 *  --mf <bool>
 *   Setting this to true will maintain the field data. Setting it to false
 *   will eliminate the field data from the output. Only relevant if
 *   requesting to transfer field data to point data
 *
 *  --pf <bool>
 *   Set to true to transfer point data to field data
 *
 *  --fp <bool>
 *   Set to true to transfer field data to point data
 *
 *  -o <string>,  --output <string>
 *   (required)  Output VTK polydata file name
 *
 *  -i <string>,  --input <string>
 *   (required)  Input VTK polydata file name
 *
 *  --,  --ignore_rest
 *   Ignores the rest of the labeled arguments following this flag.
 *
 *  --version
 *   Displays version information and exits.
 *
 *  -h,  --help
 *   Displays usage information and exits.
 *  
 *  $Date: 2013-03-25 13:23:52 -0400 (Mon, 25 Mar 2013) $
 *  $Revision: 383 $
 *  $Author: jross $
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <tclap/CmdLine.h>
#include "cipConventions.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkPolyData.h"
#include "vtkIndent.h"

int main( int argc, char *argv[] )
{
  //
  // Define arguments
  //
  std::string inFileName    = "NA";
  std::string outFileName   = "NA";
  bool        maintainField = true;
  bool        maintainPoint = true;
  bool        fieldToPoint  = true;
  bool        pointToField  = false;

  //
  // Program and argument descriptions for user help
  //
  std::string programDesc = "This program can be used to transfer the contents \
of a VTK polydata's field data to point data and vice-versa. Generally, field \
data applies to a dataset as a whole and need not have a one-to-one correspondence \
with the points. However, this may be the case in some instances (esp. with the \
particles datasets). In those cases it may be helpful to have the data contained in \
field data arrays also stored in point data arrays (e.g. for rendering purposes). \
Field data will only be transferred provided that the number of tuples in the \
field data array is the same as the number of points.";

  std::string inFileNameDesc  = "Input VTK polydata file name";
  std::string outFileNameDesc = "Output VTK polydata file name";
  std::string maintainFieldDesc  = "Setting this to true will maintain the field data. \
Setting it to false will eliminate the field data from the output. Only relevant if \
requesting to transfer field data to point data";
  std::string maintainPointDesc  = "Setting this to true will maintain the field data. \
Setting it to false will eliminate the field data from the output. Only relevant if \
requesting to transfer field data to point data";
  std::string fieldToPointDesc = "Set to true to transfer field data to point data";
  std::string pointToFieldDesc = "Set to true to transfer point data to field data";

  //
  // Parse the input arguments
  //
  try
    {
    TCLAP::CmdLine cl( programDesc, ' ', "$Revision: 383 $" );

    TCLAP::ValueArg<std::string> inFileNameArg( "i", "input", inFileNameDesc, true, inFileName, "string", cl );
    TCLAP::ValueArg<std::string> outFileNameArg( "o", "output", outFileNameDesc, true, outFileName, "string", cl );
    TCLAP::ValueArg<bool> fieldToPointArg( "", "fp", fieldToPointDesc, false, fieldToPoint, "bool", cl );
    TCLAP::ValueArg<bool> pointToFieldArg( "", "pf", pointToFieldDesc, false, pointToField, "bool", cl );
    TCLAP::ValueArg<bool> maintainFieldArg( "", "mf", maintainFieldDesc, false, maintainField, "bool", cl );
    TCLAP::ValueArg<bool> maintainPointArg( "", "mp", maintainPointDesc, false, maintainPoint, "bool", cl );

    cl.parse( argc, argv );

    inFileName    = inFileNameArg.getValue();
    outFileName   = outFileNameArg.getValue();
    fieldToPoint  = fieldToPointArg.getValue();
    pointToField  = pointToFieldArg.getValue();
    maintainField = maintainFieldArg.getValue();
    maintainPoint = maintainPointArg.getValue();
    }
  catch ( TCLAP::ArgException excp )
    {
    std::cerr << "Error: " << excp.error() << " for argument " << excp.argId() << std::endl;
    return cip::ARGUMENTPARSINGERROR;
    }

  //
  // Read the poly data
  //
  std::cout << "Reading VTK polydata..." << std::endl;
  vtkSmartPointer< vtkPolyDataReader > reader = vtkSmartPointer< vtkPolyDataReader >::New();
    reader->SetFileName( inFileName.c_str() );
    reader->Update();
  
  unsigned int numberOfPoints = reader->GetOutput()->GetNumberOfPoints();
  unsigned int numberOfFieldDataArrays = reader->GetOutput()->GetFieldData()->GetNumberOfArrays();
  unsigned int numberOfPointDataArrays = reader->GetOutput()->GetPointData()->GetNumberOfArrays();

  vtkSmartPointer< vtkPoints > outputPoints = vtkSmartPointer< vtkPoints >::New();

  //
  // Create a new polydata to contain the output
  //
  vtkSmartPointer< vtkPolyData > outPolyData = vtkSmartPointer< vtkPolyData >::New();
    outPolyData->SetPoints( reader->GetOutput()->GetPoints() );

  //
  // First create a list of the point data already stored in the input polydata.
  // We will only transfer field data provided there is not already corresponding
  // point data. Do the same for the field data array names.
  //
  std::vector< std::string > pointDataArrayNames;
  for ( unsigned int i=0; i<numberOfPointDataArrays; i++ )
    {
      std::string name = reader->GetOutput()->GetPointData()->GetArray(i)->GetName();
      pointDataArrayNames.push_back(name);
    }
  std::vector< std::string > fieldDataArrayNames;
  for ( unsigned int i=0; i<numberOfFieldDataArrays; i++ )
    {
      std::string name = reader->GetOutput()->GetFieldData()->GetArray(i)->GetName();
      fieldDataArrayNames.push_back(name);
    }

  //
  // Transfer the field data to point data if requested
  //
  if ( fieldToPoint )
    {
      bool alreadyPresent;
      for ( unsigned int i=0; i<numberOfFieldDataArrays; i++ )
	{
	  alreadyPresent = false;
	  std::string fieldDataArrayName = reader->GetOutput()->GetFieldData()->GetArray(i)->GetName();
	  for ( unsigned int j=0; j<numberOfPointDataArrays; j++ )
	    {
	      if ( fieldDataArrayName.compare(pointDataArrayNames[j]) == 0 )
		{
		  alreadyPresent = true;
		  break;
		}
	    }

	  if ( !alreadyPresent )
	    {
	      //
	      // The number of array tuples must be the same as the number of points
	      //
	      if ( reader->GetOutput()->GetFieldData()->GetArray(i)->GetNumberOfTuples() == numberOfPoints )
		{
		  outPolyData->GetPointData()->AddArray( reader->GetOutput()->GetFieldData()->GetArray(i) );
		}
	    }
	}
    }

  //
  // Transfer the point data to field data if requested. Note that this is not the
  // "preferred" direction to go in, but we support it so that code in the cip
  // that depends on field data storage can be accomodated. Eventually, all point-
  // specific data should be recorded as point data.
  //
  if ( pointToField )
    {
      bool alreadyPresent;
      for ( unsigned int i=0; i<numberOfPointDataArrays; i++ )
	{
	  alreadyPresent = false;
	  std::string pointDataArrayName = reader->GetOutput()->GetPointData()->GetArray(i)->GetName();
	  for ( unsigned int j=0; j<numberOfFieldDataArrays; j++ )
	    {
	      if ( pointDataArrayName.compare(fieldDataArrayNames[j]) == 0 )
		{
		  alreadyPresent = true;
		  break;
		}
	    }

	  if ( !alreadyPresent )
	    {
	      outPolyData->GetFieldData()->AddArray( reader->GetOutput()->GetPointData()->GetArray(i) );
	    }
	}
    }

  //
  // Add the field data to the output if requested
  //
  if ( maintainField )
    {
      for ( unsigned int i=0; i<numberOfFieldDataArrays; i++ )
	{
	  outPolyData->GetFieldData()->AddArray( reader->GetOutput()->GetFieldData()->GetArray(i) );
	}
    }

  //
  // Add the point data to the output if requested
  //
  if ( maintainPoint )
    {
      for ( unsigned int i=0; i<numberOfPointDataArrays; i++ )
	{
	  outPolyData->GetPointData()->AddArray( reader->GetOutput()->GetPointData()->GetArray(i) );
	}
    }

  //
  // Write the poly data
  //
  std::cout << "Writing VTK polydata..." << std::endl;
  vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
    writer->SetFileName( outFileName.c_str() );
    writer->SetInput( outPolyData );
    writer->Update();

  std::cout << "DONE." << std::endl;

  return 0;
}

#endif
