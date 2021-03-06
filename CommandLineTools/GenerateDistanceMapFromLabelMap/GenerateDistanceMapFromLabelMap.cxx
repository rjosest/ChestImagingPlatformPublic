/** \file
 *  \ingroup commandLineTools 
 *  \details This program can be used to compute a distance map from an
 *  input label map (that adheres to the CIP label map conventions
 *  laid out in cipConventions.h). The user must specify which
 *  structure of interest the distance map should be computed with
 *  respect to by indicating the chest region and/or chest type. If
 *  the chest type is not specified, any voxel meeting the indicated
 *  chest type will be set to foreground and vice versa. The user also
 *  has the option of downsampling the label map prior to distance map
 *  computation, which should speed computation time. The resulting
 *  distance map will by upsampled by the same amount before writing.
 *  
 * USAGE:
 *
 * GenerateDistanceMapFromLabelMap.exe  [-t <unsigned int>] 
 *                                      [-r <unsigned int>] 
 *                                      [-s <double>] -d <string> 
 *                                      -l <string> [--] [--version]
 *                                      [-h]
 *
 * Where:
 *   -t <unsigned int>,  --type <unsigned int>
 *     Specify the chest type of the object the distance map is to be
 *     computed with respect to. UNDEFINEDTYPE by default
 *
 *   -r <unsigned int>,  --region <unsigned int>
 *     Specify the chest region of the object the distance map is to be
 *     computed with respect to. UNDEFINEDREGION by default
 *
 *   -s <double>,  --downsample <double>
 *     Downsample factor. The input label map will be downsampled by the
 *     specified amount before the distance map is computed. The resulting
 *     distance map will then be scaled up by the same amount before
 *     writing.
 *
 *   -d <string>,  --distanceMap <string>
 *     (required)  Output distance map file name
 *
 *   -l <string>,  --labelMap <string>
 *     (required)  Input label map file name
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
 *  $Date: 2012-09-06 15:50:17 -0400 (Thu, 06 Sep 2012) $
 *  $Revision: 247 $
 *  $Author: rjosest $
 *
 */

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <GenerateDistanceMapFromLabelMapCLP.h>
#include "cipConventions.h"
#include "itkResampleImageFilter.h"
#include "itkImage.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkIdentityTransform.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkCIPExtractChestLabelMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"

namespace{
    typedef itk::Image< short, 3 >                                                        DistanceMapType;
    typedef itk::ImageFileWriter< DistanceMapType >                                       WriterType;
    typedef itk::ImageFileWriter< cip::LabelMapType >                                     DEBWriterType;
    typedef itk::SignedMaurerDistanceMapImageFilter< cip::LabelMapType, DistanceMapType > SignedMaurerType;
    typedef itk::DanielssonDistanceMapImageFilter< cip::LabelMapType, DistanceMapType >   DanielssonType;
    typedef itk::NearestNeighborInterpolateImageFunction< cip::LabelMapType, double >     NearestNeighborInterpolatorType;
    typedef itk::LinearInterpolateImageFunction< DistanceMapType, double >                LinearInterpolatorType;
    typedef itk::ResampleImageFilter< cip::LabelMapType, cip::LabelMapType >              LabelMapResampleType;
    typedef itk::ResampleImageFilter< DistanceMapType, DistanceMapType >                  DistanceMapResampleType;
    typedef itk::ImageRegionIteratorWithIndex< cip::LabelMapType >                        LabelMapIteratorType;
    typedef itk::ImageRegionIteratorWithIndex< DistanceMapType >                          DistanceMapIteratorType;
    typedef itk::IdentityTransform< double, 3 >                                           IdentityType;
    typedef itk::CIPExtractChestLabelMapImageFilter                                       LabelMapExtractorType;

    cip::LabelMapType::Pointer ResampleImage( cip::LabelMapType::Pointer image, float downsampleFactor )
    {
        cip::LabelMapType::SizeType inputSize = image->GetBufferedRegion().GetSize();
    
        cip::LabelMapType::SpacingType inputSpacing = image->GetSpacing();
    
        cip::LabelMapType::SpacingType outputSpacing;
        outputSpacing[0] = inputSpacing[0]*downsampleFactor;
        outputSpacing[1] = inputSpacing[1]*downsampleFactor;
        outputSpacing[2] = inputSpacing[2]*downsampleFactor;
    
        cip::LabelMapType::SizeType outputSize;
        outputSize[0] = static_cast< unsigned int >( static_cast< double >( inputSize[0] )/downsampleFactor );
        outputSize[1] = static_cast< unsigned int >( static_cast< double >( inputSize[1] )/downsampleFactor );
        outputSize[2] = static_cast< unsigned int >( static_cast< double >( inputSize[2] )/downsampleFactor );
    
        NearestNeighborInterpolatorType::Pointer interpolator = NearestNeighborInterpolatorType::New();
        
        IdentityType::Pointer transform = IdentityType::New();
        transform->SetIdentity();
    
        LabelMapResampleType::Pointer resampler = LabelMapResampleType::New();
        resampler->SetTransform( transform );
        resampler->SetInterpolator( interpolator );
        resampler->SetInput( image );
        resampler->SetSize( outputSize );
        resampler->SetOutputSpacing( outputSpacing );
        resampler->SetOutputOrigin( image->GetOrigin() );
        try
        {
            resampler->Update();
        }
        catch ( itk::ExceptionObject &excp )
        {
            std::cerr << "Exception caught down sampling:";
            std::cerr << excp << std::endl;
        }
    
        return resampler->GetOutput();
    /*
     subsampledROIImage->SetRegions( resampler->GetOutput()->GetBufferedRegion().GetSize() );
     subsampledROIImage->Allocate();
     subsampledROIImage->FillBuffer( 0 );
     subsampledROIImage->SetSpacing( outputSpacing );
     subsampledROIImage->SetOrigin( image->GetOrigin() );
     
     LabelMapIteratorType rIt( resampler->GetOutput(), resampler->GetOutput()->GetBufferedRegion() );
     LabelMapIteratorType sIt( subsampledROIImage, subsampledROIImage->GetBufferedRegion() );
     
     rIt.GoToBegin();
     sIt.GoToBegin();
     while ( !sIt.IsAtEnd() )
     {
     sIt.Set( rIt.Get() );
     
     ++rIt;
     ++sIt;
     }
     */
    } 


    DistanceMapType::Pointer ResampleImage( DistanceMapType::Pointer image, float downsampleFactor )
    {
    DistanceMapType::SizeType inputSize = image->GetBufferedRegion().GetSize();
    
    DistanceMapType::SpacingType inputSpacing = image->GetSpacing();
    
    DistanceMapType::SpacingType outputSpacing;
    outputSpacing[0] = inputSpacing[0]*downsampleFactor;
    outputSpacing[1] = inputSpacing[1]*downsampleFactor;
    outputSpacing[2] = inputSpacing[2]*downsampleFactor;
    
    DistanceMapType::SizeType outputSize;
    outputSize[0] = static_cast< unsigned int >( static_cast< double >( inputSize[0] )/downsampleFactor );
    outputSize[1] = static_cast< unsigned int >( static_cast< double >( inputSize[1] )/downsampleFactor );
    outputSize[2] = static_cast< unsigned int >( static_cast< double >( inputSize[2] )/downsampleFactor );
    
    LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    
    IdentityType::Pointer transform = IdentityType::New();
    transform->SetIdentity();
    
    DistanceMapResampleType::Pointer resampler = DistanceMapResampleType::New();
    resampler->SetTransform( transform );
    resampler->SetInterpolator( interpolator );
    resampler->SetInput( image );
    resampler->SetSize( outputSize );
    resampler->SetOutputSpacing( outputSpacing );
    resampler->SetOutputOrigin( image->GetOrigin() );
    try
    {
        resampler->Update();
    }
    catch ( itk::ExceptionObject &excp )
    {
        std::cerr << "Exception caught down sampling:";
        std::cerr << excp << std::endl;
    }
    
    
    return resampler->GetOutput();
    
    } 

    #endif
};

int main( int argc, char *argv[] )
{
   PARSE_ARGS;

  //
  // Parse the input arguments
  //
     
    
    //parse args
    
    unsigned int cipRegion = (unsigned int) cipRegionTemp;
    unsigned int cipType = (unsigned int) cipTypeTemp;

  if ( cipType == cip::UNDEFINEDTYPE && cipRegion == cip::UNDEFINEDREGION )
    {
      std::cout << "Must specify a chest region or chest type" << std::endl;
      return cip::ARGUMENTPARSINGERROR;
    }

  // Instantiate ChestConventions for convenience
  cip::ChestConventions conventions;

  //
  // Read the label map image
  //
  std::cout << "Reading label map..." << std::endl;
  cip::LabelMapReaderType::Pointer reader = cip::LabelMapReaderType::New(); 
    reader->SetFileName( labelMapFileName );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught reading image:";
    std::cerr << excp << std::endl;

    return cip::LABELMAPREADFAILURE;
    }

  std::cout << "Isolationg region and type of interest..." << std::endl;
  LabelMapExtractorType::Pointer extractLabelMap = LabelMapExtractorType::New();
    extractLabelMap->SetInput(reader->GetOutput());
  if ( cipRegion != 0 )
    {
    extractLabelMap->SetChestRegion((unsigned char)cipRegion);
    }
  if ( cipType != 0 )
    {
    extractLabelMap->SetChestType((unsigned char)cipType);
    }
    extractLabelMap->Update();
  
  LabelMapIteratorType it( extractLabelMap->GetOutput(), extractLabelMap->GetOutput()->GetBufferedRegion() );

  it.GoToBegin();
  while ( !it.IsAtEnd() )
    {
    if ( it.Get() != 0 )
      {
	if ( it.Get() != 0 )
	  {
	    it.Set( 1 );
	  }
	else
	  {
	    it.Set( 0 );
	  }
      }
    
    ++it;
    }

  cip::LabelMapType::Pointer subSampledLabelMap;

  std::cout << "Downsampling label map..." << std::endl;
  subSampledLabelMap = ResampleImage( reader->GetOutput(), downsampleFactor );
  
  std::cout << "Generating distance map..." << std::endl;
  DanielssonType::Pointer distanceMap = DanielssonType::New();
    distanceMap->SetInput( subSampledLabelMap );
    distanceMap->InputIsBinaryOn();
    // distanceMap->SetSquaredDistance( 0 );
    // distanceMap->SetUseImageSpacing( 1 );
    // distanceMap->SetInsideIsPositive( interiorIsPositive );
  try
    {
    distanceMap->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught generating distance map:";
    std::cerr << excp << std::endl;

    return cip::GENERATEDISTANCEMAPFAILURE;
    }

  DistanceMapType::Pointer upSampledDistanceMap;

  std::cout << "Upsampling distance map..." << std::endl;
  upSampledDistanceMap = ResampleImage( distanceMap->GetOutput(), 1.0/downsampleFactor );

  std::cout << "Writing to file..." << std::endl;
  WriterType::Pointer writer = WriterType::New();
    writer->SetInput( upSampledDistanceMap );
    writer->SetFileName( distanceMapFileName );
    writer->UseCompressionOn();
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught writing image:";
    std::cerr << excp << std::endl;
    
    return cip::LABELMAPWRITEFAILURE;
    }

  std::cout << "DONE." << std::endl;

  return cip::EXITSUCCESS;
}

