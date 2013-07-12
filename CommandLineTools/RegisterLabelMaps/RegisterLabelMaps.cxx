#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkCenteredTransformInitializer.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"
#include "itkLungConventions.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkExtractLungLabelMapImageFilter.h"


typedef itk::Image< unsigned short, 3 >                                             ImageType;
typedef itk::ResampleImageFilter< ImageType, ImageType >                            ResampleFilterType;
typedef itk::ImageFileReader< ImageType >                                           ImageReaderType;
typedef itk::ImageFileWriter< ImageType >                                           ImageWriterType;
typedef itk::RegularStepGradientDescentOptimizer                                    OptimizerType;
typedef itk::ImageRegistrationMethod< ImageType, ImageType >                        RegistrationType;
typedef itk::KappaStatisticImageToImageMetric< ImageType, ImageType >               MetricType;
typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >           InterpolatorType;
typedef itk::AffineTransform< double, 3 >                                           TransformType;
typedef itk::CenteredTransformInitializer< TransformType, ImageType, ImageType >    InitializerType;
typedef OptimizerType::ScalesType                                                   OptimizerScalesType;
typedef itk::ImageRegionIteratorWithIndex< ImageType >                              IteratorType;
typedef itk::RegionOfInterestImageFilter< ImageType, ImageType >                    RegionOfInterestType;
typedef itk::ResampleImageFilter< ImageType, ImageType >                            ResampleType;
typedef itk::IdentityTransform< double, 3 >                                         IdentityType;
typedef itk::ExtractLungLabelMapImageFilter                                         LabelMapExtractorType;


void ResampleImage( ImageType::Pointer, ImageType::Pointer, float );
void WriteTransformFile( TransformType::Pointer, char* );


void usage()
{
  std::cerr << "\n";
  std::cerr << "Usage: Registers <options> where <options> is one or more " << std::endl;
  std::cerr << "of the following:\n\n";
  std::cerr << "   <-h>     Display (this) usage information\n";
  std::cerr << "   <-if>    Fixed image file name\n";
  std::cerr << "   <-im>    Moving image file name\n";
  std::cerr << "   <-d>     Down sample factor. The fixed and moving images will be down sampled\n";
  std::cerr << "            by this amount before registration. (Default is 1.0, i.e. no down sampling).\n";
  std::cerr << "   <-max>   Max step length. Defaults is 1.0\n";
  std::cerr << "   <-min>   Min step length. Default is 0.001\n";
  std::cerr << "   <-n>     Number of iterations. Default is 200\n";
  std::cerr << "   <-t>     Translation scale. Default is 0.001\n";
  std::cerr << "   <-oi>    Output image file name\n";
  std::cerr << "   <-ot>    Output transform file name\n";

  exit(1);
}


int main( int argc, char *argv[] )
{
  bool ok;

  char* fixedImageFileName               = new char[512];  strcpy( fixedImageFileName, "q" );
  char* movingImageFileName              = new char[512];  strcpy( movingImageFileName, "q" );
  char* outputImageFileName              = new char[512];  strcpy( outputImageFileName, "q" ); 
  char* outputTransformFileName          = new char[512];  strcpy( outputTransformFileName, "q" );
  float maxStepLength                    = 1.0;
  float minStepLength                    = 0.001;
  int   numberOfIterations               = 200;
  float translationScale                 = 0.001;
  float downsampleFactor                 = 1.0;

  while ( argc > 1 )
    {
    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-h") == 0))
      {
      argc--; argv++;
      ok = true;
      usage();      
      }

    if ((ok == false) && (strcmp(argv[1], "-ot") == 0))
      {
      argc--; argv++;
      ok = true;

      outputTransformFileName = argv[1];

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-oi") == 0))
      {
      argc--; argv++;
      ok = true;

      outputImageFileName = argv[1];

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-t") == 0))
      {
      argc--; argv++;
      ok = true;

      translationScale = atof( argv[1] );

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-d") == 0))
      {
      argc--; argv++;
      ok = true;

      downsampleFactor = atof( argv[1] );

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-n") == 0))
      {
      argc--; argv++;
      ok = true;

      numberOfIterations = atoi( argv[1] );

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-min") == 0))
      {
      argc--; argv++;
      ok = true;

      minStepLength = atof( argv[1] );

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-max") == 0))
      {
      argc--; argv++;
      ok = true;

      maxStepLength = atof( argv[1] );

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-im") == 0))
      {
      argc--; argv++;
      ok = true;

      movingImageFileName = argv[1];

      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-if") == 0))
      {
      argc--; argv++;
      ok = true;

      fixedImageFileName = argv[1];

      argc--; argv++;
      }

    }

  ImageType::Pointer subSampledFixedImage = ImageType::New();
  {
  std::cout << "Reading fixed image..." << std::endl;
  ImageReaderType::Pointer fixedReader = ImageReaderType::New();
    fixedReader->SetFileName( fixedImageFileName );
  try
    {
    fixedReader->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while updating fixed reader:";
    std::cerr << excp << std::endl;
    }

  std::cout << "Subsampling fixed image..." << std::endl;
  ResampleImage( fixedReader->GetOutput(), subSampledFixedImage, downsampleFactor );
  }

  ImageType::Pointer subSampledMovingImage = ImageType::New();
  {
  std::cout << "Reading moving image..." << std::endl;
  ImageReaderType::Pointer movingReader = ImageReaderType::New();
    movingReader->SetFileName( movingImageFileName );
  try
    {
    movingReader->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while updating moving reader:";
    std::cerr << excp << std::endl;
    }

  std::cout << "Subsampling moving image..." << std::endl;
  ResampleImage( movingReader->GetOutput(), subSampledMovingImage, downsampleFactor );
  }

  //-------
  // Now extract the whole lung region from the resized fixed and
  // moving images
  //
  std::cout << "Extracting whole lung region from down sampled fixed image..." << std::endl;
  LabelMapExtractorType::Pointer fixedExtractor = LabelMapExtractorType::New();
    fixedExtractor->SetInput( subSampledFixedImage );
    fixedExtractor->SetLungRegion( static_cast< unsigned char >( WHOLELUNG ) );
    fixedExtractor->Update();

  std::cout << "Extracting whole lung region from down sampled moving image..." << std::endl;
  LabelMapExtractorType::Pointer movingExtractor = LabelMapExtractorType::New();
    movingExtractor->SetInput( subSampledMovingImage );
    movingExtractor->SetLungRegion( static_cast< unsigned char >( WHOLELUNG ) );
    movingExtractor->Update();

  //-------
  // Set up the registration framework
  //
  LungConventions conventions;
  
  unsigned short foregroundValue = conventions.GetValueFromLungRegionAndType( static_cast< unsigned char >( WHOLELUNG ), 
                                                                              static_cast< unsigned char >( UNDEFINEDTYPE) );

  MetricType::Pointer metric = MetricType::New();
    metric->ComplementOn();
    metric->SetForegroundValue( foregroundValue );

  TransformType::Pointer transform = TransformType::New();

  InitializerType::Pointer initializer = InitializerType::New();
    initializer->SetTransform( transform );
    initializer->SetFixedImage(  fixedExtractor->GetOutput() );
    initializer->SetMovingImage( movingExtractor->GetOutput() );
    initializer->MomentsOn();
    initializer->InitializeTransform();

  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
    optimizerScales[0] =  1.0;   optimizerScales[1] =  1.0;   optimizerScales[2] =  1.0;
    optimizerScales[3] =  1.0;   optimizerScales[4] =  1.0;   optimizerScales[5] =  1.0;
    optimizerScales[6] =  1.0;   optimizerScales[7] =  1.0;   optimizerScales[8] =  1.0;
    optimizerScales[9]  =  translationScale;
    optimizerScales[10] =  translationScale;
    optimizerScales[11] =  translationScale;

  OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetScales( optimizerScales );
    optimizer->SetMaximumStepLength( maxStepLength );
    optimizer->SetMinimumStepLength( minStepLength ); 
    optimizer->SetNumberOfIterations( numberOfIterations );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  std::cout << "Starting registration..." << std::endl;
  RegistrationType::Pointer registration = RegistrationType::New();  
    registration->SetMetric( metric );
    registration->SetOptimizer( optimizer );
    registration->SetInterpolator( interpolator );
    registration->SetTransform( transform );
    registration->SetFixedImage( fixedExtractor->GetOutput() );
    registration->SetMovingImage( movingExtractor->GetOutput() );
    registration->SetFixedImageRegion( fixedExtractor->GetOutput()->GetBufferedRegion() );
    registration->SetInitialTransformParameters( transform->GetParameters() );
  try 
    { 
    registration->StartRegistration(); 
    } 
  catch( itk::ExceptionObject &excp ) 
    { 
    std::cerr << "ExceptionObject caught while executing registration" << std::endl; 
    std::cerr << excp << std::endl; 
    } 

  OptimizerType::ParametersType finalParams = registration->GetLastTransformParameters();

  TransformType::Pointer finalTransform = TransformType::New();
    finalTransform->SetParameters( finalParams );
    finalTransform->SetCenter( transform->GetCenter() );

  if ( strcmp(outputTransformFileName, "q") != 0 )
    {
    std::cout << "Writing transform..." << std::endl;
    WriteTransformFile( finalTransform, outputTransformFileName );
    }

  if ( strcmp(outputImageFileName, "q") != 0 )
    {
    std::cout << "Resampling moving image..." << std::endl;
    ResampleFilterType::Pointer resample = ResampleFilterType::New();
      resample->SetTransform( finalTransform );
      resample->SetInput( movingExtractor->GetOutput() );
      resample->SetSize( fixedExtractor->GetOutput()->GetBufferedRegion().GetSize() );
      resample->SetOutputOrigin(  fixedExtractor->GetOutput()->GetOrigin() );
      resample->SetOutputSpacing( fixedExtractor->GetOutput()->GetSpacing() );
      resample->SetInterpolator( interpolator );
      resample->SetDefaultPixelValue( 0 );
      resample->Update();

    ImageType::Pointer upsampledImage = ImageType::New();

    std::cout << "Upsampling to original size..." << std::endl;
    ResampleImage( resample->GetOutput(), upsampledImage, 1.0/downsampleFactor );

    ImageWriterType::Pointer writer = ImageWriterType::New();
      writer->SetInput( upsampledImage );
      writer->SetFileName( outputImageFileName );
      writer->UseCompressionOn();
    try
      {
      writer->Update();
      }
    catch ( itk::ExceptionObject &excp )
      {
      std::cerr << "Exception caught writing output image:";
      std::cerr << excp << std::endl;
      }
    }    

  std::cout << "DONE." << std::endl;

  return 0;
}


void WriteTransformFile( TransformType::Pointer transform, char* fileName )
{
  itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
    transformWriter->SetInput( transform );
    transformWriter->SetFileName( fileName );
  try
    {
    transformWriter->Update();
    }
  catch ( itk::ExceptionObject &excp )
    {
    std::cerr << "Exception caught while updating transform writer:";
    std::cerr << excp << std::endl;
    }
}


void ResampleImage( ImageType::Pointer image, ImageType::Pointer subsampledROIImage, float downsampleFactor )
{
  ImageType::SizeType inputSize = image->GetBufferedRegion().GetSize();

  ImageType::SpacingType inputSpacing = image->GetSpacing();

  ImageType::SpacingType outputSpacing;
    outputSpacing[0] = inputSpacing[0]*downsampleFactor;
    outputSpacing[1] = inputSpacing[1]*downsampleFactor;
    outputSpacing[2] = inputSpacing[2]*downsampleFactor;

  ImageType::SizeType outputSize;
    outputSize[0] = static_cast< unsigned int >( static_cast< double >( inputSize[0] )/downsampleFactor );
    outputSize[1] = static_cast< unsigned int >( static_cast< double >( inputSize[1] )/downsampleFactor );
    outputSize[2] = static_cast< unsigned int >( static_cast< double >( inputSize[2] )/downsampleFactor );

  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  IdentityType::Pointer transform = IdentityType::New();
    transform->SetIdentity();

  ResampleType::Pointer resampler = ResampleType::New();
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

  subsampledROIImage->SetRegions( resampler->GetOutput()->GetBufferedRegion().GetSize() );
  subsampledROIImage->Allocate();
  subsampledROIImage->FillBuffer( 0 );
  subsampledROIImage->SetSpacing( outputSpacing );
  subsampledROIImage->SetOrigin( image->GetOrigin() );

  IteratorType rIt( resampler->GetOutput(), resampler->GetOutput()->GetBufferedRegion() );
  IteratorType sIt( subsampledROIImage, subsampledROIImage->GetBufferedRegion() );

  rIt.GoToBegin();
  sIt.GoToBegin();
  while ( !sIt.IsAtEnd() )
    {
    sIt.Set( rIt.Get() );
    
    ++rIt;
    ++sIt;
    }
} 
