
/** \file
 *  \ingroup commandLineTools
 *  \details This program registers 2 label maps, source and target, and
 * a transformation file as well as the transformed image
 *
 *  USAGE: ./RegisterLabelMaps --regionVec 1 -m

 /net/th914_nas.bwh.harvard.edu/mnt/array1/share/Processed/COPDGene/11622T/11622T_INSP_STD_HAR_COPD/11622T_INSP_STD_HAR_COPD_leftLungRightLung.nhdr
 -f
 /net/th914_nas.bwh.harvard.edu/mnt/array1/share/Processed/COPDGene/10393Z/10393Z_INSP_STD_HAR_COPD/10393Z_INSP_STD_HAR_COPD_leftLungRightLung.nhdr
 --outputImage /projects/lmi/people/rharmo/projects/dataInfo/testoutput.nrrd
 --outputTransform output_transform_file -d 12

 *
 *  $Date: $
 *  $Revision: $
 *  $Author:  $
 *
 */

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegistrationMethod.h"
#include "itkCenteredTransformInitializer.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkKappaStatisticImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkAffineTransform.h"
#include "itkTransformFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkCIPExtractChestLabelMapImageFilter.h"
#include "RegisterLabelMaps2DCLP.h"
#include "cipConventions.h"
#include "cipHelper.h"
#include <sstream>
#include "itkQuaternionRigidTransform.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include <fstream>
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkCIPExtractChestLabelMapImageFilter.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xinclude.h>
#include <libxml/xmlIO.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include "itkNormalizeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageMaskSpatialObject.h"
//#include <uuid/uuid.h>

namespace
{
#define MY_ENCODING "ISO-8859-1"
  typedef itk::Image< unsigned short, 2 >                                                           LabelMapType2D; 
  typedef itk::RegularStepGradientDescentOptimizer                                                  OptimizerType;
    
    typedef itk::GradientDescentOptimizer GradOptimizerType;
    
  typedef OptimizerType::ScalesType                                                                 OptimizerScalesType;
 
  typedef itk::IdentityTransform< double, 2 >                                                       IdentityType;
  typedef itk::CIPExtractChestLabelMapImageFilter                                                   LabelMapExtractorType;
 
  typedef itk::GDCMImageIO                                                                          ImageIOType;
  typedef itk::GDCMSeriesFileNames                                                                  NamesGeneratorType;
 
  typedef itk::Image< short, 2 >                                                                    ShortImageType;
  typedef itk::ResampleImageFilter< LabelMapType2D, LabelMapType2D >                                ResampleFilterType;
  
    
  typedef itk::AffineTransform<double, 2 >                                                          TransformType2D;
  typedef itk::QuaternionRigidTransform<double>                                                     RigidTransformType;
  typedef itk::CenteredTransformInitializer< TransformType2D, LabelMapType2D, LabelMapType2D >  InitializerType2D;
  typedef itk::CenteredTransformInitializer< TransformType2D, ShortImageType, ShortImageType >  InitializerType2DIntensity;
    
  typedef itk::ImageFileReader< LabelMapType2D >                                               LabelMap2DReaderType;
  typedef itk::ImageFileReader< ShortImageType >                                               ShortReaderType;
    
  typedef itk::ImageRegistrationMethod<LabelMapType2D,LabelMapType2D >                         RegistrationType;
  typedef itk::ImageRegistrationMethod<ShortImageType,ShortImageType >                         CTRegistrationType;
    
  typedef itk::KappaStatisticImageToImageMetric<LabelMapType2D, LabelMapType2D >               MetricType;
  typedef itk::NormalizedCorrelationImageToImageMetric<ShortImageType, ShortImageType  >       ncMetricType;
    
  typedef itk::NearestNeighborInterpolateImageFunction< LabelMapType2D, double >               InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction< ShortImageType, double >               CTInterpolatorType;
    
    
  typedef itk::ImageRegionIteratorWithIndex< LabelMapType2D >                                  IteratorType;
  typedef itk::ImageFileWriter< LabelMapType2D >                                               ImageWriterType;
  //typedef itk::RegionOfInterestImageFilter< LabelMapType2D, LabelMapType2D >                   RegionOfInterestType;
  typedef itk::ResampleImageFilter< LabelMapType2D, LabelMapType2D >                           ResampleType;
  typedef itk::ImageRegionIteratorWithIndex< LabelMapType2D >                                  LabelMapIteratorType;
    
    typedef itk::DiscreteGaussianImageFilter< ShortImageType, ShortImageType > GaussianFilterType;
    
    typedef itk::NormalizeImageFilter<ShortImageType,ShortImageType> FixedNormalizeFilterType;
    typedef itk::NormalizeImageFilter<ShortImageType,ShortImageType> MovingNormalizeFilterType;
    typedef itk::ImageMaskSpatialObject< 2 >   MaskType;
   


  struct REGIONTYPEPAIR
  {
    unsigned char region;
    unsigned char type;
  };

  struct REGISTRATION_XML_DATA
  {
    std::string registrationID;
    float similarityValue;
    std::string transformationLink;
    std::string sourceID;
    std::string destID;
    std::string similarityMeasure;
    std::string image_type;
    int transformationIndex;
  };

  void WriteTransformFile( TransformType2D::Pointer transform, char* fileName )
  {
    itk::TransformFileWriter::Pointer transformWriter =
      itk::TransformFileWriter::New();
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
 

  LabelMapType2D::Pointer ReadLabelMap2DFromFile( std::string
						   labelMapFileName )
  {
    std::cout << "Reading label map..." << std::endl;
    LabelMap2DReaderType::Pointer reader = LabelMap2DReaderType::New();
    reader->SetFileName( labelMapFileName );
    try
      {
      reader->Update();
      }
    catch ( itk::ExceptionObject &excp )
      {
          std::cerr << "Exception caught reading label map:";
          std::cerr << excp << std::endl;
      }

    return reader->GetOutput();
  }

    
  ShortImageType::Pointer ReadCTFromFile( std::string fileName )
    {
        ShortReaderType::Pointer reader = ShortReaderType::New();
        reader->SetFileName( fileName );
        try
        {
            reader->Update();
        }
        catch ( itk::ExceptionObject &excp )
        {
            std::cerr << "Exception caught reading CT image:";
            std::cerr << excp << std::endl;
            return NULL;
        }
        
        return reader->GetOutput();
    }

  void WriteRegistrationXML(const char *file, REGISTRATION_XML_DATA
			    &theXMLData)
  {
    std::cout<<"Writing registration XML file"<<std::endl;
    xmlDocPtr doc = NULL;       /* document pointer */
    xmlNodePtr root_node = NULL; /* Node pointers */
    xmlDtdPtr dtd = NULL;       /* DTD pointer */

    doc = xmlNewDoc(BAD_CAST "1.0");
    root_node = xmlNewNode(NULL, BAD_CAST "Registration");
    xmlDocSetRootElement(doc, root_node);

    dtd = xmlCreateIntSubset(doc, BAD_CAST "root", NULL, BAD_CAST
			     "RegistrationOutput_v2.dtd");

    time_t timer;
    time(&timer);
    std::stringstream tempStream;
    tempStream<<timer;
    theXMLData.registrationID.assign("Registration");
    theXMLData.registrationID.append(tempStream.str());
    theXMLData.registrationID.append("_");
    theXMLData.registrationID.append(theXMLData.sourceID.c_str());
    theXMLData.registrationID.append("_to_");
    theXMLData.registrationID.append(theXMLData.destID.c_str());


    xmlNewProp(root_node, BAD_CAST "Registration_ID", BAD_CAST
	       (theXMLData.registrationID.c_str()));

    // xmlNewChild() creates a new node, which is "attached"
    // as child node of root_node node.
    std::ostringstream similaritString;
    //std::string tempsource;
    similaritString <<theXMLData.similarityValue;
    std::ostringstream transformationIndexString;    
    transformationIndexString <<theXMLData.transformationIndex;
    
    xmlNewChild(root_node, NULL, BAD_CAST "image_type", BAD_CAST
                (theXMLData.image_type.c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "transformation", BAD_CAST
		(theXMLData.transformationLink.c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "transformation_index", BAD_CAST
                (transformationIndexString.str().c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "movingID", BAD_CAST
		(theXMLData.sourceID.c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "fixedID", BAD_CAST
		(theXMLData.destID.c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "SimilarityMeasure", BAD_CAST
		(theXMLData.similarityMeasure.c_str()));
    xmlNewChild(root_node, NULL, BAD_CAST "SimilarityValue", BAD_CAST
		(similaritString.str().c_str()));
    xmlSaveFormatFileEnc(file, doc, "UTF-8", 1);
    xmlFreeDoc(doc);

  }


} //end namespace

int main( int argc, char *argv[] )
{

  std::vector< unsigned char >  regionVec;
  std::vector< unsigned char >  typeVec;
  std::vector< unsigned char >  regionPairVec;
  std::vector< unsigned char >  typePairVec;
    
  std::cout << "about to parse args" << std::endl;  

  PARSE_ARGS;

  std::cout << "agrs parsed, new" << std::endl;

    
  //Read in region and type pair
  std::vector< REGIONTYPEPAIR > regionTypePairVec;

  for ( unsigned int i=0; i<regionVecArg.size(); i++ )
    {
      regionVec.push_back(regionVecArg[i]);
    }
  for ( unsigned int i=0; i<typeVecArg.size(); i++ )
    {
      typeVec.push_back( typeVecArg[i] );
    }
  if (regionPairVec.size() == typePairVecArg.size())
    {
      for ( unsigned int i=0; i<regionPairVecArg.size(); i++ )
      {
          REGIONTYPEPAIR regionTypePairTemp;

          regionTypePairTemp.region = regionPairVecArg[i];
          argc--; argv++;
          regionTypePairTemp.type   = typePairVecArg[i];

          regionTypePairVec.push_back( regionTypePairTemp );
      }
    }

  //Read in fixed image label map from file and subsample
  LabelMapType2D::Pointer fixedLabelMap2D = LabelMapType2D::New();
  LabelMapType2D::Pointer movingLabelMap2D =LabelMapType2D::New();
  ShortImageType::Pointer fixedCT2D = ShortImageType::New();
  ShortImageType::Pointer movingCT2D = ShortImageType::New();
    
  if(isIntensity != true)
  {
      if ( strcmp( fixedImageFileName.c_str(), "q") != 0 )
      {
          std::cout << "Reading label map from file..." << std::endl;

          fixedLabelMap2D = ReadLabelMap2DFromFile( fixedImageFileName );
          if (fixedLabelMap2D.GetPointer() == NULL)
          {
              return cip::LABELMAPREADFAILURE;
          }
      }
      else
      {
          std::cerr <<"Error: No lung label map specified"<< std::endl;
          return cip::EXITFAILURE;
      }

 
      //Read in moving image label map from file and subsample
  
      if ( strcmp( movingImageFileName.c_str(), "q") != 0 )
      {
          std::cout << "Reading label map from file..." << std::endl;
          movingLabelMap2D = ReadLabelMap2DFromFile( movingImageFileName );

          if (movingLabelMap2D.GetPointer() == NULL)
          {
              return cip::LABELMAPREADFAILURE;
          }
      }
      else
      {
          std::cerr <<"Error: No lung label map specified"<< std::endl;
          return cip::EXITFAILURE;
      }
  }
    
  else //intensity based registration, load CT images
  {
      if ( strcmp( fixedImageFileName.c_str(), "q") != 0 )
        {
            std::cout << "Reading CT from file..." << std::endl;
            
            fixedCT2D = ReadCTFromFile( fixedImageFileName );
            if (fixedCT2D.GetPointer() == NULL)
            {
                return cip::LABELMAPREADFAILURE;
            }
            
        }
        else
        {
            std::cerr <<"Error: No CT specified"<< std::endl;
            return cip::EXITFAILURE;
        }
        
        //Read in moving image label map from file and subsample
        
        if ( strcmp( movingImageFileName.c_str(), "q") != 0 )
        {
            std::cout << "Reading CT from file..." << std::endl;
            
            movingCT2D = ReadCTFromFile( movingImageFileName );
            
            if (movingCT2D.GetPointer() == NULL)
            {
                return cip::LABELMAPREADFAILURE;
            }
            
        }
        else
        {
            std::cerr <<"Error: No CT specified"<< std::endl;
            return cip::EXITFAILURE;
        }
  }
    

    MetricType::Pointer metric = MetricType::New();
    metric->SetForegroundValue( 1);    //because we are minimizing as opposed to maximizing
    ncMetricType::Pointer nc_metric = ncMetricType::New();

    typedef itk::Image< unsigned char, 2 >   ImageMaskType;
     MaskType::Pointer  spatialObjectMask = MaskType::New();
    typedef itk::ImageFileReader< ImageMaskType >    MaskReaderType;
    MaskReaderType::Pointer  maskReader = MaskReaderType::New();
    
    if ( strcmp( fixedLabelmapFileName.c_str(), "q") != 0 )
    {
        std::cout<<"reading fixed label map "<<fixedLabelmapFileName.c_str() <<std::endl;
        maskReader->SetFileName(fixedLabelmapFileName.c_str() );
        
        try
        {
            maskReader->Update();
        }
        catch( itk::ExceptionObject & err )
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }
        spatialObjectMask->SetImage(maskReader->GetOutput());
        nc_metric->SetFixedImageMask( spatialObjectMask );
        
    }

        
    
  TransformType2D::Pointer transform2D = TransformType2D::New();
  RigidTransformType::Pointer rigidTransform = RigidTransformType::New();
  std::cout<<"initializing transform"<<std::endl;
    

    
  InitializerType2D::Pointer initializer2D = InitializerType2D::New();
  InitializerType2DIntensity::Pointer initializer2DIntensity = InitializerType2DIntensity::New();
   
    
  OptimizerType::Pointer optimizer = OptimizerType::New();
  GradOptimizerType::Pointer grad_optimizer = GradOptimizerType::New();
   
  FixedNormalizeFilterType::Pointer fixedNormalizer =FixedNormalizeFilterType::New();
  MovingNormalizeFilterType::Pointer movingNormalizer =MovingNormalizeFilterType::New();
    
  GaussianFilterType::Pointer fixedSmoother  = GaussianFilterType::New();
  GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();
  fixedSmoother->SetVariance( 1.5 );
  movingSmoother->SetVariance( 1.5 );
    
  if(isIntensity != true)
  {
      initializer2D->SetTransform( transform2D );
      initializer2D->SetFixedImage(fixedLabelMap2D );
      initializer2D->SetMovingImage( movingLabelMap2D);
      initializer2D->MomentsOn();
      initializer2D->InitializeTransform(); //this makes it work
  }
  else
  {

      fixedSmoother->SetInput( fixedCT2D );
      fixedSmoother->Update();
      movingSmoother->SetInput( movingCT2D);
      movingSmoother->Update();
      
      fixedNormalizer->SetInput(  fixedSmoother->GetOutput() );
      movingNormalizer->SetInput( movingSmoother->GetOutput() );
      fixedNormalizer->Update();
      movingNormalizer->Update();
      initializer2DIntensity->SetTransform( transform2D );
      initializer2DIntensity->SetFixedImage(fixedNormalizer->GetOutput());
      initializer2DIntensity->SetMovingImage(movingNormalizer->GetOutput()  );
      initializer2DIntensity->MomentsOn();
      initializer2DIntensity->InitializeTransform();
      //optimizer->MaximizeOn();

  }

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  CTInterpolatorType::Pointer CTinterpolator = CTInterpolatorType::New();      

  std::cout << "Starting registration..." << std::endl;

  RegistrationType::Pointer registration = RegistrationType::New();
  CTRegistrationType::Pointer CTregistration = CTRegistrationType::New();
  double bestValue;
  TransformType2D::Pointer finalTransform2D = TransformType2D::New();
    
  if (isIntensity == true)
  {
      //grad_optimizer->SetLearningRate( 0.2 );
      CTregistration->SetMetric( nc_metric );
      CTregistration->SetFixedImage( fixedSmoother->GetOutput() );
      CTregistration->SetMovingImage(movingSmoother->GetOutput());
      
      CTregistration->SetOptimizer( optimizer );
      CTregistration->SetInterpolator( CTinterpolator );
      CTregistration->SetTransform( transform2D );
      
      CTregistration->SetInitialTransformParameters( transform2D->GetParameters());
      
      try
      {
          CTregistration->Initialize();
          CTregistration->Update();
      }
      catch( itk::ExceptionObject &excp )
      {
          std::cerr << "ExceptionObject caught while executing registration" <<
          std::endl;
          
          std::cerr << excp << std::endl;
      }
      
      //get all params to output to file
      numberOfIterations = optimizer->GetCurrentIteration();
      bestValue = optimizer->GetValue();
      
      std::cout <<" best similarity value = " <<bestValue<<std::endl;
      
      GradOptimizerType::ParametersType finalParams;
      finalParams =CTregistration->GetLastTransformParameters();
      finalTransform2D->SetParameters( finalParams );
  }
      
  else
  {
      registration->SetMetric( metric );
      registration->SetFixedImage( fixedLabelMap2D );
      registration->SetMovingImage(movingLabelMap2D);
      
      registration->SetOptimizer( optimizer );
      registration->SetInterpolator( interpolator );
      registration->SetTransform( transform2D );
      
      registration->SetInitialTransformParameters( transform2D->GetParameters());
      
      try
      {
          registration->Initialize();
          registration->Update();
      }
      catch( itk::ExceptionObject &excp )
      {
          std::cerr << "ExceptionObject caught while executing registration" <<
          std::endl;
          
          std::cerr << excp << std::endl;
      }
      
      //get all params to output to file
      numberOfIterations = optimizer->GetCurrentIteration();
      bestValue = optimizer->GetValue();
      
      std::cout <<" best similarity value = " <<bestValue<<std::endl;
      
      OptimizerType::ParametersType finalParams;
      finalParams =registration->GetLastTransformParameters();
      finalTransform2D->SetParameters( finalParams );
  }



    finalTransform2D->SetCenter( transform2D->GetCenter() );
    
    std::cout<<"writing final transform"<<std::endl;

  if ( strcmp(outputTransformFileName.c_str(), "q") != 0 )
    {
      std::string infoFilename = outputTransformFileName;
      int result = infoFilename.find_last_of('.');
      if (std::string::npos != result)
	  infoFilename.erase(result);
      // append extension:
      infoFilename.append(".xml");

      REGISTRATION_XML_DATA labelMapRegistrationXMLData;
      labelMapRegistrationXMLData.similarityValue = (float)(bestValue);
      const char *similarity_type = metric->GetNameOfClass();
      labelMapRegistrationXMLData.similarityMeasure.assign(similarity_type);
        
      labelMapRegistrationXMLData.image_type.assign("leftLungRightLung");
      labelMapRegistrationXMLData.transformationIndex = 0;
        
        
      //if the patient IDs are specified  as args, use them,
      //otherwise, extract from patient path

      int pathLength = 0, pos=0, next=0;

      if ( strcmp(movingImageID.c_str(), "q") != 0 )
	{
	  labelMapRegistrationXMLData.sourceID.assign(movingImageID);
	}
      else
	{
	  //first find length of path
	  next=1;
	  while(next>=1)
	    {
	      next = movingImageFileName.find("/", next+1);
	      pathLength++;
	    }
	  pos=0;
	  next=0;
       
	  std::string tempSourceID;
	  for (int i = 0; i < (pathLength-1);i++)
	    {
	      pos= next+1;
	      next = movingImageFileName.find("/", next+1);
	    }
       
	  labelMapRegistrationXMLData.sourceID.assign(movingImageFileName.c_str());
	  labelMapRegistrationXMLData.sourceID.erase(next,labelMapRegistrationXMLData.sourceID.length()-1);
	  labelMapRegistrationXMLData.sourceID.erase(0, pos);
	}
    
      if ( strcmp(fixedImageID.c_str(), "q") != 0 )
	{
	  labelMapRegistrationXMLData.destID =fixedImageID.c_str();
	}
      else
	{
	  pos=0;
	  next=0;
	  for (int i = 0; i < (pathLength-1);i++)
	    {
	      pos = next+1;
	      next = fixedImageFileName.find('/', next+1);
	    }

	  labelMapRegistrationXMLData.destID.assign(fixedImageFileName.c_str());//
	  //=tempSourceID.c_str();//movingImageFileName.substr(pos, next-1).c_str();

	  labelMapRegistrationXMLData.destID.erase(next,labelMapRegistrationXMLData.destID.length()-1);
	  labelMapRegistrationXMLData.destID.erase(0, pos);

	}

      //remove path from output transformation file before storing in xml
      std::cout<<"outputtransform filename="<<outputTransformFileName.c_str()<<std::endl;
      pos=0;
      next=0;
      for (int i = 0; i < (pathLength);i++)
        {
	  pos = next+1;
	  next = outputTransformFileName.find('/', next+1);
        }

      labelMapRegistrationXMLData.transformationLink.assign(outputTransformFileName.c_str());
      labelMapRegistrationXMLData.transformationLink.erase(0,pos);
      WriteRegistrationXML(infoFilename.c_str(),labelMapRegistrationXMLData);

      std::cout << "Writing transform..." << std::endl;
      itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
      transformWriter->SetInput( finalTransform2D );

   
      transformWriter->SetFileName( outputTransformFileName );
      transformWriter->Update();
    }

      std::cout << "DONE." << std::endl;

      return 0;
    }
