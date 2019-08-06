/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
 //  Software Guide : BeginLatex
 //
 //  Probably the most common representation of datasets in clinical
 //  applications is the one that uses sets of DICOM slices in order to compose
 //  3-dimensional images. This is the case for CT, MRI and PET scanners. It is
 //  very common therefore for image analysts to have to process volumetric
 //  images stored in a set of DICOM files belonging to a
 //  common DICOM series.
 //
 //  The following example illustrates how to use ITK functionalities in order
 //  to read a DICOM series into a volume and then save this volume in another
 //  file format.
 //
 //  The example begins by including the appropriate headers. In particular we
 //  will need the \doxygen{GDCMImageIO} object in order to have access to the
 //  capabilities of the GDCM library for reading DICOM files, and the
 //  \doxygen{GDCMSeriesFileNames} object for generating the lists of filenames
 //  identifying the slices of a common volumetric dataset.
 //
 //  \index{itk::ImageSeriesReader!header}
 //  \index{itk::GDCMImageIO!header}
 //  \index{itk::GDCMSeriesFileNames!header}
 //  \index{itk::ImageFileWriter!header}
 //
 //  Software Guide : EndLatex
 // Software Guide : BeginCodeSnippet
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h" 
#include "itkWarpImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkTransformFileReader.h"
#include "itkRegularStepGradientDescentBaseOptimizer.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkResampleImageFilter.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkImageRegistrationMethod.h"
#include "itkImageRegistrationMethodv4.h"
#include "itkBSplineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkRigid3DTransform.h"


template < class TPixel = double, unsigned int VImageDimension = 3>
class CommandIterationUpdate : public itk::Command
{
public:

	typedef  CommandIterationUpdate                         Self;
	typedef  itk::Command                                   Superclass;
	typedef  itk::SmartPointer<Self>                        Pointer;
	itkNewMacro(Self);


	using InternalPixelType = double;

	typedef itk::Image< TPixel, VImageDimension>           InternalImageType;
	typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
	typedef itk::Image<  VectorPixelType, VImageDimension > DeformationFieldType;
	using DisplacementFieldType = itk::Image< VectorPixelType, VImageDimension >;


public:
	CommandIterationUpdate() {};
	// Software Guide : EndCodeSnippet
public:

	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	using OptimizerPointer = const OptimizerType *;
	int iteration; 

	void Execute(itk::Object *caller, const itk::EventObject & event) override
	{
		Execute((const itk::Object *)caller, event);
	}

	void Execute(const itk::Object * object, const itk::EventObject & event) override
	{
    		std::cout << "*************"<< std::endl; 
			OptimizerType::ConstPointer optimizer = static_cast<const OptimizerType*>(object); 
			std::cout << optimizer->GetCurrentIteration() << "   ";
			std::cout << optimizer->GetValue() << "   ";
			std::cout << optimizer->GetCurrentPosition() << "   ";

	}

};



 // Software Guide : EndCodeSnippet
int main(void)
{

	// Define input and output image type 
	using PixelType = double;
	constexpr unsigned int Dimension = 3;
	const unsigned int SplineOrder = 3; 
	using FixedImageType = itk::Image< PixelType, Dimension >;
	using MovingImageType = itk::Image<PixelType, Dimension >;

	using ImageIOType = itk::GDCMImageIO;
	ImageIOType::Pointer dicomIO = ImageIOType::New();
	
	// Define fixed image reader
	using FixedReaderType = itk::ImageSeriesReader< FixedImageType >;
	FixedReaderType::Pointer fixedImageReader = FixedReaderType::New();
	fixedImageReader->SetImageIO(dicomIO); 

	// Define moving image reader
	using MovingReaderType = itk::ImageSeriesReader< MovingImageType >; 
	MovingReaderType::Pointer movingImageReader = MovingReaderType::New(); 
	movingImageReader->SetImageIO(dicomIO); 
	
	// Define dicom series name generator 
	using NamesGeneratorType = itk::GDCMSeriesFileNames;
	NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
	nameGenerator->SetUseSeriesDetails(true);
	nameGenerator->AddSeriesRestriction("0008|0021");

	NamesGeneratorType::Pointer fixedNameGenerator = NamesGeneratorType::New();
	NamesGeneratorType::Pointer movingNameGenerator = NamesGeneratorType::New();

	// fix directory
	std::string const fixedDirectoryName = "C:\\Registration\\breastcancer\\ha\\pre-prone\\dicom\\es";
	std::string const movingDirectoryName = "C:\\Registration\\breastcancer\\ha\\pre-supine\\dicom";

	fixedNameGenerator->SetDirectory(fixedDirectoryName);
	movingNameGenerator->SetDirectory(movingDirectoryName);

	FixedImageType::Pointer fixedVolume = FixedImageType::New(); 
	MovingImageType::Pointer movingVolume = FixedImageType::New(); 

	// read dicom volume
	try
	{
		using SeriesIdContainer = std::vector< std::string >;

		const SeriesIdContainer & fixedSeriesUID = fixedNameGenerator->GetSeriesUIDs();
		std::string fixedSeriesIdentifier = fixedSeriesUID.begin()->c_str();

		const SeriesIdContainer & movingSeriesUID = movingNameGenerator->GetSeriesUIDs();
		std::string movingSeriesIdentifier = movingSeriesUID.begin()->c_str();


		using FileNamesContainer = std::vector< std::string >;;

		FileNamesContainer fixedFileNames;
		FileNamesContainer movingFileNames;

		fixedFileNames = fixedNameGenerator->GetFileNames(fixedSeriesIdentifier);
		movingFileNames = movingNameGenerator->GetFileNames(movingSeriesIdentifier);

		fixedImageReader->SetFileNames(fixedFileNames);
		movingImageReader->SetFileNames(movingFileNames);

		try
		{
			fixedImageReader->Update();
			movingImageReader->Update();

		}
		catch (itk::ExceptionObject &ex)
		{
			std::cout << ex << std::endl;
			return EXIT_FAILURE;
		}

		fixedVolume = fixedImageReader->GetOutput();
		movingVolume = movingImageReader->GetOutput();

	}

	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return EXIT_FAILURE;
	}
	

	// Interpolators
	using FixedLinearInterpolatorType = itk::LinearInterpolateImageFunction<FixedImageType,double >;
	using MovingLinearInterpolatorType = itk::LinearInterpolateImageFunction<MovingImageType, double>; 

	// Metric
	//typedef itk::MattesMutualInformationImageToImageMetric<FixedImageType, MovingImageType> MetricType; 
	typedef itk::MeanSquaresImageToImageMetricv4<FixedImageType, MovingImageType> MetricType; 

	// Transforms
	typedef itk::BSplineTransform<double, Dimension, SplineOrder> TransformType; 
	
	// Optimizer
	using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
	
	// Registration
	typedef itk::ImageRegistrationMethod<FixedImageType, MovingImageType> RegistrationType; 
	typedef itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType> RegistrationTypev4; 

	// Instantiate Objects
	MetricType::Pointer metric = MetricType::New(); 
	OptimizerType::Pointer optimizer = OptimizerType::New(); 
	RegistrationTypev4::Pointer registrationv4 = RegistrationTypev4::New(); 

	// Create Interpolators for fixed and moving images
	FixedLinearInterpolatorType::Pointer fixedInterpolator = FixedLinearInterpolatorType::New(); 
	MovingLinearInterpolatorType::Pointer movingInterpolator = MovingLinearInterpolatorType::New();

	// Add interpolators to metric
	metric->SetFixedInterpolator(fixedInterpolator); 
	metric->SetMovingInterpolator(movingInterpolator); 


	// Create observer object to track registration 
	CommandIterationUpdate<PixelType, Dimension>::Pointer observer = CommandIterationUpdate<PixelType, Dimension>::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);


	// Add processors to registration object
	registrationv4->SetOptimizer(optimizer); 
	registrationv4->SetMetric(metric); 
	
	// Set Moving Initial Transform
	TransformType::Pointer movingInitialTransform = TransformType::New();
	TransformType::ParametersType initialParameters(
		movingInitialTransform->GetNumberOfParameters());
	initialParameters[0] = 0.0; // Initial offset in mm along X
	initialParameters[1] = 0.0; // Initial offset in mm along Y
	movingInitialTransform->SetParameters(initialParameters);
	registrationv4->SetMovingInitialTransform(movingInitialTransform);

	// Set Fixed and Moving Images
	registrationv4->SetFixedImage(fixedVolume); 
	registrationv4->SetMovingImage(movingVolume); 

	// set options for optimizer
	 optimizer->SetLearningRate( 4 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetRelaxationFactor( 0.5 );

	try
	{
		registrationv4->Update();
		std::cout << "Optimizer stop condition: "
			<< registrationv4->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	// Registration result is saved in transform 
	TransformType::ConstPointer transform = registrationv4->GetTransform(); 

	// Mapping is done with ResampleFilter 
	using ResampleFilterType = itk::ResampleImageFilter<MovingImageType, FixedImageType >;
	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(transform);
	resample->SetInput(movingVolume);
	resample->SetSize(fixedVolume->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedVolume->GetOrigin());
	resample->SetOutputSpacing(fixedVolume->GetSpacing());
	resample->SetOutputDirection(fixedVolume->GetDirection());

	try
	{
		resample->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}



	return EXIT_SUCCESS;
}