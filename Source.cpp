#include "Cartographer.h"
#include "Processor.h"
#include "Composite.h"
#include "OutputFile.h"
#include "Analysis.h"
#include <stdexcept>
#include <iostream>

double rfiMax;
bool LSSProcessing;
std::string skynet_filename;

double setRFIMaxVal(SurveyParameters sParams, int MJD)
{
	// WE NEED EXCEPTION HANDLING FOR FREQUENCIES/CONFIGURATIONS THAT DON'T ALIGN WITH THESE FREQUENCIES
	// RFI Max can't be used for multiple surveys with different dates yet.
	// RFIMax() needs exception handling for different configurations.

	double maxRFI;
	if (sParams.tele == FOURTY_FOOT)
	{
		if (sParams.frequency == 1.405)
		{
			if (sParams.channel == COMPOSITE)
			{
				maxRFI = 0.7;
			}
			else
			{
				maxRFI = 0.7;
			}
		}
	}
	else
	{
		if (sParams.frequency == 9.000)
		{
			if (sParams.channel == COMPOSITE)
			{
				maxRFI = 0.8;
			}
			else
			{
				maxRFI = 0.8;
			}
		}
		else if (sParams.frequency == 1.395 && MJD < 56870)
		{
			if (sParams.channel == COMPOSITE)
			{
				maxRFI = 0.8;
			}
			else
			{
				maxRFI = 0.7;
			}
		}
		else if (sParams.frequency == 1.550)
		{
			if (sParams.channel == COMPOSITE)
			{
				maxRFI = 0.9;
			}
			else
			{
				maxRFI = 0.8;
			}
		}
		else if ((sParams.frequency == 1.680 || sParams.frequency == 1.700) && MJD < 56870)
		{
			if (sParams.channel == COMPOSITE)
			{
				maxRFI = 1.1;
			}
			else
			{
				maxRFI = 0.9;
			}
		}
	}

	// if (sParams.simulation == true)
	//{
	//	maxRFI = 0.95;
	// }
	// else
	//{
	//	maxRFI = 0.8;
	// }
	maxRFI = 0.8;

	return maxRFI;
}

void setInputMapParams(char *argv[], MapParameters &mParams)
{
	// TEMPORARY
	double photometryOn = atof(argv[15]);
	double m10PlusCriteria = atof(argv[13]);
	double largeScaleStruct = atof(argv[21]);

	// PERMANENT
	mParams.SSSMapping = true;
	mParams.pixelSize = 0.05;
	mParams.rfiScale = atof(argv[12]);

	// WHY DO WE HAVE MULTIPLE VARIABLES FOR THE SAME VALUE ?
	mParams.processedWeightScale = atof(argv[14]);

	if (m10PlusCriteria == 1)
	{
		mParams.m10PlusProcessing = true;
	}
	else
	{
		mParams.m10PlusProcessing = false;
	}

	if (largeScaleStruct == 1)
	{
		mParams.LSSMapping = true;
	}
	else
	{
		mParams.LSSMapping = false;
	}

	if (photometryOn)
	{
		mParams.correlatedWeightMap = true;
	}
	else
	{
		mParams.correlatedWeightMap = false;
	}
}
void setInputSurveyParams(char *argv[], SurveyParameters &sParams)
{
	// TEMPORARY
	std::string inputChannel = argv[2];
	std::string inputCalibrationMethod = argv[4];
	std::string inputProcessingCoordinates = argv[5];

	// DEBUGGING
	sParams.forcedTS = 0.0;
	sParams.tele = TWENTY_METER;
	sParams.trimSize = atof(argv[20]);

	// MISC
	sParams.tracking = false;

	// SET FLUX CHANNEL
	if (inputChannel == "left")
	{
		sParams.channel = LEFT;
	}
	else if (inputChannel == "right")
	{
		sParams.channel = RIGHT;
	}
	else
	{
		sParams.channel = COMPOSITE;
	}
	// SET CALIBRATION METHOD
	if (inputCalibrationMethod == "pre")
	{
		sParams.calMethod = PRE;
	}
	else if (inputCalibrationMethod == "post")
	{
		sParams.calMethod = POST;
	}
	else if (inputCalibrationMethod == "interpolated")
	{
		sParams.calMethod = INTERPOLATED;
	}
	else
	{
		sParams.calMethod = NONE;
	}
	// SET PROCESSING COORDINATES
	if (inputProcessingCoordinates == "equatorial")
	{
		sParams.pCoordinate = EQUATORIAL;
	}
	else if (inputProcessingCoordinates == "galactic")
	{
		sParams.pCoordinate = GALACTIC;
	}
	else
	{
		sParams.pCoordinate = EQUATORIAL;
	}
}
void setInputPhotoParams(char *argv[], PhotoParams &pParams)
{
	double perform = atof(argv[15]);
	pParams.innerRadius = atof(argv[16]);
	pParams.outerRadius = atof(argv[17]);
	std::string centroidType = argv[18];

	if (perform)
	{
		pParams.perform = true;
	}
	else
	{
		pParams.perform = false;
	}
	if (centroidType == "center")
	{
		pParams.centroidType = CENTER;
	}
	else if (centroidType == "brightest")
	{
		pParams.centroidType = BRIGHTEST;
	}
	else
	{
		pParams.centroidType = COORDINATES;
	}

	// pParams.numberOfSources = atof(argv[12]);
	// pParams.coordinatesDeg = { 0, 0 };
	// pParams.coordinatesPixels = { 103, 21, 231, 61, 153, 161, 134, 182, 385, 196, 331, 214, 245, 258, 201, 267, 285, 360, 299, 409 };
}
void setInputSpectralParams(int argc, char *argv[], SpectralParameters &cParams)
{
	// IN-PROGRESS
	cParams.subScale = 0.0;
	cParams.modSubScale = 0.0;
	cParams.modSubZones = {}; //{1423.929214, 1424.051278, 1424.066536, 1424.219116, 1424.234374, 1424.310664};// {1435.0, 1447.5, 1522.5, 1630.0};
	cParams.receiver = HI;
	if (atoi(argv[3]) == 2)
		cParams.receiver = LO;

	cParams.velocity = 0.0;

	// FREQUENCY SELECTION
	double minFreq = atof(argv[9]);
	double maxFreq = atof(argv[10]);
	cParams.inclusionBand = {minFreq, maxFreq}; // MHz
	cParams.exclusionBand = {};
	for (int i = 22; i < argc; i++)
	{
		cParams.exclusionBand.push_back(atof(argv[i]));
	}

	// FILENAME(S)
	cParams.files.push_back(skynet_filename);

	/*cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.cyb.fits");
	cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.02.cyb.fits");
	cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.03.cyb.fits");*/
}
void setInputProcessingParams(char *argv[], ProcessorParameters &procParams)
{
	double photometryOn = atof(argv[15]);
	double rawMap = atof(argv[8]);

	procParams.bgScaleBW = atof(argv[11]);
	procParams.rfiScaleBW = atof(argv[12]);
	procParams.timeShift = argv[6];
	procParams.timeShiftValue = atof(argv[7]);
	procParams.wScaleBW = atof(argv[14]);

	if (photometryOn)
	{
		procParams.wCorrMap = true;
	}
	else
	{
		procParams.wCorrMap = false;
	}
	if (rawMap)
	{
		procParams.raw = true;
	}
	else
	{
		procParams.raw = false;
	}

	// HARD-CODED FOR NOW
	procParams.lssProc = false;
}

int main(int argc, char *argv[])
{
	/*
	The inputs are as follows:
	arg 1 = filename
	-> The file to be processed
	arg 2 = Processing Coordinates
	-> These are the coordinates that we process the data in
	arg 3 = Flux Channel
	-> Processing channel
	arg 4 = Calibration Method
	-> Which calibration data to be used
	arg 5 = Background Scalecoo
	-> In beamwidths
	arg 6 = RFI Scale
	-> In beamwidths
	arg 7 = Weight Scale
	-> ...
	arg 8 = Time Shifting
	-> 1 or 0 corresponds to true/on or false/off
	arg 9 = Raw
	-> Generate a raw map
	arg 10 = LSS
	-> Interested in large
	arg 11 = photometry on
	-> ...
	arg 12 = Aperture Radius
	-> The radius around the source
	arg 13 = Annulus Radius
	-> The outer radius for background measurements
	arg 14 = Centroid Finding Method
	-> Look for center, brightest, or coordinates
	arg 15 = Coordinates
	-> Location of sources
	arg 16 = Trim Size
	-> Trims the turning edges
	arg 17 = m10+ Criteria
	-> 1 or 0 corresponds to true/on or false/off
	*/

	// WELCOME TO RDP_2.0! INSTRUCTIONS IN COMMENTS BELOW
	std::vector<Survey> surveyHold;
	std::vector<Survey> surveys;
	std::ofstream outputFile;

	// FILENAME
	skynet_filename = argv[1];

	// MAP
	MapParameters mParams;
	setInputMapParams(argv, mParams);

	// SURVEY
	SurveyParameters sParams;
	setInputSurveyParams(argv, sParams);

	// SPECTRAL
	SpectralParameters cParams;
	setInputSpectralParams(argc, argv, cParams);

	// PHOTOMETRY
	PhotoParams pParams;
	setInputPhotoParams(argv, pParams);

	// PROCESSING
	ProcessorParameters procParams;
	setInputProcessingParams(argv, procParams);

	// SURVEY
	Survey survey(sParams, cParams, skynet_filename);
//	Survey survey(sParams, cParams, "/skynet/radio-cartographer/testing/Temporary/barnard2.fits");
	surveyHold.push_back(survey);

//    cParams.files.clear();
//    cParams.files.emplace_back("/skynet/radio-cartographer/testing/Temporary/barnard.fits");
//    Survey survey2(sParams, cParams, "/skynet/radio-cartographer/testing/Temporary/barnard.fits");
//    surveyHold.push_back(survey2);

	// OUTPUT
	Output output;
	output.printTelescopeInfo(sParams);
	output.printWScale(mParams.processedWeightScale);

	std::cout << "Out of Survey\n";
	rfiMax = setRFIMaxVal(sParams, (int)survey.getMJD());

	// PROCESSOR
	Processor processor(procParams);

	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.determineDataBoundaries(surveyHold[i]);
	}
	std::cout << "determined Data bapoundarie\n";
	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.characterizeData(surveyHold[i]);
		processor.performBGSubtractionMulti(surveyHold[i]);
		processor.performTimeShifting(surveyHold[i]);
		processor.set2DScatter(surveyHold[i]);
		processor.processLargeScaleStructure(surveyHold[i]);
		surveyHold[i].calculateEdgeParameters();
	}

	processor.levelLSSData(surveyHold);
	surveys = surveyHold;
	Composite composite(surveys);

	// RFI AND THETA GAP
	processor.performRFIRejectionMulti(composite, surveys);
	processor.calculateProcThetaGapMulti(composite);

	// PROCMAP
	Cartographer cartographer(composite, mParams);
	Map procMap = cartographer.generateProcMapsMulti();
	procMap.printSSSProcFluxMap();
	procMap.printProcPathMap();

	// PHOTOMETRY AND OUTPUT
	if (pParams.perform && !(procParams.raw))
	{
		Analysis analysis(composite, mParams);
		analysis.photometryMulti(procMap, pParams);

		// procMap.printSSSCorrelationMap();
	}
	if (!(pParams.perform) || procParams.raw)
	{
		output.printPhotometryHolder();
	}
	// if (!(procParams.raw))
	// {
	procMap.printSSSScaleMap();
	procMap.printSSSWeightMap();
	procMap.printSSSCorrelationMap();
	// }

	// EXIT
	std::cout << "The code has exited successfully with return code 0" << std::endl;

	return 0;
}