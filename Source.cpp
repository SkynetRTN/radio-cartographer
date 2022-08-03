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
	//WE NEED EXCEPTION HANDLING FOR FREQUENCIES/CONFIGURATIONS THAT DON'T ALIGN WITH THESE FREQUENCIES
	//RFI Max can't be used for multiple surveys with different dates yet.
	//RFIMax() needs exception handling for different configurations.

	double maxRFI;
	if (sParams.tele == FOURTY_FOOT)
	{
		if (sParams.frequency == 1.405)
		{
			if (sParams.channel == COMPOSITE) { maxRFI = 0.7; }
			else { maxRFI = 0.7; }
		}
	}
	else
	{
		if (sParams.frequency == 9.000)
		{
			if (sParams.channel == COMPOSITE) { maxRFI = 0.8; }
			else { maxRFI = 0.8; }
		}
		else if (sParams.frequency == 1.395 && MJD < 56870)
		{
			if (sParams.channel == COMPOSITE) { maxRFI = 0.8; }
			else { maxRFI = 0.7; }
		}
		else if (sParams.frequency == 1.550)
		{
			if (sParams.channel == COMPOSITE) { maxRFI = 0.9; }
			else { maxRFI = 0.8; }
		}
		else if ((sParams.frequency == 1.680 || sParams.frequency == 1.700) && MJD < 56870)
		{
			if (sParams.channel == COMPOSITE) { maxRFI = 1.1; }
			else { maxRFI = 0.9; }
		}
	}

	//if (sParams.simulation == true)
	//{
	//	maxRFI = 0.95;
	//}
	//else
	//{
	//	maxRFI = 0.8;
	//}
	maxRFI = 0.8;

	return maxRFI;
}

void setInputMapParams(char* argv[], MapParameters &mParams)
{
	// TEMPORARY
	double photometryOn = atof(argv[11]);
	double m10PlusCriteria = atof(argv[17]);
	double largeScaleStruct = atof(argv[10]);

	// PERMANENT
	mParams.SSSMapping = true;
	mParams.pixelSize = 0.05;
	mParams.rfiScale = atof(argv[6]);

	if (m10PlusCriteria == 1) {
		mParams.m10PlusProcessing = true;
	}
	else {
		mParams.m10PlusProcessing = false;
	}

	if (largeScaleStruct == 1) {
		mParams.LSSMapping = true;
	}
	else {
		mParams.LSSMapping = false;
	}

	if (photometryOn) {
		mParams.correlatedWeightMap = true;
	}
	else {
		mParams.correlatedWeightMap = false;
	}
}
void setInputSurveyParams(char* argv[], SurveyParameters &sParams)
{
	// TEMPORARY
	std::string inputChannel = argv[3];
	std::string inputCalibrationMethod = argv[4];
	std::string inputProcessingCoordinates = argv[2];

	// DEBUGGING
	sParams.forcedTS = 0.0;
	sParams.tele = TWENTY_METER;
	sParams.trimSize = atof(argv[16]);

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
	else
	{
		sParams.calMethod = INTERPOLATED;
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
void setInputPhotoParams(char* argv[], PhotoParams &pParams)
{
	double perform = atof(argv[11]);
	pParams.innerRadius = atof(argv[12]);
	pParams.outerRadius = atof(argv[13]);
	std::string centroidType = argv[14];

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

	//pParams.numberOfSources = atof(argv[12]);
	//pParams.coordinatesDeg = { 0, 0 };
	//pParams.coordinatesPixels = { 103, 21, 231, 61, 153, 161, 134, 182, 385, 196, 331, 214, 245, 258, 201, 267, 285, 360, 299, 409 };

}
void setInputSpectralParams(char* argv[], SpectralParameters &cParams)
{
	// IN-PROGRESS
	cParams.subScale = 0.0;
	cParams.modSubScale = 0.0;
	cParams.modSubZones = {};// {1435.0, 1447.5, 1522.5, 1630.0};
	cParams.receiver = XX;
	cParams.velocity = 0.0;

	//cParams.files.push_back(skynet_filename);

	/*cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.cyb.fits");
	cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.02.cyb.fits");
	cParams.files.push_back("Skynet_58397_fig16_2_36380_37171.03.cyb.fits");*/
}
void setInputProcessingParams(char* argv[], ProcessorParameters &procParams)
{
	double photometryOn = atof(argv[11]);
	double rawMap = atof(argv[9]);

	procParams.bgScaleBW = atof(argv[5]);
	procParams.rfiScaleBW = atof(argv[6]);
	procParams.timeShift = atof(argv[8]);
	procParams.wScaleBW = atof(argv[7]);

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

int main(int argc, char* argv[])
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
	arg 5 = Background Scale
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
  arg 18+ = all file names for appending
  -> just listed out
	*/

	//WELCOME TO RDP_2.0! INSTRUCTIONS IN COMMENTS BELOW
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
	setInputSpectralParams(argv, cParams);

	// PHOTOMETRY
	PhotoParams pParams;
	setInputPhotoParams(argv, pParams);

	// PROCESSING
	ProcessorParameters procParams;
	setInputProcessingParams(argv, procParams);

	// OUTPUT
	Output output;
	output.printTelescopeInfo(sParams);
	output.printWScale(mParams.processedWeightScale);

	// SURVEY
  cParams.files.push_back(skynet_filename);
  Survey survey(sParams, cParams, skynet_filename);
  surveyHold.push_back(survey);

  int i;
  for ( i = 18; i < argc; i++ ) { // can we make this less wack, eventually?
    // printf("survey: %s\n", argv[i]);
    cParams.files.clear();
	  cParams.files.push_back(argv[i]);
	  Survey surveyN(sParams, cParams, argv[i]);
	  surveyHold.push_back(surveyN);
  }

	// cParams.files.push_back("40963_oj287_jul_30_1230.fits");
	// Survey survey(sParams, cParams, "40963_oj287_jul_30_1230.fits");
	// surveyHold.push_back(survey);

	// cParams.files.clear();
	// cParams.files.push_back("40965_oj287_jul_30_1330.fits");
	// Survey survey2(sParams, cParams, "40965_oj287_jul_30_1330.fits");
	// surveyHold.push_back(survey2);

	// cParams.files.clear();
	// cParams.files.push_back("40966_oj287_jul_30_1430.fits");
	// Survey survey3(sParams, cParams, "40966_oj287_jul_30_1430.fits");
	// surveyHold.push_back(survey3);

	// cParams.files.clear();
	// cParams.files.push_back("40967_oj287_jul_30_1530.fits");
	// Survey survey4(sParams, cParams, "40967_oj287_jul_30_1530.fits");
	// surveyHold.push_back(survey4);

	// Survey survey2(sParams, "Skynet_56812_Jupiter_8889_9652.txt");
	// surveyHold.push_back(survey2);
	// "40963_oj287_jul_30_1230.fits"

	rfiMax = setRFIMaxVal(sParams, (int)survey.getMJD());

	// PROCESSOR
	Processor processor(procParams);

	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.determineDataBoundaries(surveyHold[i]);
	}

	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.characterizeData(surveyHold[i]);
		processor.performBGSubtractionMulti(surveyHold[i]);
		processor.performTimeShifting(surveyHold[i]);
		processor.set2DScatter(surveyHold[i]);
		processor.processLargeScaleStructure(surveyHold[i]);
		surveyHold[i].calculateEdgeParameters();
	}

	if (mParams.LSSMapping)
	{
		processor.levelLSSData(surveyHold);
	}
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

		procMap.printSSSCorrelationMap();
	}
	if (!(pParams.perform) || procParams.raw)
	{
		output.printPhotometryHolder();
	}
	if (!(procParams.raw))
	{
		procMap.printSSSScaleMap();
		procMap.printSSSWeightMap();
	}

	// EXIT
	std::cout << "The code has exited successfully with return code 0" << std::endl;
	
	return 0;
}