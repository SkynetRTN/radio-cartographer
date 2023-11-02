#include <stdexcept>
#include <iostream>
#include <fstream>


#include "PreProcessorNew.h"
#include "Processor.h"
#include "Cartographer.h"
#include "Composite.h"
#include "OutputFile.h"
#include "Analysis.h"
#include "FileReader.h"
#include "io\config\ConfigParser.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

bool LSSProcessing = false;
std::string skynet_filename;

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
	std::string inputCalibrationMethod = argv[3];
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
void setInputPhotoParams(char *argv[], PhotometryParams &pParams)
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

	procParams.lssProc = false;
}

int main(int argc, char *argv[])
{	
	//=================================================================
	//================== PARSE THE INPUT ARGUMENTS ====================
	//=================================================================
	skynet_filename = argv[1];

	MapParameters mParams;
	setInputMapParams(argv, mParams);

	SurveyParameters sParams;
	setInputSurveyParams(argv, sParams);

	SpectralParameters cParams;
	setInputSpectralParams(argc, argv, cParams);

	PhotometryParams pParams;
	setInputPhotoParams(argv, pParams);

	ProcessorParameters procParams;
	setInputProcessingParams(argv, procParams);

	//=================================================================
	//====================== READ IN DATA FILE ========================
	//=================================================================

	ConfigParser configuration("example.json");
	configuration.parse();

	FileReader reader(skynet_filename);
	Input input = reader.read();

	std::vector<std::string> filenamesX = {"x", "y", "z"};
	std::vector<FileReader> x = {};

	for (int i = 0; i < 3; ++i) {
		x.emplace_back(FileReader(filenamesX[i]));
	}

	/*FileReader reader2("G28_C02.raw.dcr.fits");
	Input input2 = reader2.read();

	FileReader reader3("G28_C03.raw.dcr.fits");
	Input input3 = reader3.read();*/

	/*FileReader reader20M("0088749.fits");
	Input input20M = reader20M.read();*/

	//=================================================================
	//===================== PREPROCESS THE DATA =======================
	//=================================================================
	if (reader.getFileType() == FileReader::FileType::FITS_20M) {
		PreProcessingParameters ppParams;
		//ppParams.inclusionBands.emplace_back(std::make_pair(1355.0 * pow(10.0, 6.0), 1400.0* pow(10.0, 6.0)));
		//ppParams.modifiedSubtractionBands.emplace_back(std::make_pair(1300.0 * pow(10.0, 6.0), 1500.0 * pow(10.0, 6.0)));
		//ppParams.subtractionScale = 6.0;

		PreProcessorNew preProcessor = PreProcessorNew();
		preProcessor.process(ppParams, input);
	}

	//=================================================================
	//====================== CREATE THE SURVEY ========================
	//=================================================================

	std::vector<Survey> surveyHold, surveys;
	Survey survey(sParams, input);
	surveyHold.push_back(survey);

	/*Survey survey2(sParams, input2);
	surveyHold.push_back(survey2);

	Survey survey3(sParams, input3);
	surveyHold.push_back(survey3);*/

	//Survey survey2(sParams, cParams, skynet_filename);
	//surveyHold.push_back(survey2);

	// OUTPUT
	Output output;
	output.printTelescopeInfo(sParams);
	output.printWScale(mParams.processedWeightScale);

	//=================================================================
	//====================== PROCESS THE DATA =========================
	//=================================================================
	Processor processor(procParams);

	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.determineDataBoundaries(surveyHold[i]);
	}

	// BACKGROUND SUBTRACTION
	for (int i = 0; i < surveyHold.size(); i++)
	{
		processor.characterizeData(surveyHold[i]);
		processor.performBGSubtractionMulti(surveyHold[i]);
		processor.performTimeShifting(surveyHold[i]);
		processor.set2DScatter(surveyHold[i]);
		processor.processLargeScaleStructure(surveyHold[i]);
		surveyHold[i].calculateEdgeParameters();
	}

	//processor.levelLSSData(surveyHold);
	surveys = surveyHold;
	Composite composite(surveys);

	// RFI AND THETA GAP
	processor.performRFIRejectionMulti(composite, surveys);
	processor.calculateProcThetaGapMulti(composite);

	// PROCESSED MAP
	Cartographer cartographer(composite, mParams);
	Map procMap = cartographer.generateProcMapsMulti();
	procMap.printSSSProcFluxMap();
	procMap.printProcPathMap();

	//=================================================================
	//===================== PERFORM PHOTOMETRY ========================
	//=================================================================
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

	return EXIT_SUCCESS;
}