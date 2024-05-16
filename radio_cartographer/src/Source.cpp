#include <stdexcept>
#include <iostream>
#include <fstream>

#include "io\OutputFile.h"
#include "PreProcessorNew.h"
#include "processor\Processor.h""
#include "map\Cartographer.h"
#include "map\Composite.h"
#include "io\observation\FileReader.h"
#include "photometry\photometry.h"
#include "io\config\ConfigParser.h"


/*
 * Skynet's Single-Dish Radio Mapping Algorithm. 
 * 
 * This algorithm's implementation has been published in the Astrophysical 
 * Journal and can be read here for free: https://arxiv.org/abs/1808.06128.
 * The algorithm has been integrated into the image-processing library of 
 * the Skynet Robotic Telescope Network, which includes optical telescopes
 * spanning four continents, and now also Green Bank Observatory's 20-meter 
 * diameter radio telescope in West Virginia. Skynet serves hundreds of
 * professional users, and additionally tens of thousands of students, of 
 * all ages.
 * 
 * The Single-Dish algorithm makes heavy use of Skynet's outlier rejection
 * rejection algorithm, Robust Chauvnet Rejection (RCR). The RCR algorithm
 * has been published in the Astrophysical Journal and is freely available
 * for reading at https://arxiv.org/abs/1807.05276. An additional paper can 
 * also be read for free at https://arxiv.org/abs/2301.07838.
 *
 * Source code developed by Dylan Dutton, John Martin, Michael Maples, and 
 * Travis Berger. Front-end development and integration into the Skynet 
 * Network by Reed Fu, Omar H. Shaban, and Josh Haislip. Algorithm design
 * by Dr. Dan Reichart. Additional support by James Finn and Nick Konz.
 *
 */

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
void setInputSpectralParams(int argc, char *argv[], SpectralParameters &cParams)
{
	// IN-PROGRESS
	cParams.subScale = 6.0;
	cParams.modSubScale = 0.0;
	cParams.modSubZones = {}; 
	cParams.receiver = HI;
	if (atoi(argv[3]) == 2)
		cParams.receiver = LO;

	cParams.velocity = 0.0;

	// FREQUENCY SELECTION
	double minFreq = atof(argv[9]);
	double maxFreq = atof(argv[10]);
	cParams.inclusionBand = { minFreq, maxFreq }; // MHz
	cParams.exclusionBand = {};
	for (int i = 22; i < argc; i++)
	{
		cParams.exclusionBand.push_back(atof(argv[i]));
	}

	cParams.files.push_back("0106181.fits");
}

int main(int argc, char *argv[])
{	
	// ================ READ IN THE CONFIGURATION FILE ================
	ConfigParser config = ConfigParser("example.json");
	config.parse();

	//// ================= READ IN THE OBSERVATION DATA =================
	std::vector<Observation> observations;

	// testing
	/*SpectralParameters cParams;
	setInputSpectralParams(argc, argv, cParams);
	SurveyParameters sParams;
	setInputSurveyParams(argv, sParams);
	Survey surveyOld = Survey(sParams, cParams, "0117690.fits");*/
	//surveys.emplace_back(surveyOld);

	for (const std::string& filename : config.obsFilenames)
	{
		FileReader file = FileReader(filename);
		Observation obs = file.read();

		// If file contained multiple maps, split them up into individual observations
		if (file.getFileType() == FileReader::FileType::FITS_GBT && obs.hasMultipleMaps)
		{
			for (int i = 1; i < obs.mapStartLocations.size(); ++i)
			{
				int start = obs.mapStartLocations[i - 1];
				int end = obs.mapStartLocations[i];

				observations.emplace_back(obs.splitObservation(start, end));

				if (i == obs.mapStartLocations.size() - 1) 
				{
					observations.emplace_back(obs.splitObservation(end));
				}
			}
		}
		else 
		{
			if (file.getFileType() == FileReader::FileType::FITS_20M)
			{
				PreProcessorNew preProcessor(config.preProcessor);
				preProcessor.process(obs);
			}

			observations.emplace_back(obs);
		}
	}

	if (config.runMethod == "merge") {
		observations = { FileReader::merge(observations) };
	}

	//// ================== TRANSFORM TO SURVEY DATA ====================
	std::vector<Survey> surveys;
	//surveys.emplace_back(surveyOld);

	for (const Observation& obs : observations ){
		surveys.emplace_back(Survey(config.survey, obs));
	}

	Output output;
	output.printTelescopeInfo(config.survey);
	output.printWScale(config.map.processedWeightScale);

	//// ===================== PROCESS THE SURVEYS ======================
	Processor processor(config.processor);

	for (Survey& survey : surveys) {
		processor.determineDataBoundaries(survey);
	}

	for (Survey& survey : surveys)
	{
		processor.characterizeData(survey);
		processor.performBGSubtractionMulti(survey);
		processor.performTimeShifting(survey);
		processor.set2DScatter(survey);
		processor.processLargeScaleStructure(survey);
		survey.calculateEdgeParameters();
	}

	Composite composite(surveys);

	processor.performRFIRejectionMulti(composite, surveys);
	processor.calculateProcThetaGapMulti(composite, config.map.LSSMapping);

	// ====================== GENERATE THE MAPS =======================
	Cartographer cartographer(composite, config.map);
	Map processedMap = cartographer.generateProcMapsMulti();

	processedMap.printSSSProcFluxMap();
	processedMap.printProcPathMap();

	// ====================== PERFORM PHOTOMETRY ======================
	if (config.photometry.perform && !config.processor.raw)
	{
		Photometry photometry(composite, config.map);
		photometry.photometerMulti(processedMap, config.photometry);
	}

	// =============== GENERATE THE REMAINING MAPS ====================
	if (!config.photometry.perform || config.processor.raw)
	{
		output.printPhotometryHolder();
	}

	processedMap.printSSSScaleMap();
	processedMap.printSSSWeightMap();
	processedMap.printSSSCorrelationMap();

	return EXIT_SUCCESS;
}