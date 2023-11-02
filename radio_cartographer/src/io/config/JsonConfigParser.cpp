#include <fstream>
#include <iostream>
#include "io\config\JsonConfigParser.h"


using json = nlohmann::json;

JsonConfigParser::JsonConfigParser(ConfigParser& config) {
	JsonConfigParser(config, false);
}

JsonConfigParser::JsonConfigParser(ConfigParser& config, bool logToConsole) {
	std::ifstream f(config.filename);
	this->data = json::parse(f);
	this->configuration = &config;
	this->logToConsole = logToConsole;
}

void JsonConfigParser::parse() {
	parseMapParameters();
	parseSurveyParameters();
	parsePhotometryParameters();
	parseProcessingParameters();
	parsePreProcessingParameters();

	if (logToConsole)
		std::cout << "[ SUCCESS ] CONFIGURATION" << std::endl;
}

/**
 * Parses the configuration parameters in the 'map' element.
 */
void JsonConfigParser::parseMapParameters() {
	configuration->map.pixelSize = data["map"]["pixelSize"].template get<double>();
	configuration->map.m10PlusProcessing = data["map"]["m10PlusCriteria"].template get<bool>();
	configuration->map.LSSMapping = data["map"]["largeScaleMapping"].template get<bool>();
}

/**
 * Parses the configuration parameters in the 'survey' element.
 */
void JsonConfigParser::parseSurveyParameters() {
	std::string fluxChannel = data["survey"]["fluxChannel"].template get<std::string>();
	configuration->survey.channel = configuration->enumFromString<Channel>(fluxChannel);

	std::string calMethod = data["survey"]["calibrationChannel"].template get<std::string>();
	configuration->survey.calMethod = configuration->enumFromString<CalMethods>(calMethod);

	std::string coordinate = data["survey"]["processingCoordinate"].template get<std::string>();
	configuration->survey.pCoordinate = configuration->enumFromString<Coordinates>(coordinate);
}

/**
 * Parses the configuration parameters in the 'photometry' element.
 */
void JsonConfigParser::parsePhotometryParameters() {
	configuration->photomtery.perform = data["photometry"]["perform"].template get<bool>();
	configuration->photomtery.innerRadius = data["photometry"]["aperture"].template get<double>();;
	configuration->photomtery.outerRadius = data["photometry"]["annulus"].template get<double>();;
	
	std::string centroid = data["photometry"]["centroid"].template get<std::string>();
	configuration->photomtery.outerRadius = configuration->enumFromString<CentroidMethods>(centroid);
}

/**
 * Parses the configuration parameters in the 'processing' element. This includes 
 * parameters for Background Subtraction, RFI Subtraction, and Time Shifting.
 */
void JsonConfigParser::parseProcessingParameters() {
	configuration->processor.raw = configuration->intFromBool(data["processor"]["raw"].template get<bool>());
	configuration->processor.wScaleBW = data["processor"]["weightScale"].template get<double>();
	
	//configuration->processor.performBG = data["processor"]["background"]["perform"].template get<bool>();
	configuration->processor.bgScaleBW = data["processor"]["background"]["scale"].template get<double>();

	//configuration->processor.performRFI = data["processor"]["rfi"]["perform"].template get<bool>();
	configuration->processor.rfiScaleBW = data["processor"]["rfi"]["scale"].template get<double>();

	//configuration->processor.performTS = data["processor"]["timeShift"]["perform"].template get<bool>();
	configuration->processor.timeShift = configuration->intFromBool(data["processor"]["timeShift"]["perform"].template get<bool>());
	configuration->processor.timeShiftValue = data["processor"]["timeShift"]["strength"].template get<double>();
}

/**
 * Parses the configuration parameters in the 'preprocessing' element.
 */
void JsonConfigParser::parsePreProcessingParameters() {
	//configuration->preProcessor.perform = data["preProcessor"]["perform"].template get<bool>();
	configuration->preProcessor.subtractionScale = data["preProcessor"]["scale"].template get<double>();
	configuration->preProcessor.modifiedSubtractionScale = data["preProcessor"]["modifiedScale"].template get<double>();
	configuration->preProcessor.velocity = data["preProcessor"]["velocity"].template get<double>();
	//configuration->preProcessor. = data["preProcessor"]["minFrequency"].template get<double>();
	
	
	
	
	configuration->preProcessor.velocity = data["preProcessor"]["velocity"].template get<double>();

}


JsonConfigParser::~JsonConfigParser() {

}