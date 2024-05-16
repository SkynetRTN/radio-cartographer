#include <fstream>
#include <iostream>
#include "io\config\JsonConfigParser.h"


using json = nlohmann::json;

JsonConfigParser::JsonConfigParser(ConfigParser& config) : JsonConfigParser(config, false) {}

JsonConfigParser::JsonConfigParser(ConfigParser& config, bool logToConsole) {
	this->jsonFilename = config.configFilename;
	this->configuration = &config;
	this->logToConsole = logToConsole;
}

void JsonConfigParser::parse() {

	parseFile();
	parseMetaData();
	parseMapParameters();
	parseSurveyParameters();
	parsePhotometryParameters();
	parseProcessingParameters();
	parsePreProcessingParameters();
	parseScatterParameters();

	setEntangledParameters();

	if (logToConsole)
		std::cout << "[ SUCCESS ] JSON CONFIGURATION" << std::endl;
}


/**
 * Parses a JSON configuration file.
 */
void JsonConfigParser::parseFile() {
	std::ifstream f(jsonFilename);
	json j = json::parse(f, nullptr, false);

	if (j.is_discarded()) {
		throw std::runtime_error("[ Json Config Parser ] Invalid JSON format for " + jsonFilename);
	}
	else {
		this->data = j;
	}
}

/**
 * Parses the configuration parameters in the 'meta' element.
 *
 * Meta is a catch all element that contains information that is not easily
 * organized into it's own element.
 */
void JsonConfigParser::parseMetaData() {
	try { // try to parse a vector of filenames
		configuration->obsFilenames = data["meta"]["filename"].template get<std::vector<std::string>>();
	}
	catch (const json::type_error& e1) {
		try { // no vector so try to parse a single string
			configuration->obsFilenames = { data["meta"]["filename"].template get<std::string>() };
		}
		catch (const json::type_error& e2) {
			printErrorMessage("filename", "meta");
			throw;
		}
	}
	try { // parse the filename
		configuration->runMethod = data["meta"]["runMethod"].template get<std::string>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("runMethod", "meta");
		throw;
	}

	if (configuration->obsFilenames.empty()) {
		throw std::invalid_argument("[ Json Config Parser ] No input filenames provided.");
	}
	else if (configuration->obsFilenames.size() > 1 && configuration->runMethod == "single") {
		throw std::invalid_argument("[ Json Config Parser ] Cannot run 'single' on multiple observations.");
	}
}


/**
 * Parses the configuration parameters in the 'map' element.
 */
void JsonConfigParser::parseMapParameters() {
	try { // Read pixelSize element
		configuration->map.pixelSize = data["map"]["pixelSize"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("pixelSize", "map");
		throw;
	}
	try { // Read m10PlusCriteria element
		configuration->map.m10PlusProcessing = data["map"]["m10PlusCriteria"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("m10PlusCriteria", "map");
		throw;
	}
	try { // Read largeScaleMapping element
		configuration->map.LSSMapping = data["map"]["largeScaleMapping"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("largeScaleMapping", "map");
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'survey' element.
 */
void JsonConfigParser::parseSurveyParameters() {
	try { // Read fluxChannel element
		std::string fluxChannel = data["survey"]["fluxChannel"].template get<std::string>();
		configuration->survey.channel = configuration->enumFromString<Channel>(fluxChannel);
	}
	catch (const json::type_error& e) {
		printErrorMessage("fluxChannel", "survey");
		throw;
	}
	try { // Read calibrationChannel element
		std::string calMethod = data["survey"]["calibrationChannel"].template get<std::string>();
		configuration->survey.calMethod = configuration->enumFromString<CalMethods>(calMethod);
	}
	catch (const json::type_error& e) {
		printErrorMessage("calibrationChannel", "survey");
		throw;
	}
	try { // Read processingCoordinate element
		std::string coordinate = data["survey"]["processingCoordinate"].template get<std::string>();
		configuration->survey.pCoordinate = configuration->enumFromString<Coordinates>(coordinate);
	}
	catch (const json::type_error& e) {
		printErrorMessage("processingCoordinate", "survey");
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'photometry' element.
 */
void JsonConfigParser::parsePhotometryParameters() {
	try { // Read perform element
		configuration->photometry.perform = data["photometry"]["perform"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("perform", "photometry");
		throw;
	}
	try { // Read aperture element
		configuration->photometry.innerRadius = data["photometry"]["aperture"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("aperture", "photometry");
		throw;
	}
	try { // Read annulus element
		configuration->photometry.outerRadius = data["photometry"]["annulus"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("annulus", "photometry");
		throw;
	}
	try { // Read centroid element
		std::string centroid = data["photometry"]["centroid"].template get<std::string>();
		configuration->photometry.centroidType = configuration->enumFromString<CentroidMethods>(centroid);
	}
	catch (const json::type_error& e) {
		printErrorMessage("centroid", "photometry");
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'processing' element.
 *
 * This includes parameters for Background Subtraction, RFI Subtraction,
 * and Time Shifting.
 */
void JsonConfigParser::parseProcessingParameters() {

	try {
		configuration->processor.raw = configuration->intFromBool(data["processor"]["raw"].template get<bool>());
	}
	catch (const json::type_error& e) {
		printErrorMessage("raw", "processor");
		throw;
	}
	try {
		configuration->processor.wScaleBW = data["processor"]["weightScale"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("weightScale", "processor");
		throw;
	}

	parseRFIParameters();
	parseBackgroundParameters();
	parseTimeShiftingParameters();
}

/**
 * Parses the configuration parameters in the 'processor/background' element.
 */
void JsonConfigParser::parseBackgroundParameters() {
	try {
		configuration->processor.performBGS = data[kProcesser][kBackground]["perform"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("perform", kProcesser + "/" + kBackground);
		throw;
	}
	try {
		configuration->processor.bgScaleBW = data[kProcesser][kBackground]["scale"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("scale", kProcesser + "/" + kBackground);
		throw;
	}
	try {
		configuration->processor.maxBGThreads = data[kProcesser][kBackground]["maxThreads"].template get<int>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("maxThreads", kProcesser + "/" + kBackground);
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'processor/rfi' element.
 */
void JsonConfigParser::parseRFIParameters() {
	try {
		configuration->processor.performRFI = data[kProcesser][kRfi]["perform"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("perform", kProcesser + "/" + kRfi);
		throw;
	}
	try {
		configuration->processor.rfiScaleBW = data[kProcesser][kRfi]["scale"].template get<double>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("scale", kProcesser + "/" + kRfi);
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'processor/timeShift' element.
 */
void JsonConfigParser::parseTimeShiftingParameters() {
	try {
		configuration->processor.performTS = data[kProcesser][kTimeShift]["perform"].template get<bool>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("perform", kProcesser + "/" + kTimeShift);
		throw;
	}
	try {
		configuration->processor.timeShiftValue = getOptionalDouble(data[kProcesser][kTimeShift]["strength"]);
	}
	catch (const json::type_error& e) {
		printErrorMessage("strength", kProcesser + "/" + kTimeShift);
		throw;
	}
}

/**
 * Parses the configuration parameters in the 'preProcessing' element.
 */
void JsonConfigParser::parsePreProcessingParameters() {

	// @TODO: Add exception handling and clean up

	configuration->preProcessor.perform = data[kPreProcesser]["perform"].template get<bool>();
	configuration->preProcessor.subtractionScale = data[kPreProcesser]["scale"].template get<double>();
	configuration->preProcessor.modifiedSubtractionScale = data[kPreProcesser]["modifiedScale"].template get<double>();

	std::string receiver = data[kPreProcesser]["receiver"].template get<std::string>();
	configuration->preProcessor.receiver = configuration->enumFromString<Receiver2>(receiver);

	configuration->preProcessor.velocity = getOptionalDouble(data[kPreProcesser]["velocity"]);

	// Initilaize optional minimum and maximum inclusion frequency bands
	std::pair<double, double> inclusionBand;

	double minFrequency = getOptionalDouble(data[kPreProcesser]["minFrequency"]);
	double maxFrequency = getOptionalDouble(data[kPreProcesser]["maxFrequency"]);

	inclusionBand.first = minFrequency != NULL ? minFrequency : DBL_MIN;
	inclusionBand.second = maxFrequency != NULL ? maxFrequency : DBL_MAX;
	
	configuration->preProcessor.inclusionBands = { inclusionBand };

	// Initialize optional exlcusion frequency bands
	std::pair<double, double> exclusionBands = std::make_pair(DBL_MAX, DBL_MIN);
	std::vector<double> exclusionFrequencies = getOptionalDoubleVector(data[kPreProcesser]["exclusionBands"]);

	if (exclusionFrequencies.size() % 2 != 0) {
		throw std::runtime_error("[ Json Config Parser ] Exclusion bands consist of a start and stop frequency.");
	}

	if (exclusionFrequencies.size() == 0) {
		configuration->preProcessor.exclusionBands = { exclusionBands };
	}
	else {
		std::pair<double, double> band;

		for (int i = 0; i < exclusionFrequencies.size(); i+=2) {
			band = std::make_pair(exclusionFrequencies[i], exclusionFrequencies[i + 1]);
			configuration->preProcessor.exclusionBands.emplace_back(band);
		}
	}
}


/**
 * Parses the scatter parameters in the 'processor/scatter' element.
 */
void JsonConfigParser::parseScatterParameters() {
	try {
		configuration->processor.maxScatterThreads = data[kProcesser][kScatter]["maxThreads"].template get<int>();
	}
	catch (const json::type_error& e) {
		printErrorMessage("maxThreads", kProcesser + "/" + kScatter);
		throw;
	}
}

/**
 * Sets parameters that are dependent on parameters from other elements.
 * 
 * This is an unfortunate consequence of poor structure design. Setting 
 * variables here prevents the need for the user to redefine the input
 * parameters multiple times.
 */
void JsonConfigParser::setEntangledParameters() {
	configuration->map.rfiScale = configuration->processor.rfiScaleBW;
	configuration->map.correlatedWeightMap = configuration->photometry.perform;
	configuration->map.processedWeightScale = configuration->processor.wScaleBW;

	configuration->processor.wCorrMap = configuration->photometry.perform;
}


/**
 * Returns the value for a JSON double that is allowed to be null.
 */
double JsonConfigParser::getOptionalDouble(json element) {
	if (element.is_null()) {
		return NULL;
	}
	try {
		return element.template get<double>();
	}
	catch (const json::type_error& e) {
		std::cerr << "[ Json Config Parser ] Invalid double value for 'preProcessor' element." << std::endl;
		throw;
	}
}

/**
 * Returns a vector of values for a JSON list that is allowed to be null.
 */
std::vector<double> JsonConfigParser::getOptionalDoubleVector(json element) {
	if (element.is_null()) {
		return {};
	}
	else {
		try {
			return element.template get<std::vector<double>>();
		}
		catch (const json::type_error& e) {
			std::cerr << "[ Json Config Parser ] Invalid list value for 'preProcessor' element." << std::endl;
			throw;
		}
	}
}

/**
 * Outputs an error message to the console given a json::type_error when reading values.
 */
void JsonConfigParser::printErrorMessage(std::string value, std::string element) {
	std::cerr << "[ Json Config Parser ] Invalid '" + value + "' value in '" + element + "' element." << std::endl;
}

JsonConfigParser::~JsonConfigParser() {

}