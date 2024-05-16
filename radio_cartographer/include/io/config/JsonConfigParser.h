#pragma once
#include <nlohmann/json.hpp>
#include "io\config\ConfigParserTemplate.h"
#include "io\config\ConfigParser.h"


class JsonConfigParser : ConfigParserTemplate {
public:
	JsonConfigParser(ConfigParser&);
	JsonConfigParser(ConfigParser&, bool);

	void parse();
	void parseFile();
	void parseMetaData();
	void parseMapParameters();
	void parseSurveyParameters();
	void parsePhotometryParameters();
	void parseProcessingParameters();
	void parsePreProcessingParameters();
	
	void parseBackgroundParameters();
	void parseTimeShiftingParameters();
	void parseScatterParameters();
	void parseRFIParameters();

	void setEntangledParameters();

	~JsonConfigParser();

private:
	std::string jsonFilename;

	std::string kProcesser = "processor";
	std::string kPreProcesser = "preProcessor";
	std::string kBackground = "background";
	std::string kTimeShift = "timeShift";
	std::string kScatter = "scatter";
	std::string kRfi = "rfi";

	bool logToConsole;
	nlohmann::json data;
	ConfigParser *configuration;

	std::vector<double> getOptionalDoubleVector(nlohmann::json);
	double getOptionalDouble(nlohmann::json);

	void printErrorMessage(std::string, std::string);
};