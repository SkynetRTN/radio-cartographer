#pragma once
#include <nlohmann/json.hpp>
#include "io\config\ConfigParserTemplate.h"
#include "io\config\ConfigParser.h"


class JsonConfigParser : ConfigParserTemplate {
public:
	JsonConfigParser(ConfigParser&);
	JsonConfigParser(ConfigParser&, bool);

	void parse();
	void parseMapParameters();
	void parseSurveyParameters();
	void parsePhotometryParameters();
	void parseProcessingParameters();
	void parsePreProcessingParameters();

	~JsonConfigParser();

private:
	bool logToConsole;
	nlohmann::json data;
	ConfigParser *configuration;
};