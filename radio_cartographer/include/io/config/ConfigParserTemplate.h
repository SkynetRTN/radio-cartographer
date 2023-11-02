#pragma once


class ConfigParserTemplate {
public:
	virtual void parse() = 0;
	virtual void parseMapParameters() = 0;
	virtual void parseSurveyParameters() = 0;
	virtual void parsePhotometryParameters() = 0;
	virtual void parseProcessingParameters() = 0;
	virtual void parsePreProcessingParameters() = 0;

	virtual ~ConfigParserTemplate() {};
};