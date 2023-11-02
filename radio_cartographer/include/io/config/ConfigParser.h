#pragma once
#include <string>
#include "Structures.h"

class ConfigParser {
public:
	ConfigParser();
	ConfigParser(std::string);

	enum class FileType {
		JSON,
	};

	void parse();
	bool exists(std::string);

	std::string filename;
	MapParameters map;
	SurveyParameters survey;
	PhotometryParams photomtery;
	ProcessorParameters processor;
	PreProcessingParameters preProcessor;

	template <typename Enum>
	Enum enumFromString(const std::string&);

	int intFromBool(const bool&);

private:
	FileType fileType;

	void setFileType(std::string);

};