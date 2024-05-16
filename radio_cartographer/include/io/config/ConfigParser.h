#pragma once
#include <string>
#include "utils\Structures.h"

class ConfigParser {
public:
	ConfigParser();
	ConfigParser(std::string);

	enum class FileType {
		JSON,
	};

	std::string runMethod;
	std::string configFilename;
	std::vector<std::string> obsFilenames;

	void parse();
	bool exists(std::string);

	MapParameters map;
	SurveyParameters survey;
	PhotometryParameters photometry;
	ProcessorParameters processor;
	PreProcessingParameters preProcessor;

	template <typename Enum>
	Enum enumFromString(const std::string&);

	int intFromBool(const bool&);

private:
	FileType fileType;

	void setFileType(std::string);
};