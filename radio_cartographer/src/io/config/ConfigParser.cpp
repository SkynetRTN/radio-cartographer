#include <iostream>
#include <fstream>

#include "io\config\ConfigParser.h"
#include "io\config\JsonConfigParser.h"


using json = nlohmann::json;


ConfigParser::ConfigParser() {}

ConfigParser::ConfigParser(std::string filename) {
	if (!exists(filename)) {
		throw "The provided file does not exsit: " + filename;
	}

	setFileType(filename);
	this->filename = filename;
}

/**
 * Determines if the file exists.
 *
 * @return true if the file exists, else false.
 */
bool ConfigParser::exists(std::string filename) {
	std::string mode = "r"; // opens non-text file formats

	FILE* file = std::fopen(filename.c_str(), mode.c_str());

	if (!file) {
		return false;
	}
	std::fclose(file);

	return true;
}

/**
 * Determines the type of the file. Supported formats: JSON.
 *
 * @return true if the file exists, else false.
 */
void ConfigParser::setFileType(std::string filename) {
	if (filename.find(".json") != std::string::npos) {
		fileType = FileType::JSON;
	}
	else {
		throw "Unsupported file type provided: " + filename + "\n";
	}
}

/**
 * Parses a configuration file. For required fields, see the documentation.
 *
 * For developers: If you would like to add your own configuration file
 * format, write a new derived parser class, e.g., XmlConfigParser, and
 * add an 'else if' statement here.
 */
void ConfigParser::parse() {
	if (fileType == FileType::JSON) {
		JsonConfigParser jsonParser = JsonConfigParser(*this, true);
		jsonParser.parse();
	}
	else {
		throw "Unsupported config file received.";
	}
}

int ConfigParser::intFromBool(const bool& b) {
	return b ? 1 : 0;
}

template<>
CalMethods ConfigParser::enumFromString(const std::string& s) {
	if (s == "pre" || s == "PRE") return CalMethods::PRE;
	if (s == "post" || s == "POST") return CalMethods::POST;
	if (s == "interpolated" || s == "INTERPOLATED") return CalMethods::INTERPOLATED;
	throw "Unexpected Enum value for CalMethods: " + s + ". Valid CalMethods: 'pre', 'post', or 'interpolated'";
}

template<>
Channel ConfigParser::enumFromString(const std::string& s) {
	if (s == "left" || s == "LEFT") return Channel::LEFT;
	if (s == "right" || s == "RIGHT") return Channel::RIGHT;
	if (s == "composite" || "COMPOSITE") return Channel::COMPOSITE;
	throw "Unexpected Enum value for Channel: " + s + ". Valid Channels: 'left', 'right', or 'composite'";
}

template<>
Coordinates ConfigParser::enumFromString(const std::string& s) {
	if (s == "equatorial" || s == "EQUATORIAL") return Coordinates::EQUATORIAL;
	if (s == "galactic" || s == "GALACTIC") return Coordinates::GALACTIC;
	throw "Unexpected Enum value for Coordinate: " + s + ". Valid Coordinates: 'equatorial' or 'galactic'";
}

template<>
CentroidMethods ConfigParser::enumFromString(const std::string& s) {
	if (s == "center" || s == "CENTER") return CentroidMethods::CENTER;
	if (s == "brightest" || s == "BRIGHTEST") return CentroidMethods::BRIGHTEST;
	if (s == "coordinates" || s == "COORDINATES") return CentroidMethods::COORDINATES;
	throw "Unexpected Enum value for CentroidMethods: " + s + ". Valid CentroidMethods: 'center', 'brightest', or 'coordinates'";
}

template<>
Receiver2 ConfigParser::enumFromString(const std::string& s) {
	if (s == "high" || s == "HIGH") return Receiver2::HIGH;
	if (s == "low" || s == "LOW") return Receiver2::LOW;
	throw "Unexpected Enum value for Receiver: " + s + ". Valid Receviers: 'high', or 'low'";
}