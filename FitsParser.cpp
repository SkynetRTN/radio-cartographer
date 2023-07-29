#include "FitsParser.h"

#include <CCfits/CCfits>
#include <iostream>


/**
 * Parses the Single Dish Flexible Image Transport System (SDFITS) file format.
 * See: https://fits.gsfc.nasa.gov/registry/sdfits.html for documentaion.
 */

FitsParser::FitsParser(std::string filename) {
	this->filename = filename;
}

std::vector<std::vector<double>> FitsParser::parse() {
	parseMetaData();
	parseExtensionTable();
}

/**
 * Opens and parses data from the FITS primary header.
 */
int FitsParser::parseMetaData() {

	// Open the primary header
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename, CCfits::Read, true));
	CCfits::PHDU& header = pInfile->pHDU();

	header.readKey(kTelescope, telescope);
	header.readKey(kCoordinate, coordinate);
	header.readKey(kObsDate, obsDate);
	header.readKey(kObsMode, obsMode);

	// Format coordinate system to something more readable
	if (coordinate.find("RA_DEC") != coordinate.npos) coordinate = "equatorial";
	else if (coordinate.find("LNG_LAT") != coordinate.npos) coordinate = "galactic";
	else {
		std::cerr << "Unexpected coordinate system: " + coordinate << std::endl;
		std::exit(EXIT_FAILURE);
	}

	parseHistoryData(header.getHistory());

	return 0;
}

/**
 * Parses the HISTORY section of the FITS header.
 * 
 * The HISTORY section is a single string stored in the primary header.
 * Each line begins with the word HISTORY and is delimited by a line break.
 * 
 * Break the string into substring containing each line and search for 
 * the relevant keywords in each line.
 */
void FitsParser::parseHistoryData(std::string history) {

	size_t pos = 0;
	std::string line, delimeter = "\n";

	while ((pos = history.find(delimeter)) != std::string::npos) {
		line = history.substr(0, pos);

		if (line.find(kIntFrequencyRange) != line.npos) setIntegratedIndexRange(line);
		else if (line.find(kDataMode) != line.npos) parseDataResolution(line);

		history.erase(0, pos + delimeter.length());
	}
}

/**
 * Parses the binary table of the FITS file.
 */
int FitsParser::parseExtensionTable() {

	// Open the binary table
	std::string tableName = "SINGLE DISH"; // Name of binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename, CCfits::Read, tableName, false));
	CCfits::ExtHDU& table = pInfile->extension(tableName);

	// CCFits requires data to be read into a vector, but we want to store as a pair.
	std::vector<double> obsFrequencies;
	table.column(kObsFrequency).read(obsFrequencies, 1, 2);
	observedFrequencyRange = std::make_pair(obsFrequencies[0], obsFrequencies[1]);

	//
	// HELLO, might be able to loop through and emplace into a vector of vectors instead of a valarray.
	// https://heasarc.gsfc.nasa.gov/fitsio/CCfits/html/classCCfits_1_1Column.html#a8409427a876646e0f1846f55146872bc
	// might be able to get num of rows from dimen method
	//

	// Read in the observed data
	std::vector<valarray<double>> vasTemp;
	table.column(kSpectrumData).readArrays(vasTemp, 1, 999999);

	// Read all other relevant table data
	for (const std::string column : columns)
	{
		table.column(column).read(inputTemp[i], 1, 999999);
	}
}

/**
 * Parses the resolution that the observation was observed in.
 */
void FitsParser::parseDataResolution(std::string line) {

	if (line.find("LOWRES") != line.npos) {
		setDataResolution(Resolution::LOW);
	}
	else if (line.find("HIRES") != line.npos) {
		setDataResolution(Resolution::HIGH);
	}
	else {
		std::cerr << "Unable to find the observation's resolution.\n";
		std::exit(EXIT_FAILURE);
	}
}

/**
 * Parses the start and stop frequency from the HISTORY string.
 * 
 * The HISTORY section of the primary header for an SDFITS is
 * poorly formatted which makes parsing un-necessarily annoying.
 * 
 * Assumes that the format of the HISTORY line is: 
 * HISTORY START,STOP channels  XXX, YYY   / range for continuum
 * where XXX is the start integer and YYY is the stop integer.
 */
void FitsParser::setIntegratedIndexRange(std::string frequencyLine) {
	int k = 0, whiteSpace = 0;
	
	// Find the position after the intial text
	while (k < frequencyLine.size() && whiteSpace < 3)
	{
		if (frequencyLine.at(k) == ' ') whiteSpace++;
		k++;
	}

	int startIndex = k + 1; // Exclude the initial text
	std::string ranges = frequencyLine.substr(startIndex, 8);

	int32_t start, stop;
	try {
		size_t pos = ranges.find(',');
		start = std::stoi(ranges.substr(0, pos));
		stop = std::stoi(ranges.substr(pos + 2, ranges.size()));
	}
	catch (const char* e) {
		std::cerr << e << "\n\n";
		std::cerr << "An error occured while reading the frequency ";
		std::cerr << "range from the primary header.\n";
		std::exit(EXIT_FAILURE);
	}
	integratedIndexRange = std::make_pair(start, stop);
}

void FitsParser::setDataResolution(Resolution resolution) {
	this->resolution = resolution;
}

FitsParser::~FitsParser() {

}