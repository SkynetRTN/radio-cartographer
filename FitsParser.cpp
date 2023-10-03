#include "FitsParser.h"
#include "HighResInput.h"

#include <CCfits/CCfits>
#include <iostream>


/**
 * Parses the Single Dish Flexible Image Transport System (SDFITS) file format.
 * See: https://fits.gsfc.nasa.gov/registry/sdfits.html for documentaion.
 */

FitsParser::FitsParser(std::string filename) {
	this->input.filename = filename;
	this->logToConsole = false;
}

FitsParser::FitsParser(std::string filename, bool logToConsole) {
	this->input.filename = filename;
	this->logToConsole = logToConsole;
}

Input FitsParser::parse() {

	parseMetaData();
	parseExtensionTable();

	if (input.getDataResolution() == Input::Resolution::HIGH) {
		sortHighResoultionData();
	}
	else if (input.getDataResolution() == Input::Resolution::LOW) {
		sortLowResolutionData();
	}
	return input;
}

/**
 * Opens and parses data from the FITS primary header.
 */
int FitsParser::parseMetaData() {
	if (logToConsole) std::cout << "[ INFO ] Parsing metadata" << std::endl;

	// Open the primary header
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, true));
	CCfits::PHDU& header = pInfile->pHDU();

	std::string coordinate;
	header.readKey(kTelescope, input.telescope);
	header.readKey(kObsMode, input.mapPattern);
	header.readKey(kCoordinate, coordinate);

	// Format coordinate system to something more readable
	if (coordinate.find("RA_DEC") != coordinate.npos) input.coordinate = "equatorial";
	else if (coordinate.find("LNG_LAT") != coordinate.npos) input.coordinate = "galactic";
	else throw "Unexpected coordinate system: " + coordinate + "\n";

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
	if (logToConsole) std::cout << "[ INFO ] Parsing history information" << std::endl;

	size_t pos = 0;
	std::string line, delimeter = "\n";

	while ((pos = history.find(delimeter)) != std::string::npos) {
		line = history.substr(0, pos);

		if (line.find(kIntFrequencyRange) != line.npos) parseIntegratedIndexRange(line);
		else if (line.find(kDataMode) != line.npos) parseDataResolution(line);

		history.erase(0, pos + delimeter.length());
	}
}

/**
 * Parses the binary table of the FITS file.
 */
int FitsParser::parseExtensionTable() {
	if (logToConsole) std::cout << "[ INFO ] Parsing binary table" << std::endl;

	// Open the binary table
	std::string tableName = "SINGLE DISH"; // Name of binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, tableName, false));
	CCfits::ExtHDU& table = pInfile->extension(tableName);

	// GBO alternates frequencies storing the higher one first
	table.column(kObsFrequency).read(input.observedFrequencies, 1, 2);

	// Read in the observed data
	std::vector<valarray<double>> vaSpectra;
	table.column(kSpectrumData).readArrays(vaSpectra, 1, table.rows());

	// Reverse the data since it is intuitvely backwards
	for (int i = 0; i < vaSpectra.size(); i++) {
		std::reverse(std::begin(vaSpectra[i]), std::end(vaSpectra[i]));
		spectra.emplace_back(std::begin(vaSpectra[i]), std::end(vaSpectra[i]));
	}

	// Reverse integrated ranges to match data
	int start = spectra[0].size() - input.getHigherIntegratedIndex();
	int stop = spectra[0].size() - input.getLowerIntegratedIndex();
	input.setIntegratedIndexRange(std::vector<int> {start, stop});

	// Read in all other relevant data
	table.column(kUtcSecs).read(times, 1, table.rows());
	table.column(kRaDegs).read(ras, 1, table.rows());
	table.column(kDecDegs).read(decs, 1, table.rows());
	table.column(kAzimuth).read(azimuths, 1, table.rows());
	table.column(kElevation).read(elevations, 1, table.rows());
	table.column(kDataDump).read(dataDumps, 1, table.rows());
	table.column(kCalibrationState).read(calibrations, 1, table.rows());
	table.column(kSweepIndex).read(sweepIndices, 1, table.rows());
	table.column(kSweepValid).read(valids, 1, table.rows());

	// Frequency step is used to convert from indices to frequencies later
	std::vector<double> frequencySteps;
	table.column(kFrequencyStep).read(frequencySteps, 1, 2);

	if (!frequencySteps.empty()) input.frequencyStep = frequencySteps[0];
	else throw "Input file is missing frequency step information.\n";

	return 0;
}

/**
 * Sort the FITS data for a low resolution observation.
 *
 * Low resolution files have data stored as alternating channels. Separate
 * these channels and store the information in the Input object.
 */
void FitsParser::sortLowResolutionData() {
	int numberOfRows = times.size() / 2;

	for (int i = 0; i < numberOfRows; ++i) {
		input.times.emplace_back(times[2 * i]);
		input.ras.emplace_back(ras[2 * i] / 15.0);
		input.decs.emplace_back(decs[2 * i]);
		input.azimuths.emplace_back(azimuths[2 * i]);
		input.elevations.emplace_back(elevations[2 * i]);
		input.calibrations.emplace_back(calibrations[2 * i]);
		input.dataDumps.emplace_back(dataDumps[2 * i]);
		input.valids.emplace_back(valids[2 * i]);
		input.sweepIndices.emplace_back(sweepIndices[2 * i]);

		// Separate into left and right channels
		input.lSpectra.emplace_back(spectra[2 * i]);
		input.rSpectra.emplace_back(spectra[2 * i + 1]);
	}
}

void FitsParser::sortHighResoultionData() {
	int numberOfRows = times.size() / 4;

	for (int i = 0; i < numberOfRows; ++i) {
		input.times.emplace_back(times[4 * i]);
		input.ras.emplace_back(ras[4 * i] / 15.0);
		input.decs.emplace_back(decs[4 * i]);
		input.azimuths.emplace_back(azimuths[4 * i]);
		input.elevations.emplace_back(elevations[4 * i]);
		input.calibrations.emplace_back(calibrations[4 * i]);
		input.dataDumps.emplace_back(dataDumps[4 * i]);
		input.valids.emplace_back(valids[4 * i]);
		input.sweepIndices.emplace_back(sweepIndices[4 * i]);

		// Separate into right channel
		input.rHighFrequencySpectra.emplace_back(spectra[4 * i]);
		input.lHighFrequencySpectra.emplace_back(spectra[4 * i + 2]);

		// Separate into left channel
		input.rLowFrequencySpectra.emplace_back(spectra[4 * i + 1]);
		input.lLowFrequencySpectra.emplace_back(spectra[4 * i + 3]);
	}
}

/**
 * Parses the resolution that the observation was observed in.
 */
void FitsParser::parseDataResolution(std::string line) {

	if (line.find("LOWRES") != line.npos) {
		input.setDataResolution(Input::Resolution::LOW);
	}
	else if (line.find("HIRES") != line.npos) {
		input.setDataResolution(Input::Resolution::HIGH);
	}
	else {
		throw "Unable to find the observation's resolution.\n";
	}
}

/**
 * Parses the start and stop integration index from the HISTORY string.
 * 
 * Assumes that the format of the HISTORY line is: 
 * HISTORY START,STOP channels  XXX, YYY   / range for continuum
 * where XXX is the start integer and YYY is the stop integer.
 */
void FitsParser::parseIntegratedIndexRange(std::string frequencyLine) {
	int k = 0, whiteSpace = 0;
	
	// Find the position after the intial text
	while (k < frequencyLine.size() && whiteSpace < 3)
	{
		if (frequencyLine.at(k) == ' ') whiteSpace++;
		k++;
	}

	// Exclude the initial text
	std::string ranges = frequencyLine.substr(k, 8);

	int start, stop;
	try {
		size_t pos = ranges.find(',');
		start = std::stoi(ranges.substr(0, pos));
		stop = std::stoi(ranges.substr(pos + 2, ranges.size()));
	}
	catch (const char* e) {
		std::cerr << e << "\n\n";
		throw "An error occured while reading the integration range from the primary header.\n";
	}
	if (start < 0 || stop < 0) {
		throw "Parsed a negative index range: " + std::to_string(start) + ", " + std::to_string(stop) + "\n";
	}
	input.setIntegratedIndexRange(std::vector<int> {start, stop});
}

FitsParser::~FitsParser() {

}