#include "GBTParser.h"

#include <CCfits/CCfits>
#include <iostream>
#include <cstring>


/**
 * Parses the Single Dish Flexible Image Transport System (SDFITS) file format.
 * See: https://fits.gsfc.nasa.gov/registry/sdfits.html for documentaion.
 */

GBTParser::GBTParser(std::string filename) {
	this->input.filename = filename;
	this->logToConsole = false;
}

GBTParser::GBTParser(std::string filename, bool logToConsole) {
	this->input.filename = filename;
	this->logToConsole = logToConsole;
}

Input GBTParser::parse() {

	parseMetaData();
	parseExtensionTable();
	sortLowResolutionData();

	return input;
}

/**
 * Opens and parses data from the FITS primary header.
 */
int GBTParser::parseMetaData() {
	if (logToConsole) std::cout << "[ INFO ] Parsing metadata" << std::endl;

	// Open the primary header
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, true));
	CCfits::PHDU& header = pInfile->pHDU();

	std::string coordinate;
	header.readKey(kTelescope, input.telescope);

	input.setIntegratedIndexRange(std::vector<int> {0, 0});
	input.setDataResolution(Input::Resolution::LOW);

	return 0;
}

/**
 * Parses the binary table of the FITS file.
 */
int GBTParser::parseExtensionTable() {
	if (logToConsole) std::cout << "[ INFO ] Parsing binary table" << std::endl;

	// Open the binary table
	std::string tableName = "SINGLE DISH"; // Name of binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, tableName, false));
	CCfits::ExtHDU& table = pInfile->extension(tableName);

	std::vector<std::string> obsMode;
	table.column("OBSMODE").read(obsMode, 1, 1);
	parseObservationMode(obsMode[0]);

	std::vector<std::string> dates;
	table.column("DATE-OBS").read(dates, 1, table.rows());
	convertHMStoUTC(dates);

	// GBO alternates frequencies storing the higher one first
	table.column(kObsFrequency).read(input.observedFrequencies, 1, 2);

	// Read in all other relevant data
	table.column(kContinuumData).read(continuum, 1, table.rows());
	table.column(kRaDegs).read(ras, 1, table.rows());
	table.column(kDecDegs).read(decs, 1, table.rows());
	table.column(kAzimuth).read(azimuths, 1, table.rows());
	table.column(kElevation).read(elevations, 1, table.rows());
	table.column(kSweepIndex).read(sweepIndices, 1, table.rows());

	// Dummies for non-GBT specific values
	this->dataDumps = std::vector<double>(table.rows(), 1.0);
	this->calibrations = std::vector<double>(table.rows(), 0.0);
	this->valids = std::vector<double>(table.rows(), 1.0);

	input.coordinate = "equatorial"; // @TODO: Where to get this?

	return 0;
}


void GBTParser::convertHMStoUTC(std::vector<std::string> dates) {

	std::string hours, minutes, seconds;
	times.resize(dates.size());

	for (int i = 0; i < dates.size(); i++)
	{
		// Convert hh:mm::ss to utc
		hours = dates[i][11];
		hours += dates[i][12];

		minutes = dates[i][14];
		minutes += dates[i][15];

		for (int j = 17; j < 22; j++)
		{
			seconds += dates[i][j];
		}

		times[i] = ::atof(hours.c_str()) * 3600.00 + ::atof(minutes.c_str()) * 60.0 + ::atof(seconds.c_str());

		hours.clear(), minutes.clear(), seconds.clear();
	}
}

void GBTParser::parseObservationMode(std::string mode) {

	size_t pos = mode.find(':');
	if (pos != std::string::npos) {
		input.mapPattern = mode.substr(0, pos);
	}

	// convert to lower for a later comparison
	for (char& c : input.mapPattern) {
		c = tolower(c);
	}
}


/**
 * Sort the FITS data for a low resolution observation.
 *
 * Low resolution files have data stored as alternating channels. Separate
 * these channels and store the information in the Input object.
 */
void GBTParser::sortLowResolutionData() {
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
		input.sweepIndices.emplace_back(sweepIndices[2 * i] - 1);

		// Separate into left and right channels
		input.lContinuum.emplace_back(continuum[2 * i]);
		input.rContinuum.emplace_back(continuum[2 * i + 1]);
	}
}

GBTParser::~GBTParser() {

}