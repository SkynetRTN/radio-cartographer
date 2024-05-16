#include "io\observation\HundredMeterParser.h"

#include <CCfits/CCfits>
#include <iostream>
#include <cstring>


/**
 * Parses the Single Dish Flexible Image Transport System (SDFITS) file format.
 * See: https://fits.gsfc.nasa.gov/registry/sdfits.html for documentaion.
 */

HundredMeterParser::HundredMeterParser(std::string filename) : HundredMeterParser(filename, false) {}

HundredMeterParser::HundredMeterParser(std::string filename, bool logToConsole) {
	this->input.filename = filename;
	this->logToConsole = logToConsole;
}

Observation HundredMeterParser::parse()
{
	parseMetaData();
	parseExtensionTable();
	sortLowResolutionData();

	return input;
}

/**
 * Opens and parses data from the FITS primary header.
 */
int HundredMeterParser::parseMetaData() {
	if (logToConsole) std::cout << "[ INFO ] Parsing metadata" << std::endl;

	// Open the primary header
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, true));
	CCfits::PHDU& header = pInfile->pHDU();

	std::string coordinate;
	header.readKey(kTelescope, input.telescope);

	input.setIntegratedIndexRange(std::vector<int> {0, 0});
	input.setDataResolution(Observation::Resolution::LOW);

	return 0;
}

/**
 * Parses the binary table of the FITS file.
 */
int HundredMeterParser::parseExtensionTable() {
	if (logToConsole) std::cout << "[ INFO ] Parsing binary table" << std::endl;

	// Open the binary table
	std::string tableName = "SINGLE DISH"; // Name of binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input.filename, CCfits::Read, tableName, false));
	CCfits::ExtHDU& table = pInfile->extension(tableName);

	// Read and parse observing mode(s)
	std::vector<std::string> obsMode;
	table.column(kObsMode).read(obsMode, 1, 1);
	input.mapPattern = parseObservationMode(obsMode[0]);
	table.column(kObsMode).read(obsModes, 1, table.rows());

	for (int i = 0; i < obsModes.size(); i += 2)
	{
		input.observingModes.emplace_back(parseObservationMode(obsModes[i]));
	}

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
	table.column(kCalibration).read(calibrationFlags, 1, table.rows());

	// Dummies for non-GBT specific values
	this->dataDumps = std::vector<double>(table.rows(), 1.0);
	this->calibrations = std::vector<double>(table.rows(), 0.0);
	this->valids = std::vector<double>(table.rows(), 1.0);

	input.coordinate = "equatorial"; // @TODO: Where to get this?

	return 0;
}


void HundredMeterParser::convertHMStoUTC(std::vector<std::string> dates) {

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

std::string HundredMeterParser::parseObservationMode(std::string mode) 
{
	std::string pattern;

	size_t pos = mode.find(':');
	if (pos != std::string::npos) {
		pattern = mode.substr(0, pos);
	}

	// convert to lower for a later comparison
	for (char& c : pattern) {
		c = tolower(c);
	}

	return pattern;
}


/**
 * Sort the FITS data for a low resolution observation.
 *
 * GBT files have a single polarization with alternating calibrations. 
 * This is contrasted with 20m files that have alternating polarizations.
 */
void HundredMeterParser::sortLowResolutionData()
{
	const int numberOfRows = times.size();

	input.mapStartLocations.emplace_back(0);

	int nextSweepIndex, currSweepIndex;
	for (int i = 0; i < numberOfRows; ++i)
	{
		currSweepIndex = sweepIndices[i];
		if (i < numberOfRows - 1) {
			nextSweepIndex = sweepIndices[i + 1];
		}
		else {
			nextSweepIndex = sweepIndices[i];
		}
		
		// Check if the sweep index reset back to 1
		if (currSweepIndex > nextSweepIndex) {
			input.hasMultipleMaps = true;
			input.mapStartLocations.emplace_back((i + 1) / 2);
		}

		if (calibrationFlags[i] == "F")
		{
			input.times.emplace_back(times[i]);
			input.ras.emplace_back(ras[i] / 15.0);
			input.decs.emplace_back(decs[i]);
			input.azimuths.emplace_back(azimuths[i]);
			input.elevations.emplace_back(elevations[i]);
			input.calibrations.emplace_back(calibrations[i]);
			input.dataDumps.emplace_back(dataDumps[i]);
			input.valids.emplace_back(valids[i]);
			input.sweepIndices.emplace_back(sweepIndices[i] - 1);

			input.lContinuum.emplace_back(continuum[i]);
		}
		else if (calibrationFlags[i] == "T")
		{   
			// Store calibration ON flux for calibration
			input.rContinuum.emplace_back(continuum[i]);
		}
	}
}

HundredMeterParser::~HundredMeterParser() {

}