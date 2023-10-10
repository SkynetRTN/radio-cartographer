#pragma once
#include <string>
#include <vector>
#include <utility>

#include "Input.h"

class GBTParser {
public:
	GBTParser(std::string);
	GBTParser(std::string, bool);

	Input input;
	Input parse();

	~GBTParser();

private:
	bool logToConsole;

	// FITS Primary Header Keys
	std::string kObsDate = "DATE";
	std::string kObsMode = "OBSMODE";
	std::string kTelescope = "TELESCOP";
	std::string kCoordinate = "COORDREF";
	std::string kDataMode = "DATAMODE";
	std::string kIntFrequencyRange = "START,STOP";
	std::string kObsFrequency = "OBSFREQ";
	std::string kContinuumData = "DATA";
	std::string kMJulianDate = "MJD";

	// FITS Binary Table Keys
	std::string kUtcSecs = "UTSECS";
	std::string kRaDegs = "CRVAL2";
	std::string kDecDegs = "CRVAL3";
	std::string kAzimuth = "AZIMUTH";
	std::string kElevation = "ELEVATIO";
	std::string kSweepIndex = "PROCSEQN";

	// FITS Binary Table Data
	std::vector<double> times;
	std::vector<double> ras;
	std::vector<double> decs;
	std::vector<double> azimuths;
	std::vector<double> elevations;
	std::vector<double> dataDumps;
	std::vector<double> calibrations;
	std::vector<double> sweepIndices;
	std::vector<double> valids;

	std::vector<double> continuum;

	std::string filename;

	int parseMetaData();
	int parseExtensionTable();

	void convertHMStoUTC(std::vector<std::string>);
	void parseObservationMode(std::string);
	void sortLowResolutionData();
};