#pragma once
#include <string>
#include <vector>
#include <utility>

#include "io\Observation.h"

class TwentyMeterParser {
public:
	TwentyMeterParser(std::string);
	TwentyMeterParser(std::string, bool);

	Observation input;
	Observation parse();

	~TwentyMeterParser();

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
	std::string kSpectrumData = "DATA";
	std::string kMJulianDate = "MJD"; 

	// FITS Binary Table Keys
	std::string kUtcSecs = "UTSECS";
	std::string kRaDegs = "CRVAL2";
	std::string kDecDegs = "CRVAL3";
	std::string kAzimuth = "AZIMUTH";
	std::string kElevation = "ELEVATIO";
	std::string kDataDump = "NSAMPS";
	std::string kCalibrationState = "CALSTATE";
	std::string kSweepIndex = "SWPINDEX";
	std::string kSweepValid = "SWPVALID";
	std::string kFrequencyStep = "CDELT1";

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

	std::vector<std::vector<double>> spectra;

	std::string filename;

	int parseMetaData();
	int parseExtensionTable();

	void sortLowResolutionData();
	void sortHighResoultionData();
	void parseHistoryData(std::string);
	void parseDataResolution(std::string);
	void parseIntegratedIndexRange(std::string);
};