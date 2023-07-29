#pragma once
#include <string>
#include <vector>
#include <utility>

class FitsParser {
public:
	FitsParser(std::string);
	std::vector<std::vector<double>> parse();

private:
	enum class Resolution {
		HIGH,
		LOW,
	};

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
	std::vector<std::string> columns = { "UTSECS", "CRVAL2", "CRVAL3", "AZIMUTH", "ELEVATIO", 
										 "NSAMPS", "CALSTATE", "SWPINDEX", "SWPVALID", "CDELT1" };

	// Inputs
	std::string filename;

	// Extracted Data
	std::string telescope;
	std::string obsDate;
	std::string obsMode;
	std::string coordinate;

	Resolution resolution;
	std::pair<std::int32_t, std::int32_t> integratedIndexRange;
	std::pair<std::int32_t, std::int32_t> observedFrequencyRange;

	// Outputs
	std::vector<double> times;
	std::vector<double> ras;
	std::vector<double> decs;
	std::vector<double> azimuths;
	std::vector<double> elevations;
	std::vector<double> fluxesL;
	std::vector<double> fluxesR;
	std::vector<double> fluxesComp;
	std::vector<double> calibrations;
	std::vector<double> dataDumps;

	// Methods
	int parseMetaData();
	void parseHistoryData(std::string);
	int parseExtensionTable();

	void parseDataResolution(std::string);

	void setIntegratedIndexRange(std::string);
	void setDataResolution(Resolution);

	std::pair<double, double> getIntegratedFrequencyRange();
	std::pair<double, double> getObserveddFrequencyRange();
};