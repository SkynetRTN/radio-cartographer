#pragma once
#include <string>
#include <vector>
#include <utility>

class Input {
public:
	Input();
	Input(std::string);
	
	enum class Resolution {
		HIGH,
		LOW,
	};

	enum class FileType {
		MD2,
		TEXT,
		FITS
	};

	std::string filename;

	// Observation information
	double frequencyStep;
	std::string telescope;
	std::string coordinate;
	std::string mapPattern;
	std::string scanType;

	std::vector<double> observedFrequencies;

	// Observation data
	std::vector<double> times;
	std::vector<double> ras;
	std::vector<double> decs;
	std::vector<double> azimuths;
	std::vector<double> elevations;
	std::vector<double> calibrations;
	std::vector<double> dataDumps;
	std::vector<double> valids;
	std::vector<double> sweepIndices;

	std::vector<double> lContinuum;
	std::vector<double> rContinuum;

	// Low resolution observation Data
	std::vector<std::vector<double>> lSpectra;
	std::vector<std::vector<double>> rSpectra;

	// High resolution observation data
	std::vector<std::vector<double>> lLowFrequencySpectra;
	std::vector<std::vector<double>> rLowFrequencySpectra;
	std::vector<std::vector<double>> lHighFrequencySpectra;
	std::vector<std::vector<double>> rHighFrequencySpectra;

	void setDataResolution(Resolution);
	void setIntegratedIndexRange(std::vector<int>);

	Resolution getDataResolution();
	int getLowerIntegratedIndex();
	int getHigherIntegratedIndex();

	void convertFromColumnData(std::vector<std::vector<double>>);


private:
	Resolution resolution;
	std::vector<int> integratedIndexRange;
};
