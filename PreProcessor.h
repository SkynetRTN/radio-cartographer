#pragma once
#include <string>
#include <sstream>
#include <vector>
#include <CCfits/CCfits>

enum Receiver { XX, YY };

struct SpectralParameters
{
	double subScale;
	double modSubScale;
	double velocity;

	std::vector<std::string> files;
	std::vector<double> modSubZones;
	Receiver receiver;
};

class Spectra
{
public:
	Spectra();

	void setScatter(double, int);
	void setFlux(std::vector<std::vector<double>>);
	void setSpectra20(std::vector<std::vector<double>>&);
	void determineFrequency(double, double, double);
	void determineBaselines(SpectralParameters);
	void exciseDropOuts(std::vector<std::vector<double>>&);
	void interpolateDropOuts(std::vector<std::vector<double>>&);
	void insertPoints(std::vector<std::vector<double>>&, double, int, int, std::string);
	void setWeights();

	double getScatter(int);
	std::vector<double> getWeight(int);
	std::vector<double> getFreqDist(int);
	std::vector<double> getBaseline(int);
	std::vector<double> getFlux(int);

	std::vector<std::vector<double>> getFreqDists();
	std::vector<std::vector<double>> getWeights();
	std::vector<std::vector<double>> baselines; // Should live in SpectralParameters
	
	std::vector<int> excisedIndices;
	void keepBackground(std::vector<std::vector<double>>&, double);

	~Spectra();

private:
	void velocityToFrequency(SpectralParameters);

	double freqEnd;
	double freqStart;
	double freqDistEnd;

	std::vector<double> scatters;
	std::vector<std::vector<double>> flux;
	std::vector<std::vector<double>> freqDists;
	std::vector<std::vector<double>> frequencies;
	std::vector<std::vector<double>> weights;
};

class PreProcessor
{
public:
	PreProcessor();

	// Main
	std::vector<std::vector<double>> sdfitsReader(SpectralParameters);

	// Input functions
	int accessExtTable20m(std::vector<valarray<double>>&, int);
	int accessExtTable40f(int);
	int accessExtTableGBT(int);
	int primaryHeader();
	void accessFITS();

	void getHistoryInfo(const std::string);	

	// Scatter
	void determineScatterMulti(Spectra&);
	double calculateScatter(int);

	// Getters
	int getMJDate();
	bool getScanningDirection();
	double getYear();
	double getFrequency();
	std::string getTelescope();
	std::string getMapPattern();
	std::string getCoordinateSystem();

	// Spectral Cleaning
	void performCleaningMulti(Spectra, double);
	void setBG(std::vector<double>, int);
	void buildSpectra(Spectra&, SpectralParameters);
	void removeBG();

	~PreProcessor();

private:

	// Validating Functions
	void validator();
	void validateYear();
	void validateSystem();
	void validatePattern();
	void validateFrequency();

	// Determining Functions
	void determineTelescope();

	// Formatting Functions
	void formatCoordinateSystem();
	void formatObservationYear(std::string);
	void formatSpectra20(std::vector<valarray<double>>);
	void formatHistoryInput(std::string);
	void formatGBTInput(std::string, std::vector<string>, std::vector<double>&);

	// Setter Functions
	void setMJD(std::vector<double>);
	void setResolution(std::string);
	void setFilename(SpectralParameters);
	void setScatterParams(Spectra);

	// Misc Functions
	void sortLowResData();
	void sortHiResData(Receiver);
	void appendColumnData();
	int averageSpectra();

	// Misc Variables
	double centerOffset;
	std::vector<int> intRanges;
	std::vector<std::vector<double>> fdists;
	std::vector<std::vector<double>> weights;

	// Observation Parameters
	std::string telescope;
	double frequency;
	double year;
	bool hiRes;
	int MJDate;

	// Mapping Parameters
	bool scanningDirection;
	std::string pattern;
	std::string system;

	// Miscellaneous Variables
	std::vector<std::string> filename;

	// Data Containers
	std::vector<std::vector<double>> outputData;
	std::vector<std::vector<double>> inputData;

	// Flux Data
	std::vector<std::vector<double>> spectra20;
	std::vector<double> continuum;

	// Spectral Cleaning Vectors
	std::vector<std::vector<double>> bgOutput;
};
