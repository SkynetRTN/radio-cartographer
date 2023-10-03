#pragma once
#include <vector>
#include "Input.h"


enum class Receiver2 {
	HIGH,
	LOW,
};

struct PreProcessingParameters
{
	double velocity = 0.0;

	double subtractionScale = 0.0;
	double modifiedSubtractionScale = 0.0;;

	std::vector<double> exclusionBand;
	std::vector<std::pair<double, double>> inclusionBands;
	std::vector<std::pair<double, double>> modifiedSubtractionBands;

	Receiver2 receiver;
};

class PreProcessorNew {
public:

	PreProcessorNew();

	double observingFrequency;
	std::vector<std::vector<double>> leftRawSpectra, rightRawSpectra;

	void process(PreProcessingParameters, Input&);

	void constructWeights();
	void constructFrequencies(Input);
	void constructBaselines(PreProcessingParameters);
	void interpolateFluxDropouts();

	bool isInRange(double, std::vector<std::pair<double, double>>);

	std::vector<std::vector<double>> getLeftCleanSpectra();
	std::vector<std::vector<double>> getRightCleanSpectra();

private:
	
	int dataSize;
	std::vector<double> frequencies;
	std::vector<double> frequencyDistances;

	std::vector<double> weights;
	std::vector<double> baselines;

	std::vector<std::vector<double>> leftCleanSpectra, rightCleanSpectra;

	void constructCleaningParams(PreProcessingParameters, Input);
	void setPreProcessingParams(PreProcessingParameters, Input);
	void averageSpectra(PreProcessingParameters, Input&);

	double calculateScatter(std::vector<double>);
	std::vector<double> calculateScatterMulti(std::vector<std::vector<double>>);

	std::vector<std::vector<double>> performCleaningMulti(std::vector<std::vector<double>>, std::vector<double>, double);
};