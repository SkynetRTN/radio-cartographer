#pragma once
#include <vector>
#include "io\Observation.h"
#include "utils\Structures.h"


class PreProcessorNew {
public:

	PreProcessorNew(PreProcessingParameters);
	PreProcessorNew(PreProcessingParameters, int, int);

	double observingFrequency;
	std::vector<std::vector<double>> leftRawSpectra, rightRawSpectra;

	void process(Observation&);

	void constructWeights();
	void constructBaselines();
	void interpolateFluxDropouts();
	void constructFrequencies(Observation);

	std::vector<std::vector<double>> getLeftCleanSpectra();
	std::vector<std::vector<double>> getRightCleanSpectra();

	~PreProcessorNew();

private:
	PreProcessingParameters params;

	int _MAX_BACKGROUND_THREADS;
	int _MAX_SCATTER_THREADS;

	int dataSize;
	int originalDataSize;
	
	std::vector<double> weights;
	std::vector<double> baselines;
	std::vector<double> frequencies;
	std::vector<double> frequencyDistances;

	std::vector<std::vector<double>> leftCleanSpectra, rightCleanSpectra;

	void setRawData(Observation);
	void averageSpectra(Observation&);

	double calculateScatter(std::vector<double>);
	std::vector<double> calculateScatterMulti(std::vector<std::vector<double>>);
	std::vector<std::vector<double>> performCleaningMulti(std::vector<std::vector<double>>, std::vector<double>, double);
};