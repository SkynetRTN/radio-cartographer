#include <iostream>
#include <future>

#include "PreProcessorNew.h"
#include "BackgroundCUDA.h"
#include "Tools.h"
#include "RCR.h"


int _MAX_BACKGROUND_THREADS = 8; //50
int _MAX_SCATTER_THREADS = 8; //125

class PointToPointFunc : public NonParametric
{
public:
	PointToPointFunc(std::vector<double>, std::vector<double>);

	void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);
	void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::vector<double> workingChannel, angDist;
	std::vector<int> indices;
	~PointToPointFunc();

};

PreProcessorNew::PreProcessorNew() 
{

}

void PreProcessorNew::process(PreProcessingParameters params, Input& input) 
{
	std::vector<double> scatters;
	std::vector<std::vector<double>> background;

	setPreProcessingParams(params, input);
	
	if (params.subtractionScale != 0.0) {
		constructCleaningParams(params, input);

		scatters = calculateScatterMulti(leftRawSpectra);
		leftCleanSpectra = performCleaningMulti(leftRawSpectra, scatters, params.subtractionScale);

		scatters = calculateScatterMulti(rightRawSpectra);
		rightCleanSpectra = performCleaningMulti(rightRawSpectra, scatters, params.subtractionScale);
	}
	else {
		leftCleanSpectra = leftRawSpectra;
		rightCleanSpectra = rightRawSpectra;
	}

	averageSpectra(params, input);
}

/**
 * Sets the to-be-cleaned data based on the data resolution and user-provided receiver.
 */
void PreProcessorNew::setPreProcessingParams(PreProcessingParameters params, Input input) 
{
	if (input.observedFrequencies.size() != 2) {
		throw "Expected two observing frequencies (equal for low res data).";
	}

	// Default to lower frequency, override in next step if necessary
	observingFrequency = input.observedFrequencies[0];

	if (input.getDataResolution() == Input::Resolution::LOW) {
		leftRawSpectra = input.lSpectra;
		rightRawSpectra = input.lSpectra;
	}
	else if (input.getDataResolution() == Input::Resolution::HIGH) {
		if (params.receiver == Receiver2::LOW) {
			leftRawSpectra = input.lLowFrequencySpectra;
			rightRawSpectra = input.rLowFrequencySpectra;

			observingFrequency = input.observedFrequencies[1];
		}
		else if (params.receiver == Receiver2::HIGH) {
			leftRawSpectra = input.lHighFrequencySpectra;
			rightRawSpectra = input.rHighFrequencySpectra;
		}
	}
	else {
		throw "Unsupported data resolution provided.";
	}

	if (leftRawSpectra.empty() || rightRawSpectra.empty()) {
		throw "Empty data set in either the left or right channel.";
	}

	if (leftRawSpectra.size() != rightRawSpectra.size()) {
		throw "Mismatch in the number of spectra for the left and right channels.";
	}

	// Trim the data to the integrated points only
	for (int i = 0; i < leftRawSpectra.size(); ++i) {
		leftRawSpectra[i].erase(leftRawSpectra[i].begin(), leftRawSpectra[i].begin() + input.getLowerIntegratedIndex());
		leftRawSpectra[i].erase(leftRawSpectra[i].begin() + (input.getHigherIntegratedIndex() - input.getLowerIntegratedIndex() + 1), leftRawSpectra[i].end());

		rightRawSpectra[i].erase(rightRawSpectra[i].begin(), rightRawSpectra[i].begin() + input.getLowerIntegratedIndex());
		rightRawSpectra[i].erase(rightRawSpectra[i].begin() + (input.getHigherIntegratedIndex() - input.getLowerIntegratedIndex() + 1), rightRawSpectra[i].end());
	}

	dataSize = leftRawSpectra[0].size();
}

/**
 * Averages each spectra to a single continuum data point.
 */
void PreProcessorNew::averageSpectra(PreProcessingParameters params, Input& input)
{
	if (params.exclusionBand.empty())
		params.exclusionBand.emplace_back(DBL_MAX);
		params.exclusionBand.emplace_back(DBL_MIN);

	int numberOfPoints = leftCleanSpectra.size();

	for (int i = 0; i < numberOfPoints; ++i) {

		if (leftCleanSpectra[i].size() != rightCleanSpectra[i].size())
			throw "Mismatch in the number of data points for spectra " + std::to_string(i);

		int pointsInSum = 0;
		double leftSum = 0, rightSum = 0;

		for (int j = 0; j < dataSize; ++j) {
			leftSum += leftCleanSpectra[i][j];
			rightSum += rightCleanSpectra[i][j];

			pointsInSum += 1;
		}

		input.lContinuum.emplace_back(leftSum / pointsInSum);
		input.rContinuum.emplace_back(rightSum / pointsInSum);
	}
}


void PreProcessorNew::constructCleaningParams(PreProcessingParameters params, Input input)
{
	constructWeights();
	constructFrequencies(input);
	constructBaselines(params);
	interpolateFluxDropouts();
}

std::vector<std::vector<double>> PreProcessorNew::performCleaningMulti(std::vector<std::vector<double>> data, std::vector<double> scatters, double subtractionScale)
{
	int counter = 0, completedThreads = 0, liveThreads = 0; 
	std::vector<std::future<std::vector<double>>> futureVec;
	std::vector<double> results;

	BackgroundCUDA cuda;
	futureVec.resize(data.size());

	std::vector<std::vector<double>> background;
	background.resize(data.size());

	for (int i = 0; i < data.size(); i++)
	{
		cuda = BackgroundCUDA(baselines, frequencyDistances, weights, data[i], scatters[i]);

		futureVec[i] = std::async(std::launch::async, &BackgroundCUDA::calculateBG, cuda, subtractionScale * pow(10, 6));
		counter++;
		liveThreads++;

		if (liveThreads >= _MAX_BACKGROUND_THREADS)
		{
			for (int i = completedThreads; i < counter; i++)
			{
				results = futureVec[i].get();
				background[i] = results;
			}
			completedThreads += liveThreads;
			liveThreads = 0;
		}
	}

	for (int i = completedThreads; i < data.size(); i++)
	{
		results = futureVec[i].get();
		background[i] = results;
	}

	std::vector<std::vector<double>> cleanData(background.size(), std::vector<double>(dataSize));

	// Keep cleaned value if the baseline is not zero
	for (size_t i = 0; i < data.size(); ++i) {
		for (size_t j = 0; j < data[i].size(); ++j) {
			if (baselines[j] != 0.0) {
				cleanData[i][j] = background[i][j];
			}
			else {
				cleanData[i][j] = data[i][j];
			}
		}
	}
	
	return background;
}

/**
 * Multithreaded calling routine for scatter calculations. Max threads is defined by the 
 * the variable _MAX_SCATTER_THREADS at the top of this file.
 */
std::vector<double> PreProcessorNew::calculateScatterMulti(std::vector<std::vector<double>> data)
{
	int counter = 0, completedThreads = 0, liveThreads = 0;
	std::vector<std::future<double>> futureVec;
	std::vector<double> scatters;
	double result;

	futureVec.resize(data.size());
	scatters.resize(data.size());

	for (int i = 0; i < data.size(); i++)
	{
		futureVec[i] = std::async(std::launch::async, &PreProcessorNew::calculateScatter, this, data[i]);
		counter++;
		liveThreads++;

		if (liveThreads >= _MAX_SCATTER_THREADS)
		{
			for (int i = completedThreads; i < counter; i++)
			{
				result = futureVec[i].get();
				scatters[i] = result;
			}
			completedThreads += liveThreads;
			liveThreads = 0;
		}
	}

	for (int i = completedThreads; i < data.size(); i++)
	{
		result = futureVec[i].get();
		scatters[i] = result;
	}

	return scatters;
}

/**
 * Calculates the scatter for the given vector of data.
 */
double PreProcessorNew::calculateScatter(std::vector<double> data)
{
	RCR rcr = RCR(SS_MEDIAN_DL);
	PointToPointFunc ptpf(data, frequencyDistances);

	rcr.setNonParametricModel(ptpf);
	rcr.performRejection(weights, data);

	return 0.8197*rcr.result.sigma;
}

void PreProcessorNew::constructWeights()
{
	weights = std::vector<double>(dataSize, 1.0);
}

/**
 * Converts index values to their coresponding frequency values.
 */
void PreProcessorNew::constructFrequencies(Input input)
{
	double step = std::abs(input.frequencyStep);
	double zerothFrequency = observingFrequency - (step * (dataSize) / 2.0) + step;

	for (int i = input.getLowerIntegratedIndex(); i < input.getHigherIntegratedIndex() + 1; ++i)
	{
		frequencyDistances.emplace_back(step * (i - input.getLowerIntegratedIndex()));
		frequencies.emplace_back(zerothFrequency + (step * i));
	}
}

void PreProcessorNew::constructBaselines(PreProcessingParameters params)
{
	baselines.resize(dataSize);
	std::pair<double, double> velocityRange;

	// Convert user-provided velocity to frequency range
	if (params.velocity != 0.0) {
		const double c = 299792.458;
		const double f = 1420405751.7667;
		const double r = params.velocity / c;

		velocityRange = std::make_pair(f * (1.0 - r), f * (1.0 + r));
	}

	for (int i = 0; i < frequencies.size(); ++i) {
		if (isInRange(frequencies[i], params.modifiedSubtractionBands)) {
			baselines[i] = params.modifiedSubtractionScale;
		}
		else if (params.velocity != 0.0 && frequencies[i] >= velocityRange.first && frequencies[i] <= velocityRange.second) {
			baselines[i] = 0.0;
		}
		else
		{
			baselines[i] = params.subtractionScale;
		}
	}
}

/**
 * Replaces invalid values using an interpolated value from the two nearest valid values.
 */
void PreProcessorNew::interpolateFluxDropouts() 
{
	double minFluxValue = 0.0;

	for (int i = 0; i < leftRawSpectra.size(); ++i) {
		bool interpolateLeft = false;
		bool interpolateRight = false;

		// right and left data should always contain same number of points
		for (int j = 0; j < dataSize; ++j) {
			if (leftRawSpectra[i][j] < minFluxValue)
				interpolateLeft = true;
			if (rightRawSpectra[i][j] < minFluxValue)
				interpolateRight = true;
		}

		if (interpolateLeft)
			Tools::interpolate(leftRawSpectra[i], minFluxValue);
		if (interpolateRight)
			Tools::interpolate(rightRawSpectra[i], minFluxValue);
	}
}

bool PreProcessorNew::isInRange(double value, std::vector<std::pair<double, double>> ranges)
{
	for (const std::pair<double, double> &range : ranges) {
		if (value > range.first && value < range.second) {
			return true;
		}
	}
	return false;
}

std::vector<std::vector<double>> PreProcessorNew::getLeftCleanSpectra() {
	return leftCleanSpectra;
}

std::vector<std::vector<double>> PreProcessorNew::getRightCleanSpectra() {
	return rightCleanSpectra;
}