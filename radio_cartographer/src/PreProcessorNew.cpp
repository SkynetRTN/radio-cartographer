#include <iostream>
#include <future>

#include "PreProcessorNew.h"
#include "processor\BackgroundCUDA.h"
#include "utils\Tools.h"
#include "utils\RCR.h"


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

PreProcessorNew::PreProcessorNew(PreProcessingParameters params) : PreProcessorNew(params, 8, 8) {}

PreProcessorNew::PreProcessorNew(PreProcessingParameters params, int bgThreads, int scatterThreads)
{
	this->params = params;
	this->_MAX_BACKGROUND_THREADS = bgThreads;    // unc server = 50
	this->_MAX_SCATTER_THREADS = scatterThreads;  // unc server = 125
}

void PreProcessorNew::process(Observation& input) 
{
	std::vector<double> scatters;
	std::vector<std::vector<double>> background;

	setRawData(input);
	constructFrequencies(input);
	
	if (params.perform) 
	{
		constructWeights();
		constructBaselines();
		interpolateFluxDropouts();

		scatters = calculateScatterMulti(leftRawSpectra);
		leftCleanSpectra = performCleaningMulti(leftRawSpectra, scatters, params.subtractionScale);

		scatters = calculateScatterMulti(rightRawSpectra);
		rightCleanSpectra = performCleaningMulti(rightRawSpectra, scatters, params.subtractionScale);
	}
	else {
		leftCleanSpectra = leftRawSpectra;
		rightCleanSpectra = rightRawSpectra;
	}

	averageSpectra(input);
}

/**
 * Sets the to-be-cleaned data based on the data resolution and user-provided receiver.
 */
void PreProcessorNew::setRawData(Observation input) 
{
	if (input.observedFrequencies.size() != 2) {
		throw "Expected two observing frequencies (equal for low res data).";
	}

	// Default to lower frequency, override in next step if necessary
	observingFrequency = input.observedFrequencies[0];

	if (input.getDataResolution() == Observation::Resolution::LOW) {
		leftRawSpectra = input.lSpectra;
		rightRawSpectra = input.rSpectra;
	}
	else if (input.getDataResolution() == Observation::Resolution::HIGH) {
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
		throw std::invalid_argument("Unsupported data resolution provided.");
	}

	if (leftRawSpectra.empty() || rightRawSpectra.empty()) {
		throw std::runtime_error("Empty data set in either the left or right channel.");
	}

	if (leftRawSpectra.size() != rightRawSpectra.size()) {
		throw std::runtime_error("Mismatch in the number of spectra for the left and right channels.");
	}

	// Trim the data to the integrated points only
	originalDataSize = leftRawSpectra[0].size();
	for (int i = 0; i < leftRawSpectra.size(); ++i) {
		leftRawSpectra[i].erase(leftRawSpectra[i].begin() + input.getHigherIntegratedIndex(), leftRawSpectra[i].end());
		leftRawSpectra[i].erase(leftRawSpectra[i].begin(), leftRawSpectra[i].begin() + input.getLowerIntegratedIndex());

		rightRawSpectra[i].erase(rightRawSpectra[i].begin() + input.getHigherIntegratedIndex(), rightRawSpectra[i].end());
		rightRawSpectra[i].erase(rightRawSpectra[i].begin(), rightRawSpectra[i].begin() + input.getLowerIntegratedIndex());
	}
	dataSize = leftRawSpectra[0].size();
}

/**
 * Averages each spectra to a single continuum data point.
 */
void PreProcessorNew::averageSpectra(Observation& input)
{
	int numberOfPoints = leftCleanSpectra.size();

	for (int i = 0; i < numberOfPoints; ++i) 
	{
		if (leftCleanSpectra[i].size() != rightCleanSpectra[i].size())
			throw std::runtime_error("Mismatch in the number of data points for spectra " + std::to_string(i));

		int pointsInSum = 0;
		double leftSum = 0, rightSum = 0;

		for (int j = 0; j < dataSize; ++j) {
			if (!Tools::isInRange(frequencies[j], params.exclusionBands) && Tools::isInRange(frequencies[j], params.inclusionBands)) 
			{
				leftSum += leftCleanSpectra[i][j];
				rightSum += rightCleanSpectra[i][j];
				pointsInSum += 1;
			}
		}
		input.lContinuum.emplace_back(leftSum / pointsInSum);
		input.rContinuum.emplace_back(rightSum / pointsInSum);
	}
}

/**
 * Performs a multi-threaded background subtraction routine. The measured
 * background is saved as the cleaned data.
 */
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

	return 0.8197 * rcr.result.sigma;
}

/**
 * Constructs the weight vector for background subtraction.
 *
 * There are no weights associated with the input data. However, a weight vector
 * is needed for point rejection, so we construct a dummy vector here.
 */
void PreProcessorNew::constructWeights()
{
	weights = std::vector<double>(dataSize, 1.0);
}

/**
 * Converts index values to their coresponding frequency values. 
 *
 * Constructs the frequency vector for the original data size and the truncates
 * it to the trimmed data size. This is done to prevent messy logic that results
 * from data files only storing the center frequency value and the possibility 
 * of integrated ranges being before, intersecting or after the center.
 */
void PreProcessorNew::constructFrequencies(Observation input)
{
	double step = std::abs(input.frequencyStep);

	frequencies.resize(originalDataSize);

	// Construct the lower half
	for (int i = 0; i <= originalDataSize / 2; ++i)
		frequencies[(originalDataSize / 2) - i] = observingFrequency + (i * step);

	// Construct the upper half
	for (int i = originalDataSize / 2; i < originalDataSize; ++i)
		frequencies[i] = observingFrequency - ((i - (originalDataSize / 2)) * step);

	// Keep only the values in the integration range
	frequencies.erase(frequencies.begin() + input.getHigherIntegratedIndex(), frequencies.end());
	frequencies.erase(frequencies.begin(), frequencies.begin() + input.getLowerIntegratedIndex());

	frequencyDistances.resize(frequencies.size());
	for (int i = 0; i < frequencyDistances.size(); ++i) {
		frequencyDistances[i] = i * step;
	}

	// Check if the inclusion and exclusion bands excludes all data points
	int includedCount = 0;
	for (int i = 0; i < frequencies.size(); ++i)
	{
		if (!Tools::isInRange(frequencies[i], params.exclusionBands) &&
			Tools::isInRange(frequencies[i], params.inclusionBands))
		{
			includedCount += 1;
		}
	}

	if (includedCount == 0) {
		throw std::runtime_error("The inclusion or exclusion bands provided excluded all data.");
	}
}

/**
 * Constructs the baseline vector used for background subtraction. Accounts
 * for user input that modifies the baselines.
 */
void PreProcessorNew::constructBaselines()
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
		if (Tools::isInRange(frequencies[i], params.modifiedSubtractionBands)) {
			baselines[i] = params.modifiedSubtractionScale * pow(10, 6);
		}
		else if (params.velocity != NULL && Tools::isInRange(frequencies[i], { velocityRange })) {
			baselines[i] = 0.0;
		}
		else{
			baselines[i] = params.subtractionScale * pow(10, 6);
		}
	}
}

/**
 * Replaces invalid (negative) flux values using an interpolated value from the two 
 * nearest valid (non-negative) flux values.
 */
void PreProcessorNew::interpolateFluxDropouts() 
{
	double minFluxValue = 0.0;

	for (int i = 0; i < leftRawSpectra.size(); ++i) {
		bool interpolateLeft = false;
		bool interpolateRight = false;

		// Right and left data should always contain same number of points
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

std::vector<std::vector<double>> PreProcessorNew::getLeftCleanSpectra() {
	return leftCleanSpectra;
}

std::vector<std::vector<double>> PreProcessorNew::getRightCleanSpectra() {
	return rightCleanSpectra;
}

PreProcessorNew::~PreProcessorNew()
{
}