#pragma once

#include <vector>

class SpectralBaseline
{
public:
	SpectralBaseline(bool forward, double scatter, double xAnchor, double yAnchor, std::vector<bool>& checks, std::vector<double>& flux, std::vector<double>& dataDumps, std::vector<double>& angDist);
	SpectralBaseline();
	double rejectPoints();
	int returnPoints();
	std::vector<double> getResults();
	void setLocalModels();

	~SpectralBaseline();

private:
	bool forward;
	double scatter;
	double xAnchor;
	double yAnchor;
	std::vector<double> BLFlux;
	std::vector<double> BLAngDist;
	std::vector<double> BLDataDumps;
	std::vector<bool> BLChecks;
	std::vector<double> BLResults;

	bool findDuplicate(double, std::vector<double>);
	std::vector<double> autoFixedWRegression();
	std::vector<double> autoWRegression();
	std::vector<double> applyModel(double, std::vector<double> &);

	int sufficentPointCheck();
};

class SpectralBackground
{
public:
	SpectralBackground();
	SpectralBackground(const std::vector<double>& inputAngDist, const std::vector<double>& inputFlux, const std::vector<double>& inputDataDumps, double inputScatter);
	std::vector<double> calculateBG(double baseline);
	std::vector<double> calculateBGMulti(double baseline);

	std::vector<double> baselineVec;

	~SpectralBackground();

private:
	double scatter;
	int size;
	std::vector<SpectralBaseline> baselineArray;
	std::vector<bool> checks;
	std::vector<double> angDist;
	std::vector<double> dataDumps;
	std::vector<double> flux;
	std::vector<std::vector<double>> bgData, bgWeights;

	std::vector<double> setBackground(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
	std::vector<std::vector<int>> findStartEndIndices(std::vector<double> &, double);
	void buildBaselines(bool, std::vector<int>&);
	void loadLocalModels(std::vector<double> &, int);
};
