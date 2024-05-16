#include <vector>
#include "Scan.h"
#include "PreProcessor.h"

class Baseline
{
public:
	Baseline(bool, double, double, double, std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	Baseline();
	double rejectPoints();
	int returnPoints();
	std::vector<double> getResults();
	void setLocalModels();

	~Baseline();

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

class BackgroundCUDA
{
public:
	BackgroundCUDA();
	BackgroundCUDA(Scan&, bool);
	BackgroundCUDA(Spectra&, int);
	BackgroundCUDA(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, double);
	std::vector<double> calculateBG(double);
	std::vector<double> calculateBGMulti(double);

	std::vector<double> baselineVec;

	~BackgroundCUDA();

private:
	double scatter;
	int size;
	std::vector<Baseline> baselineArray;
	std::vector<bool> checks;
	std::vector<double> angDist;
	std::vector<double> dataDumps;
	std::vector<double> dec;
	std::vector<double> flux;
	std::vector<double> ra;
	std::vector<std::vector<double>> bgData, bgWeights;

	std::vector<double> setBackground(std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
	std::vector<std::vector<int>> findStartEndIndices(std::vector<double> &, double);
	void buildBaselines(bool, std::vector<int>&);
	void loadLocalModels(std::vector<double> &, int);
};
