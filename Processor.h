#pragma once
#include <vector>
#include <string>
#include <future>
#include "Survey.h"
#include "Composite.h"

extern bool LSSProcessing;

class Processor
{
public:
	Processor();
	Processor(Survey &);
	Processor(ProcessorParameters &);

	void characterizeSurvey(Survey &);

	void performTimeShifting(Survey &);
	void forceThetaGap(Composite &, double);
	void processLargeScaleStructure(Survey &);

	void set2DScatter(Survey &);
	void repartitionScans(Composite &composite, std::vector<Survey> &surveys);

	void determineDataBoundaries(Survey &);
	void characterizeData(Survey &);

	void performBGSubtractionMulti(Survey &);
	void performRFIRejectionMulti(Composite &, std::vector<Survey> &surveys);
	void calculateProcThetaGapMulti(Composite &);

	void levelLSSData(std::vector<Survey> &);

	~Processor();

private:
	Survey surveyHold;
	std::vector<Scan> scans;
	MapTypes mapType;
	Coordinates pCoordinate;
	PartitionSet partSetProcSSS;
	PartitionSet partSetProcLSS;
	int telescope;
	bool scansInRa;

	double medianDiffAlongSweeps;
	double psfFWHM;

	double bgScaleBW;
	double rfiScaleBW;
	double wScaleBW;
	std::string timeShift;
	bool wCorrMap;
	bool lssProc;

	double timeShiftValue;
	double globalMaxDec;
	double globalMinDec;
	double globalMaxRa;
	double globalMinRa;
	double globalCenterDec;
	double globalCenterRa;
	double globalCenterLatProc;
	double globalCenterLongProc;
	double medianLatiMapped;
	double medianLongMapped;
	bool trackingForEdges;

	std::string mCoordinate;
	std::vector<double> centroidLocations;
	std::vector<double> standardGaps;
	std::vector<std::vector<double>> scatter2dLSS;
	std::vector<std::vector<double>> scatter2dSSS;
	std::vector<std::vector<std::vector<int>>> classificationsSSS;
	std::vector<std::vector<std::vector<int>>> classificationsLSS;

	// data pre-processing
	void switchChannels(Channel, Survey &);
	void calculateScatter();

	// time-shifting
	std::vector<bool> shiftRejection(std::vector<double>);
	std::vector<std::vector<double>> getScanToScanShifts(double);
	double calculateShift(std::vector<std::vector<double>> &);

	// misc
	int findPossSSS(double, double, double, std::vector<int> &);
	int findPossLSS(double, double, double, std::vector<int> &);
	void classifySSS(Survey &);
	void classifySSS(Composite &);
	void classifyLSS(Survey &);
	void classifyLSS(Composite &);

	PartitionSet determineProcSurveyDimensions(Survey &, bool);
	PartitionSet determineProcCompositeDimensions(Composite &, bool);

	// large scale structure
	std::vector<double> setDropDeltas();
	void setDropFlags(std::vector<double> &);
	std::vector<double> determinePercentDrop(std::vector<int> &, std::vector<int> &);
	void correctDrop(std::vector<double>);

	void lssDropCorrection(Survey &);
	void lssPolynomialRFISubtraction(Survey &, double, double);
	void lssElevationSubtraction(Survey &);
	void lssThetaGapCalc(Survey &);
	void lssSet2DScatter(Survey &);

	void calculateSurveyElevationScatter();
	void calculateElevationBG(int, double);

	// double lssCircleThetaGap(int, int, std::vector<double> &, std::vector<bool> &);

	void rejectModelPoints(int, int, double, double, double, std::vector<double> &, std::vector<bool> &);
	double lssRecursivePolyFit(Survey &, int, int, double, double);
	bool planeApplication(int, int, int, std::vector<std::vector<int>>, int, MapTypes);
	int correctIndex(int, int, int, int);
	std::vector<std::vector<int>> calculateScanNumbers(int &, std::vector<double> &, std::vector<bool> &);
	// int lssFindPoss(double, double, double, double, std::vector<int> &);
	bool verifyOutlierState(SurfaceType, int, int, double, std::vector<double> &, std::vector<double> &);
	double calculateSecondDerivitive(SurfaceType, int, int, double, std::vector<double> &);
	double calculateThetaPerp(int, int);
	double applySurfaceModel(double, double, double, std::vector<double> &);
	std::vector<double> LSSCollectResults(int, int, double, std::vector<double> &, std::vector<double> &, std::vector<bool> &, std::vector<double> &);

	SurfaceType determineSurfaceType(std::vector<std::vector<int>> &, MapTypes);
	std::vector<double> determineModelCoef(SurfaceType, int, std::vector<double> &, std::vector<bool> &, double, double, double);
	std::vector<double> surfaceFit3(std::vector<double> &, std::vector<bool> &, double, double, double);
	std::vector<double> surfaceFit6(std::vector<double> &, std::vector<bool> &, double, double, double);
	std::vector<double> surfaceFit10(std::vector<double> &, std::vector<bool> &, double, double, double);

	bool checkEdgeCriteria(Coordinates, PartitionSet, double, double);
	std::vector<double> determineMappedEdges(double, double, bool, PartitionSet, Coordinates, MapTypes);

	double surfaceModel(std::vector<Scan>, PartitionSet, double, double, double, double);

	void subtractDeviations(std::vector<Survey> &, std::vector<double> &);
};