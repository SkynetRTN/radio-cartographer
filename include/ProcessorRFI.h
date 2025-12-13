#pragma once
#include <vector>
#include <future>
#include "Scan.h"

class ProcessorRFI
{
public:
	ProcessorRFI(std::vector<Scan> &, RFIParameters &,  std::vector<std::vector<double>>&, std::vector<std::vector<std::vector<int>>> &);
	std::vector<double> getCentroidLocations();
	std::vector<std::vector<double>> getRFISubtracted();
	~ProcessorRFI();

private:
	double correlatedWeightMap;
	double psfFWHM;
	double rfiScaleBW;
	double rfiScaleDeg;
	double standardGap;
	PartitionSet partSetProcSSS;
	std::vector<double> centroidLocations;
	std::vector<Scan> scans;
	std::vector<std::vector<double> > rfiSubtracted;
	std::vector<std::vector<double>> scatter2dSSS;
	std::vector<std::vector<std::vector<double> > > rfiValues, rfiWeights, rfiCounts, rfiWeights2, rfiNs, rfiDists;
	std::vector<std::vector<std::vector<int>>> classificationsSSS;

	//dylan
	void setGlobalModelValues(std::vector<Scan> &, std::vector<Scan>);

	int findPossSSS(double, double, double, std::vector<int>&);
	std::vector<double> determineCenters(int, std::vector<double>&, double, double, double);
	std::vector<double> rfiCollectResults(double, double, std::vector<double> &, std::vector<int> &, std::vector<double> &);
	std::vector<double> rfiFit(std::vector<int> &, std::vector<double> &);
	std::vector<double> rfiRemovalMulti(double, double, double, double);
	std::vector<std::vector<double>> rfiBuildGlobalModels();
	void autoCentroid(std::vector<std::vector<std::vector<double> > >);
	void rfiBuildExtraModels(double, double);
	void rfiBuildLocalModels();
	void rfiRemovePoints(double, std::vector<int> &, std::vector<double> &);
	void rfiSetFinalValues(std::vector<std::future<std::vector<double>>> &);
};


