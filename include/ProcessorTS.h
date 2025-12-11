#pragma once
#include <vector>
#include "Scan.h"


class ProcessorTS
{
public:
	ProcessorTS();
	ProcessorTS(std::vector<Scan> &, MapTypes, double, double, bool);
	double find_tInt();
	~ProcessorTS();

private:
	std::vector<Scan> scans;
	PartitionSet partSetProcSSS;
	MapTypes mapType;
	double psfFWHM;
	double medianDiffAlongSweeps;
	bool scansInRa;
	std::vector<bool> shiftRejection(std::vector<double>);
	std::vector<std::vector<double>> getScanToScanShifts(double);
	double calculateShift(std::vector<std::vector<double>>&);



};



