#pragma once
#include "Survey.h"
#include "Composite.h"
#include "Map.h"


class Cartographer
{
public:

	Cartographer();
	Cartographer(Composite, MapParameters);

	Survey survey;

	Map procMap;

	double processingWeightingFunction;
	double resolution;
	double psfFWHM;
	bool m10PlusProcessing;
	bool SSSMapping;
	bool LSSMapping;
	bool correlatedWeightMap;
	double rfiScaleDeg;//DYLAN
	std::vector<double> standardGaps;//DYLAN
	std::vector<double> centroidCenters;//DYLAN
	std::string mCoordinate;
	Coordinates pCoordinate;//DYLAN

	PartitionSet compPartSetProcSSS;
	PartitionSet compPartSetProcLSS;
	std::vector<PartitionSet> partSetVecSSS;
	std::vector<PartitionSet> partSetVecLSS;

	std::vector<std::vector<std::vector<int> > > classificationsSSS;
	std::vector<std::vector<std::vector<int> > > classificationsLSS;

	std::vector<Scan> scans;

	Map generateProcMapsMulti();


	//pixel interpolation
	std::vector<double> processedGridInterpolationMulti(double, double);

	std::vector<double> collectPointsSSS(double, double);
	std::vector<double> collectPointsLSS(double, double);
	std::vector<std::vector<int> > determineScanNumbersSSS(std::vector<double> &);
	std::vector<std::vector<int> > determineScanNumbersLSS(std::vector<double> &);


	double getPixelThetaGapSSS(double, double, std::vector<double>&);
	double getPixelThetaGapLSS(double, double, std::vector<double>&);
	double getCorrMapValueSSS(double, double, std::vector<double>&, double);
	double getCorrMapValueLSS(double, double, std::vector<double>&, double);

	//edge calculations
	bool checkEdgeCriteria(double, double);
	std::vector<double> determineMappedEdges(double, double, bool, PartitionSet, Coordinates, MapTypes);//DYLAN

	//surface modeling
	std::vector<double> surfaceFit3(bool, std::vector<double>&, double, double, double);
	std::vector<double> surfaceFit6(bool, std::vector<double>&, double, double, double);
	std::vector<double> surfaceFit10(bool, std::vector<double>&, double, double, double);


	//misc
	SurfaceType determineSurfaceType(std::vector<std::vector<int>> &, MapTypes);
	std::vector<double> applySurfaceModel(SurfaceType, bool, bool, std::vector<double>&, double, double, double);
	int findPossSSS(double, double, double, std::vector<int>&);
	int findPossLSS(double, double, double, std::vector<int>&);
	bool planeApplication(int, int, int, std::vector<std::vector<int>>, int, MapTypes);
	void setMapValues(int, int, std::vector<double> &);


	~Cartographer();
};

