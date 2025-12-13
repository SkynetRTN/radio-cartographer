#pragma once
#include <utility>
#include <vector>
#include <string>
#include "Pixel.h"

enum FitsFile {
        MAIN_SMALL,
        MAIN_LARGE,
        PATH,
        SCALE,
        WEIGHT,
        WEIGHT_TEMP,
        CORRELATION
};

class Map
{
public:
	Map();
	
	//MAP PREPROCESSING
	void reserveGrid(double, double, double, double, double, double, double);

	//MAP PRINTING
	void printSSSProcFluxMap();
	void printSSSWeightMap();
	void printSSSWeight2Map();
	void printSSSCorrelationMap();
	void printSSSScaleMap();

	void printLSSProcFluxMap();
	void printLSSWeightMap();
	void printLSSWeight2Map();
	void printLSSCorrelationMap();
	void printLSSScaleMap();

	void printProcPathMap();
	void printLSSMap();
	void printScanNumberMap();
	void printIndexNumberMap();

        void printFits(const std::string &fitsName, FitsFile fitsType,
                       const std::vector<std::pair<std::string, double>> &headerDoubles = {},
                       const std::vector<std::pair<std::string, std::string>> &headerStrings = {});

	//MAP GETTERS
	double getMaxDec();
	double getMinDec();
	double getMaxRa();
	double getMinRa();
	double getCenterDec();
	double getCenterRa();
	double getResolution();
	
	//GRID GETTERS
	int getSize(int);
	bool getCentroidFlag(int, int);
	double getDec(int, int);
	double getRa(int, int);
	double getLSSFlux(int, int);
	double getProcPath(int, int);

	double getSSSProcFlux(int, int);
	double getSSSScale(int, int);
	double getSSSCorrelation(int, int); //DYLAN
	double getSSSWeight(int, int);
	double getSSSWeight2(int, int);

	double getLSSProcFlux(int, int);
	double getLSSScale(int, int);
	double getLSSCorrelation(int, int); //DYLAN
	double getLSSWeight(int, int);
	double getLSSWeight2(int, int);


	int getScanNumber(int, int);
	int getIndexNumber(int, int);


	//GRID SETTERS
	void setCentroidFlag(int, int, bool);
	void setDec(int, int, double);
	void setRa(int, int, double);
	void setLSSFlux(int, int, double);
	void setProcPath(int, int, double);

	void setSSSProcFlux(int, int, double);
	void setSSSScale(int, int, double);
	void setSSSCorrelation(int, int, double); //DYLAN
	void setSSSWeight(int, int, double);
	void setSSSWeight2(int, int, double);

	void setLSSProcFlux(int, int, double);
	void setLSSScale(int, int, double);
	void setLSSCorrelation(int, int, double); //DYLAN
	void setLSSWeight(int, int, double);
	void setLSSWeight2(int, int, double);
	
	void setScanNumber(int, int, int);
	void setIndexNumber(int, int, int);

	~Map();

private:
	std::vector<std::vector<Pixel> > grid;

	double pixelSize;
	double maxRa;
	double minRa;
	double maxDec;
	double minDec;
	double centerRa;
	double centerDec;
	double resolution;

};

