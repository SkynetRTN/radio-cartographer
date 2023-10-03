#pragma once
#define _USE_MATH_DEFINES
#include <vector>
#include "Structures.h"



class Scan
{
public:
	Scan(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>);

	//functions
	std::vector<bool> determineTSEdgePoints(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, bool);
	void updateAngDist();
	void updateAngDistTemp(double);
	void cosDecTransform(double, double, double, double);
	void undoCosTransform(double, double, double);
	void dynamicCosDecTransform(double, double, double, std::vector<double>);//DYLAN
	void undoDynamicCosDecTransform(double, double, double, std::vector<double>);//DYLAN
	void calculateScatter();
	void calculateElevationScatter();
	void calculateBackground(double);

	//raw data
	bool getCentroidFlag(int);
	bool getLastScan();
	bool getScanInRa();
	bool getRFIFlags(int);
	double getAngDist(int);
	double getDataDumps(int);
	double getDec(int);
	double getDrop2dValues(int);
	double getDumpSum();
	double getElevation(int);
	double getElevationScatter();
	double getFlux(int);
	double getInterScanGap();
	double getIntraScanGap();
	double getSSSCorrelation(int);
	double getLSSData(int);
	double getLSSDec(int);
	double getLSSRa(int);
	double getLSSThetaGap(int);
	double getLSSThetaWPrime(int);
	double getOrigDec(int);//DYLAN
	double getOrigRa(int); //DYLAN
	double getRa(int);
	double getRawDec(int);
	double getRawRa(int);
	double getRCRMinThetaGap(int);
	double getScatter();
	double getGMWeight(int);
	double getThetaGap(int);
	double getTime(int);
	double getTSDec(int);
	double getTSRa(int);
	int getCenter();
	int getDropFlag(int);
	int getEdgePointFlag(int);
	int getExtraLocalModelCount(int);//DYLAN
	int getIntHolder(int);
	int getLocalModelCount(int);//DYLAN
	int getLocalModelInstances(int);
	int getLocalModelRejections(int);
	int getLSSSize();
	int getMaxIndex();
	int getMaxIndexRaw();
	int getMinIndex();
	int getMinIndexRaw();
	int getRawSize();
	int getScanNumberInSurvey();
	int getSize();
	int getSurveyNumber();
	int getTurningPointFlag(int);//DYLAN

	std::vector<double> getFlux();
	std::vector<double> getLFlux();
	std::vector<double> getRFlux();
	std::vector<double> getDataDumps();

	std::vector<double> getAngDist();
	std::vector<double> getDec();
	std::vector<double> getDrop2dValues();
	std::vector<double> getElevation();
	std::vector<double> getLSSData();
	std::vector<double> getLSSDec();
	std::vector<double> getLSSRa();
	std::vector<double> getRa();
	std::vector<double> getRawDec();
	std::vector<double> getRawRa();
	std::vector<double> getThetaGap();
	std::vector<double> getTime();
	std::vector<double> getTSDec();
	std::vector<double> getTSRa();

	void removeBG(std::vector<double>);
	void removeElevation(std::vector<double>);
	void removeRFI(std::vector<double>);
	void removeLSSPoints(std::vector<bool>);
	void removePoints(std::vector<bool>);
	void clipRFI(std::vector<double>);

	void setBG(std::vector<double>);
	void setCenter(int);
	void setCentroidFlag(int, bool);
	void setCentroidFlag(std::vector<bool>);
	void setDec(std::vector<double>);
	void setDrop2dValues(int, double);
	void setDropFlag(int, int);
	void setEdgePointFlag(int, bool);
	void setEdgePointFlag(std::vector<bool>);
	void setElevation(int, double);
	void setElevationScatter(double);
	void setExtraLocalModelCount(int, double);
	void setFlux(int, double);
	void setFlux(std::vector<double>);
	void setInterScanGap(double);
	void setIntHolder(int, int);
	void setIntraScanGap(double);
	void setLastScan(bool);
	void setSSSCorrelation(int, double); //Sets theta_corr.
	void setLocalModelCount(int, double);
	void setLocalModelInstances(int, int);
	void setLocalModelRejections(int, int);
	void setLSSData(int, double);
	void setLSSData(std::vector<double>);
	void setLSSStruct();
	void setLSSThetaGap(std::vector<double>);
	void setLSSThetaGap(int, double);
	void setLSSThetaWPrime(int, double);
	void setMaxIndex(int);
	void setMinIndex(int);
	void setMinRCRThetaGap(int, double);
	void setRa(std::vector<double>);
	void setRFIFlags(std::vector<bool>);
	void setScanInRa(bool);
	void setScanNumberInSurvey(int);
	void setScatter(double);
	void setSize(int);//DYLAN
	void setGMWeight(int, double);
	void setSurveyNumber(int);
	void setThetaGap(int, double);
	void setThetaGap(std::vector<double>);
	void setTSDec(std::vector<double>);
	void setTSRa(std::vector<double>);
	void setTurningPointFlag(std::vector<bool>);//DYLAN
	void switchChannels(Channel);

	virtual ~Scan();
protected:
private:

	struct coordinateStruct
	{
		double globalCenterDec;
		double globalCenterRa;
		std::vector<double> dec;
		std::vector<double> ra;
		std::vector<double> rawDec;
		std::vector<double> rawRa;
		std::vector<double> TSDec;
		std::vector<double> TSRa;
		std::vector<double> workingDec;
		std::vector<double> workingRa;
	} coordinates;

	struct fluxStruct
	{
		std::vector<double> lChannel;
		std::vector<double> rChannel;
		std::vector<double> compChannel;
		std::vector<double> workingChannel;
	} fluxes;

	struct flagStruct
	{
		std::vector<bool> centroidFlag;
		std::vector<bool> edgePointFlag;
		std::vector<bool> rfiFlags;
		std::vector<bool> turningPointFlag;//DYLAN
		std::vector<int> dropFlag;
	} flags;

	struct scanPropertiesStruct
	{
		bool lastScan;
		bool scanInRa;
		bool twoChannels;
		double cleanDumpSum;
		double elevationScatter;
		double interScanGap;
		double intraScanGap;
		double scatter;
		int center;
		int LSSSize;
		int maxIndex;
		int maxIndexRaw;
		int minIndex;
		int minIndexRaw;
		int preRFISize;
		int rawSize;
		int scanNumberInSurvey;
		int size;
		int surveyNumber;
	} scanProperties;

	struct dataPropertiesStruct
	{
		std::vector<double> angDist;
		std::vector<double> time;
		std::vector<double> elevation;
		std::vector<double> background;
		std::vector<double> dataDumps;
	} dataProperties;

	struct photometryPropertiesStruct
	{
		std::vector<int> localModelInstances;
		std::vector<int> localModelRejections;
		std::vector<double> GMWeights;
		std::vector<double> localModelCount;		//DYLAN
		std::vector<double> extraLocalModelCount;   //DYLAN
		std::vector<double> SSSCorrelation;  //DYLAN
		std::vector<double> centroidLocations;
	} photometryProperties;

	std::vector<double> LSSThetaGap;
	std::vector<double> LSSData;
	std::vector<double> LSSRa;
	std::vector<double> LSSDec;
	std::vector<double> LSSDataDumps;
	std::vector<double> LSSCorrelation;
	std::vector<double> LSSSummedWeights;
	std::vector<double> LSSThetaWPrime;
	std::vector<double> LSSLocalModelCorrelation;

	std::vector<double> minRCRThetaGap;
	std::vector<double> minRCRThetaGapRaw;
	std::vector<double> thetaGap;
	std::vector<double> postRFIThetaGap;
	std::vector<double> preRFIThetaGap;

	std::vector<int> intHolder;

	std::vector<double> drop2DValues;

	void pointToPointDiff(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);
};

