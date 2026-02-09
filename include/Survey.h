#pragma once
#include <string>
#include <vector>
#include <math.h>
#include "GBParser.h"
#include "FourtyParser.h"
#include "PreProcessor.h"
#include "Scan.h"


class Survey
{
public:

	//constructors
	Survey();
	Survey(SurveyParameters &, SpectralParameters, std::string);
	Survey(SurveyParameters, std::string, std::string);

	void calculateEdgeParameters();

	std::vector<std::vector<double> > gbData; //DYLAN

	//coordinate transforms
	void convertToGalacticInitial(); //DYLAN
	void zeroCrossCheck(); // DYLAN
    // RA/DEC meridian edge case
    void zeroRaCheck();

	//data pre-processing
	void setParams(SurveyParameters);
	void setTwentyParams(GBParser &);    //DYLAN
	void setSdfitsParams(SurveyParameters &, PreProcessor);//DYLAN
	void determineInputFile(std::string); //DYLAN
	void setFortyParams(FourtyParser &); // DYLAN
	void process(bool, bool, bool, double, double);
	void switchChannels(Channel);

	//setters
	void set2DScatterVec(std::vector<double>);
	void setCentroidLocations(std::vector<double> &);
	void setClassificationsLSS(std::vector<std::vector<std::vector<int> > > &);
	void setClassificationsSSS(std::vector<std::vector<std::vector<int> > > &);
	void setEdgeRadius(double);
	void setLSS2DScatterVec(std::vector<double>);
	void setMedianDec(double);
	void setMedianRa(double);
	void setPartSetProcLSS(PartitionSet&);
	void setPartSetProcSSS(PartitionSet&);
	void setRFIScale(double);
	void setMappingCoordinate(std::string);
	void setScans(std::vector<Scan>&);
	void setScansInRa(bool);
	void setStandardThetaGap();
	void setSurveyNumber(int);
	void setTimeShift(double);

	//getters
	bool getScanDirection();
	bool getTracking();
	Channel getChannel();
	Coordinates getProcessingCoordinate();//DYLAN
	std::string getMappingCoordinate();
	double getDiffAlongSweeps();
	double getForcedTS();
	double getMedianDec();
	double getMedianDiffAlongSweeps();
	double getMedianLatiMap(); //DYLAN
	double getMedianLongMap(); //DYLAN
	double getMedianRa();
	double getMinGapThreshold();
	double getMJD();
	double getProcEdgeRadius();
	double getPSFFWHM();
	double getRFIScale();
	double getStandardGap();
	double getTimeShift();
	double getTrimSize();
	int getSurveyNumber();
	MapTypes getMapType();
	PartitionSet getPartSetProcLSS();
	PartitionSet getPartSetProcSSS();
	std::vector<double> getLSS2DScatter();
	std::vector<double> getScatter2d();
	std::vector<Scan> getScans();
	std::vector<std::vector<std::vector<int> > > getClassificationsLSS();
	std::vector<std::vector<std::vector<int> > > getClassificationsSSS();

	virtual ~Survey();

protected:
private:

	//constants
	MapTypes mapType;
	PartitionSet partSetProcSSS;
	PartitionSet partSetProcLSS;
	Channel channel;
	Coordinates pCoordinate;//DYLAN

	//processing information
	double t_int;
	double rfiScale;

	//survey information
	int telescope;
	int numberOfPetals;
	int scanCount;
	int surveyNumber;
	int MJD;
	int yearOfObs;

	bool scansInRa = false;
	bool tracking;
	bool conversion;//DYLAN
	bool debugging;
	bool ASCII;

	double trimSize;
	double forcedTS;
	double psfFWHM;
	double telescopeFrequency;
	double StandardGap;
	double minGapThreshold;
	double petalPeriod;
	double medianDiffAlongSweeps = 0;

	std::vector<std::vector<std::vector<int> > > classificationsSSS;
	std::vector<std::vector<std::vector<int> > > classificationsLSS;

	CalMethods calMethod;
	std::string dataFile;
	std::string interpolationMethod;

	std::vector<double> scatter2d;
	std::vector<double> LSS2DScatter;
	std::vector<double> centroidLocations;

	std::vector<Scan> scans;

	std::vector<std::vector<double> > times;
	std::vector<std::vector<double> > ras;
	std::vector<std::vector<double> > decs;
	std::vector<std::vector<double> > azimuths;
	std::vector<std::vector<double> > elevations;
	std::vector<std::vector<double> > fluxL;
	std::vector<std::vector<double> > fluxR;
	std::vector<std::vector<double> > fluxComp;
	std::vector<std::vector<double> > calibrationFlags;
	std::vector<std::vector<double> > dataDumps;
	std::vector<std::vector<std::vector<int> > > classifications;

	//data pre-processing
	void dataProc(std::vector<std::vector<double> >&);
	void dataProc40(std::vector<std::vector<double> >&);
	void correctDec40(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, int);
	void formatData11(std::vector<std::vector<double> >&);
	void initializeData(std::vector<std::vector<double> >);
	void gainCalibration(Channel, Channel);
	void janskyCalibration(double, Channel);

	//edge calculations
	void determineEdgeFlags();

	//misc
	std::string mCoordinate;
	void daisySweepBreaker();
	void findCenter();
	void daisyPrelim();
	void findCenters();
	void set2DScatter();

	double gauss(double, double);
	void addBackgroundSignal(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addElevationSignal(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addEnRouteDrift(double, std::vector<double>&, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addLongRFI(double, int[], int[], double[], double[], std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addPointSources(int randomScans[], int randomPoints[], std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addPointSourcesCosSq(int randomScans[], int randomPoints[], std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void addShortRFI(double, int[], int[], double[], double[], std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void determineEnRouteDriftRandomNumbers(std::vector<double>&, std::vector<std::vector<double>> &);
	void determineLongRFIRandomNumbers(double, int[], int[], double[], double[], std::vector<std::vector<double>> &);
	void determinePointSourceRandomNumbers(int randomScans[], int randomPoints[], std::vector<std::vector<double>> &);
	void determineShortRFIRandomNumbers(double, int[], int[], double[], double[], std::vector<std::vector<double>> &);
	void insertSimulation(double, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
	void replaceDataWithSimulation(std::vector<std::vector<double>> &);
	void replaceDataWithSimulationGauss(std::vector<std::vector<double>> &);
	void writeFluxToCSV(std::string, std::vector<std::vector<double>> &);
};

