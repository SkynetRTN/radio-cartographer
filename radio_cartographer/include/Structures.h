#pragma once
#include <vector>
#include <string>
#include <type_traits>

enum Receiver2 {
	HIGH,
	LOW,
};
enum Channel
{
	LEFT,
	RIGHT,
	COMPOSITE
};
enum CalMethods
{
	PRE,
	POST,
	INTERPOLATED
};
enum Coordinates
{
	EQUATORIAL,
	GALACTIC
};
enum MapTypes
{
	RASTER,
	NODDING,
	DAISY
};
enum SurfaceType
{
	M10,
	M6,
	M3,
	M0
};
enum Telescopes
{
	FOURTY_FOOT,
	TWENTY_METER,
	GBT
};
enum Quadrant
{
	TOP,
	BOTTOM,
	LEFT_QUAD,
	RIGHT_QUAD,
	DIAG_TOP_RIGHT,
	DIAG_BOTTOM_RIGHT,
	DIAG_TOP_LEFT,
	DIAG_BOTTOM_LEFT
};

enum CentroidMethods 
{ 
	CENTER, 
	BRIGHTEST,
	COORDINATES 
};

struct PartitionSet
{
	// THE PARTITION SET STORES INFORMATION ABOUT A SURVEY'S GEOMETRY WHILE ALSO
	// PARTITIONING DATA POINTS ACCORDING TO SAID GEOMETRY

	// THERE ARE TWO PARTITION SETS FOR EACH SURVEY: PROCESSED (PROC) AND RAW

	// partSetProc ONLY REFERS TO THE GEOMETRY OF THE PROCESSED DATA

	MapTypes mapType;

	double subDecInc; // = subDecInc;
	double subRaInc;  // = subRaInc;
	double subDecRes; // = subDecRes;
	double subRaRes;  // = subRaRes;

	double maxDec; // = maxDec;
	double maxRa;  // = maxRa;
	double minDec; // = minDec;
	double minRa;  // = minRa;

	// ORIGINAL DEGREE VALUES OF THE CENTER OF THE OBSERVATION [NOT INTERNAL COORDINATES]
	// USED FOR APPLYING AND INVERTING THE COSINE TRANSFORM
	double centerDecDeg;
	double centerRaDeg;

	double centerLatProcDeg;  // DYLAN
	double centerLongProcDeg; // DYLAN

	// MEDIAN COORDINATES ARE MEASURED BY PERFORMING RCR ON ALL RA AND DEC VALUES
	// THEY REFLECT THE OBSERVATION'S CENTER IN INTERNAL COORDINATES CENTERED ABOUT THE MAP'S (0,0) COORDINATE
	double medianDec;
	double medianRa;

	double medianLatiMap; // DYLAN
	double medianLongMap; // DYLAN

	std::vector<double> edgeOneParameters;
	std::vector<double> edgeTwoParameters;
	std::vector<double> edgeThreeParameters;
	std::vector<double> edgeFourParameters;

	std::vector<std::vector<double>> edgeLocations; // DYLAN

	bool tracking; // DYLAN
	double trimSize;
	double edgeRadius;
};

struct PreProcessingParameters
{
	double velocity = 0.0;

	double subtractionScale = 0.0;
	double modifiedSubtractionScale = 0.0;;

	std::vector<double> exclusionBand;
	std::vector<std::pair<double, double>> inclusionBands;
	std::vector<std::pair<double, double>> modifiedSubtractionBands;

	Receiver2 receiver;
};

struct PhotometryParams
{
	int perform;
	int numberOfSources;
	double innerRadius;
	double outerRadius;
	CentroidMethods centroidType;
	std::vector<double> coordinatesDeg;
	std::vector<double> coordinatesPixels;
};

struct SurveyParameters
{
	bool tracking;
	bool debug;

	double forcedTS;
	double frequency;
	double trimSize;
	double weightScale;

	Channel channel;
	Telescopes tele;
	MapTypes mapType;
	CalMethods calMethod;
	Coordinates pCoordinate;
};

struct ProcessorParameters
{
	double rfiScaleBW;
	double bgScaleBW;
	double wScaleBW;
	std::string timeShift;
	double timeShiftValue;
	bool wCorrMap;
	bool lssProc;
	int raw;
};

struct MapParameters
{
	bool LSSMapping;
	bool SSSMapping;
	bool standardGapWScale;
	bool rawProcessing;
	bool m10PlusProcessing;
	bool photometry;
	bool correlatedWeightMap;
	double pixelSize;
	double rawWeightScale;
	double processedWeightScale;
	double edgeThreshold;
	double rfiScale;
};

struct RFIParameters
{
	bool correlatedWeightMap;
	double rfiScaleBW;
	double psfFWHM;
	double standardGap;
	PartitionSet partSetProcSSS;
};

//template<typename Enum>
//Enum mapStringToEnum(const std::string& enumStr) {
//	//std::transform(enumStr.begin(), enumStr.end(), enumStr.begin(), ::toupper);
//
//	if (std::is_same<Enum, CentroidMethods>::value) {
//		if (enumStr == "CENTER") return Enum::CENTER;
//		if (enumStr == "BRIGHTEST") return Enum::BRIGHTEST;
//		if (enumStr == "COORDINATES") return Enum::COORDINATES;
//	}
//	else if (std::is_same<Enum, Coordinates>::value) {
//		if (enumStr == "EQUATORIAL") return Enum::EQUATORIAL;
//		if (enumStr == "GALACTIC") return Enum::GALACTIC;
//	}
//	else if (std::is_same<Enum, CalMethods>::value) {
//		if (enumStr == "PRE") return Enum::PRE;
//		if (enumStr == "POST") return Enum::POST;
//		if (enumStr == "INTERPOLATED") return Enum::INTERPOLATED;
//	}
//	else if (std::is_same<Enum, Channel>::value) {
//		if (enumStr == "left") return Enum::LEFT;
//		if (enumStr == "right") return Enum::RIGHT;
//		if (enumStr == "composite") return Enum::COMPOSITE;
//	}
//
//	throw "Invalid enum string parsed from configuration file: " + enumStr;
//}
