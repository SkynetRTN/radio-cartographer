#pragma once
#include <vector>

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
	string timeShift;
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