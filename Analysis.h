#include "Map.h"
#include "Scan.h"
#include "Survey.h"
#include "Composite.h"
#include "Cartographer.h"

enum CentroidMethods {CENTER, BRIGHTEST, COORDINATES};

struct PhotoParams
{
	int perform;
	int numberOfSources;
	double innerRadius;
	double outerRadius;
	CentroidMethods centroidType;
	std::vector<double> coordinatesDeg;
	std::vector<double> coordinatesPixels;
};

struct Dimensions
{
	double minRa;
	double maxRa;
	double minDec;
	double maxDec;

	double res;
	double rIn;
	double rOut;
	double rSearch;
};


class Analysis
{
public:
	Analysis(Composite, MapParameters);

	PartitionSet compPartSetProc;
	std::vector<Scan> scans;

	// Callers
	void photometry(Map &, PhotoParams);
	std::vector<double> photometer(int, std::vector<double>, Dimensions);

	// Multi-Threaded Callers
	void photometryMulti(Map &, PhotoParams);

	// Annulus & Aperture
	std::vector<double> annulusCalculations(Dimensions, std::vector<double>);
	std::vector<double> apertureCalculations(Dimensions, std::vector<double>, double);

	// Calculations
	void photometryMetrics(Dimensions dim, std::vector<double>, std::vector<double>, std::vector<double> &);
	std::vector<double> sigmaPlus(double, std::vector<std::vector<double>>);
	std::vector<double> determineCorrFactor(double, double, double, double);

	// Centroid
	std::vector<double> autoCentroid(int, int, double);
	std::vector<double> findBrightestPixel();
	std::vector<double> determineCenters(int, std::vector<double>&, double, double, double);


	~Analysis();

private:
	// Defining Functions
	void definePixelParameters(double &, double &, double &, double &, Map);

	// Determining Functions
	void determinePixelParameters(double &, double &, double &, double &);
	void determineDimensions(std::vector<double>, Dimensions &, PhotoParams);
	double determineAnnulusSigma(double, double);
	double determineAnnulusMu(std::vector<double>, std::vector<double>);
	std::vector<std::vector<double>> determineNRValues(std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<bool>);

	// Map Parameters
	Map map;
	bool M10Processing;
	double weightScale;
	double pixelSize;

	// Miscellaneous
	double psfFWHM;
	double rfiScaleBW;
};

