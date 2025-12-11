#include <vector>
#include "Scan.h"
class ProcessorThetaGap
{
public:
	ProcessorThetaGap();
	ProcessorThetaGap(std::vector<Scan> & scans, double psfFWHM, PartitionSet & partSetSSS, PartitionSet & partSetLSS, std::vector<std::vector<std::vector<int> > > & classSSS, std::vector<std::vector<std::vector<int> > > &classLSS);
	std::vector<std::vector<double>> calculateThetaGapSSS();
	std::vector<std::vector<double>> calculateThetaGapLSS();

	~ProcessorThetaGap();

private:
	double psfFWHM;
	PartitionSet partSetProcSSS;
	PartitionSet partSetProcLSS;
	std::vector<Scan> scans;
	std::vector<std::vector<std::vector<int> > > classificationsSSS;
	std::vector<std::vector<std::vector<int> > > classificationsLSS;

	int findPossSSS(double, double, double, std::vector<int> &);
	int findPossLSS(double, double, double, std::vector<int> &);

	std::vector<int> quadrantSort(double, double, std::vector<int>&, Quadrant);
	void determineCircleParams(double, double, Quadrant, std::vector<int>&, std::vector<double> &);
	//void safeguardEdgeThetaGaps(Composite &, double);
	double circleThetaGap(int, int, bool);
	double maxGapQuadrant(std::vector<double>&);
};

