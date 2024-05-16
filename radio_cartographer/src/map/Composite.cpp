#include "map\Composite.h"
#include "utils\RCR.h"
#include "utils\Tools.h"
#include <iostream>
#include <math.h>
#include <iostream>


static double toRad = M_PI / 180.0;
static double toDeg = 180.0 / M_PI;

class LinearModel : public FunctionalForm
{
public:
	double m, b, xBar;

	std::vector<double> w, x, y;
	LinearModel(std::vector<double>, std::vector<double>);
	LinearModel(std::vector<double>, std::vector<double>, std::vector<double>);

	void printData();
	void buildModelSpace();
	std::vector<double> getErrors(std::vector<double>);
	std::vector<double> regression();
	void getAverage();
	~LinearModel();

};

Composite::Composite()
{

}
Composite::Composite(std::vector<Survey> & surveys)
{
	std::vector<MapTypes> mapTypes;
	std::vector<Scan> scansAll;
	std::vector<Scan> scansHold;

	this->minGapThreshold = 999999;
	this->psfFWHM = surveys[0].getPSFFWHM();
	this->rfiScale = surveys[0].getRFIScale();
	this->pCoordinate = surveys[0].getProcessingCoordinate();
	this->mCoordinate = surveys[0].getMappingCoordinate();

	scatter2D.resize(surveys.size());

	for (int i = 0; i < surveys.size(); i++)
	{
		surveyMapTypes.push_back(surveys[i].getMapType());
		standardGaps.push_back(surveys[i].getStandardGap());
		surveys[i].setSurveyNumber(i);
		scatter2D[i] = surveys[i].getScatter2d();

		if (surveys[i].getMinGapThreshold() < minGapThreshold)
		{
			minGapThreshold = surveys[i].getMinGapThreshold();
		}
		scansHold = surveys[i].getScans();

		// Trims the turning point edges.
		if (surveys[0].getTrimSize() != 0.0)
		{
			truncateTurningEdges(surveys[i], scansHold);
		}

		// Collects all partition sets
		partSetVecSSS.push_back(surveys[i].getPartSetProcSSS());
		partSetVecLSS.push_back(surveys[i].getPartSetProcLSS());

		// Converts to processing coordinates
		if (mCoordinate == "equatorial" && pCoordinate == GALACTIC)
		{
			convertToGalactic(scansHold, surveys[i]);
		}
		else if (mCoordinate == "galactic" && pCoordinate == EQUATORIAL)
		{
			convertToEquatorial(scansHold, surveys[i]);
		}

		scansHold[scansHold.size() - 1].setLastScan(true);
		scansAll.insert(scansAll.end(), scansHold.begin(), scansHold.end());
	}
	this->scans = scansAll;
	assignRCRThetaGapMin(surveys);

}

//setters
void Composite::setScans(std::vector<Scan> scansHold)
{
	this->scans = scansHold;
}
void Composite::setCompPartSetProcSSS(PartitionSet partSetHold)
{
	this->compPartSetProcSSS = partSetHold;
}
void Composite::setCompPartSetProcLSS(PartitionSet partSetHold)
{
	this->compPartSetProcLSS = partSetHold;
}
void Composite::setClassificationsSSS(std::vector<std::vector<std::vector<int> > > & classVec)
{
	this->classificationsSSS = classVec;
}
void Composite::setClassificationsLSS(std::vector<std::vector<std::vector<int> > > & classVec)
{
	this->classificationsLSS = classVec;
}



//coordinate transforms
void Composite::convertToGalactic(std::vector<Scan> &scansHold, Survey &survey)
{
	PartitionSet proc = survey.getPartSetProcSSS();

	std::vector<double> holdB, holdL;
	double toB, toL, medianB, medianL;
	
	// == UNDOES THE COSINE TRANSFORM == //
	if (survey.getTracking())
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].undoCosTransform(0, proc.medianLatiMap, proc.medianLongMap);
		}
	}
	else
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].undoCosTransform(0, proc.centerDecDeg, proc.centerRaDeg);
		}
	}
	// == CONVERTS EQUATORIAL TO GALACTIC == //
	for (int i = 0; i < scansHold.size(); i++)
	{
		holdB.clear();
		holdL.clear();
		for (int j = 0; j < scansHold[i].getSize(); j++)
		{
			toB = Tools::convertToB(scansHold[i].getRa(j), scansHold[i].getDec(j));
			toL = Tools::convertToL(scansHold[i].getRa(j), scansHold[i].getDec(j));
			holdB.push_back(toB);
			holdL.push_back(toL);
		}
		scansHold[i].setDec(holdB);
		scansHold[i].setRa(holdL);
	}
	// == COSINE TRANSFORM THE GALACTIC COORDINATES == //
	if (survey.getTracking())
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			medianB = Tools::convertToB(survey.getMedianLongMap(), survey.getMedianLatiMap());
			medianL = Tools::convertToL(survey.getMedianLongMap(), survey.getMedianLatiMap());

			scansHold[i].dynamicCosDecTransform(0, medianB, medianL, scansHold[i].getDec());
			scansHold[i].updateAngDistTemp(medianB);
		}
	}
	else
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].dynamicCosDecTransform(0, proc.centerLatProcDeg, proc.centerLongProcDeg, scansHold[i].getDec());
		}
	}

}
void Composite::convertToEquatorial(std::vector<Scan> &scansHold, Survey &survey)
{
	PartitionSet proc = survey.getPartSetProcSSS();

	std::vector<double> holdDec, holdRa;
	double toDec, toRa, medianDec, medianRa;
	
	// == UNDOES THE COSINE TRANSFORM == //
	if (survey.getTracking())
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].undoDynamicCosDecTransform(0, proc.medianLatiMap, proc.medianLongMap, scansHold[i].getDec());
		}
	}
	else
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].undoDynamicCosDecTransform(0, proc.centerDecDeg, proc.centerRaDeg, scansHold[i].getDec());
		}
	}
	// == CONVERTS GALACTIC TO EQUATORIAL == //
	for (int i = 0; i < scansHold.size(); i++)
	{
		holdDec.clear();
		holdRa.clear();
		for (int j = 0; j < scansHold[i].getSize(); j++)
		{
			toDec = Tools::convertToDec(scansHold[i].getRa(j), scansHold[i].getDec(j));
			toRa  = Tools::convertToRa(scansHold[i].getRa(j), scansHold[i].getDec(j));
			holdDec.push_back(toDec);
			holdRa.push_back(toRa);
		}
		scansHold[i].setDec(holdDec);
		scansHold[i].setRa(holdRa);
	}
	// == COSINE TRANSFORM THE EQUATORIAL COORDINATES == //
	if (survey.getTracking())
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			medianDec = Tools::convertToDec(survey.getMedianLongMap(), survey.getMedianLatiMap());
			medianRa = Tools::convertToRa(survey.getMedianLongMap(), survey.getMedianLatiMap());

			scansHold[i].cosDecTransform(0, medianDec, medianRa, medianDec);
			scansHold[i].updateAngDistTemp(medianDec);
		}
	}
	else
	{
		for (int i = 0; i < scansHold.size(); i++)
		{
			scansHold[i].cosDecTransform(0, proc.centerLatProcDeg, proc.centerLongProcDeg, proc.centerLatProcDeg);
		}
	}

}

//general functions
void Composite::assignRCRThetaGapMin(std::vector<Survey> & surveys)
{
	bool withinBounds = false;
	MapTypes mapTypeHold;

	std::vector<double> edgeOneParameters, edgeTwoParameters, edgeThreeParameters, edgeFourParameters;
	double dec, ra;

	for (int i = 0; i < partSetVecSSS.size(); i++)
	{
		mapType = partSetVecSSS[i].mapType;

		edgeOneParameters = partSetVecSSS[i].edgeOneParameters;
		edgeTwoParameters = partSetVecSSS[i].edgeTwoParameters;
		edgeThreeParameters = partSetVecSSS[i].edgeThreeParameters;
		edgeFourParameters = partSetVecSSS[i].edgeFourParameters;

		if (mapType != DAISY)
		{
			bool withinBoundsRa = false, withinBoundsDec = false;
			double raBoundRight, raBoundLeft, decBoundTop, decBoundBottom;

			double mOne, mTwo, mThree, mFour;
			double bOne, bTwo, bThree, bFour;
			double ybarOne, xbarTwo, ybarThree, xbarFour;

			mOne = edgeOneParameters[0];
			mTwo = edgeTwoParameters[0];
			mThree = edgeThreeParameters[0];
			mFour = edgeFourParameters[0];

			bOne = edgeOneParameters[1];
			bTwo = edgeTwoParameters[1];
			bThree = edgeThreeParameters[1];
			bFour = edgeFourParameters[1];

			ybarOne = edgeOneParameters[2];
			xbarTwo = edgeTwoParameters[2];
			ybarThree = edgeThreeParameters[2];
			xbarFour = edgeFourParameters[2];

			for (int j = 0; j < scans.size(); j++)
			{
				for (int k = 0; k < scans[j].getSize(); k++)
				{
					dec = scans[j].getDec(k);
					ra = scans[j].getRa(k);

					raBoundRight = (mOne * (dec - ybarOne) + bOne);
					raBoundLeft = (mThree * (dec - ybarThree) + bThree);

					decBoundTop = (mTwo * (ra - xbarTwo) + bTwo);
					decBoundBottom = (mFour * (ra - xbarFour) + bFour);

					if (raBoundRight < ra  && ra < raBoundLeft)
					{
						withinBoundsRa = true;
					}

					if (decBoundTop > dec && dec > decBoundBottom)
					{
						withinBoundsDec = true;
					}

					if (withinBoundsRa && withinBoundsDec)
					{
						if (scans[j].getRCRMinThetaGap(k) > surveys[i].getMinGapThreshold())
						{
							scans[j].setMinRCRThetaGap(k, surveys[i].getMinGapThreshold());
						}
					}
				}
			}


		}
		else
		{
			double edgeRadius;
			double distance;
			double medianDec, medianRa;
			edgeRadius = partSetVecSSS[i].edgeRadius;
			medianDec = partSetVecSSS[i].medianDec;
			medianRa = partSetVecSSS[i].medianRa;
				

			for (int j = 0; j < scans.size(); j++)
			{
				for (int k = 0; k < scans[j].getSize(); k++)
				{
					dec = scans[j].getDec(k);
					ra = scans[j].getRa(k);

					distance = Tools::getModGCDistance(dec, ra, medianDec, medianRa)*toDeg;

					if ((distance / psfFWHM) < (edgeRadius / psfFWHM))
					{
						if (scans[j].getRCRMinThetaGap(k) > surveys[i].getMinGapThreshold())
						{
							scans[j].setMinRCRThetaGap(k, surveys[i].getMinGapThreshold());
						}	
					}
				}
			}
		}
	}
	
	this->scans;
}
void Composite::truncateTurningEdges(Survey &survey, std::vector<Scan> &scans)
{
	PartitionSet proc = survey.getPartSetProcSSS();
	std::vector<bool> removalFlags;
	std::vector<double> turningLati, turningLong;
	double nerfRange = survey.getTrimSize()*psfFWHM;
	double distance, i_0, j_0;
	bool stop = true;

	turningLong = proc.edgeLocations[0];
	turningLati = proc.edgeLocations[1];

	for (int k = 0; k < scans.size(); k++)
	{
		i_0 = turningLati[k];
		j_0 = turningLong[k];

		removalFlags.clear();
		removalFlags.resize(scans[k].getSize(), true);

		for (int l = 0; l < scans[k].getSize(); l++)
		{
			distance = Tools::getGCDistance(i_0, j_0, scans[k].getDec(l), scans[k].getRa(l), proc.centerDecDeg)*toDeg;

			if (distance < nerfRange)  // distance is measured from turning point inward.
			{
				removalFlags[l] = false;
			}
		}

		scans[k].removePoints(removalFlags);
		scans[k].updateAngDistTemp(proc.centerDecDeg);

		if (stop && k == scans.size() - 1)  // switches to the other turning edge's parameters.
		{									// I defined daisies to have two turning point edges as well so that this works.
			k = -1;
			turningLong = proc.edgeLocations[2];
			turningLati = proc.edgeLocations[3];
			stop = false;
		}
	}

	setScans(scans);
}




//getters
MapTypes Composite::getMapType()
{
	return mapType;
}
MapTypes Composite::getMapTypes(int index)
{
	return surveyMapTypes[index];
}
Coordinates Composite::getProcessingCoordinate()
{
	return pCoordinate;
}
PartitionSet Composite::getPartSetProcSSS(int index)
{
	return partSetVecSSS[index];
}
PartitionSet Composite::getPartSetProcLSS(int index)
{
	return partSetVecLSS[index];
}
PartitionSet Composite::getCompPartSetProcSSS()
{
	return compPartSetProcSSS;
}
PartitionSet Composite::getCompPartSetProcLSS()
{
	return compPartSetProcLSS;
}
std::vector<PartitionSet> Composite::getPartSetVecSSS()
{
	return partSetVecSSS;
}
std::vector<PartitionSet> Composite::getPartSetVecLSS()
{
	return partSetVecLSS;
}
std::vector<Scan> Composite::getScans()
{
	return scans;
}
std::vector<std::vector<std::vector<int> > > Composite::getClassificationsLSS()
{
	return classificationsLSS;
}
std::vector<std::vector<std::vector<int> > > Composite::getClassificationsSSS()
{
	return classificationsSSS;
}


double Composite::getPSFFWHM()
{
	return psfFWHM;
}
double Composite::getMinGapThreshold()
{
	return minGapThreshold;
}
std::string Composite::getMappingCoordinate()
{
	return this->mCoordinate;
}
std::vector<double> Composite::getStandardGaps()
{
	return standardGaps;
}

std::vector<std::vector<double> > Composite::getScatter2D()
{
	return scatter2D;
}


Composite::~Composite()
{

}