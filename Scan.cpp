#include "Scan.h"
#include "Tools.h"
#include "RCR.h"
#include <iostream>
#include <cmath>
#include <math.h>

static double toRad = M_PI / 180.0;
static double toDeg = 180.0 / M_PI;

Scan::Scan(std::vector<double> time, std::vector<double> dec, std::vector<double> ra, std::vector<double> elevations, std::vector<double> dataDumps, std::vector<double> lChannel, std::vector<double> rChannel, std::vector<double> compChannel)
{
	this->scanProperties.twoChannels = false;
	this->scanProperties.lastScan = false;
	this->scanProperties.size = time.size();
	this->scanProperties.rawSize = ra.size();
	
	this->dataProperties.time = time;
	this->dataProperties.dataDumps = dataDumps;
	this->dataProperties.elevation = elevations;
	
	this->coordinates.dec = dec;
	this->coordinates.ra = ra;
	this->coordinates.workingDec = dec;
	this->coordinates.workingRa = ra;

	this->fluxes.lChannel = lChannel;
	this->fluxes.rChannel = rChannel;
	this->fluxes.compChannel = compChannel;

	std::vector<int> intFiller, intFiller2;
	std::vector<bool> boolFiller;
	std::vector<double> dubFiller, dubFiller2, dubFiller3, dubFiller4;
	intFiller.resize(time.size(), 0);
	intFiller2.resize(time.size(), -1);
	boolFiller.resize(time.size(), true);
	dubFiller.resize(time.size(), -999999);
	dubFiller2.resize(time.size(), 0.0);
	dubFiller3.resize(time.size(), 999999);
	dubFiller4.resize(time.size(), 1.0);

	this->photometryProperties.localModelInstances = intFiller;
	this->photometryProperties.localModelRejections = intFiller;
	this->photometryProperties.GMWeights = dubFiller4;
	this->photometryProperties.localModelCount = dubFiller2;//Dylan
	this->photometryProperties.extraLocalModelCount = dubFiller2;//Dylan
	this->photometryProperties.SSSCorrelation = dubFiller4;//Dylan

	this->minRCRThetaGap = dubFiller3;
	this->minRCRThetaGapRaw = dubFiller3;
	this->drop2DValues = dubFiller;
	this->thetaGap = dubFiller;
	this->flags.rfiFlags = boolFiller;
	this->flags.edgePointFlag = boolFiller;
	this->intHolder = intFiller;
	this->flags.dropFlag = intFiller2;

	fluxes.workingChannel = compChannel;
	//updateAngDistTemp(0);

	//raw map data
	
	this->coordinates.rawDec = dec;
	this->coordinates.rawRa = ra;

	this->LSSRa = ra;
	this->LSSDec = dec;
	this->LSSThetaGap = dubFiller3;
	this->LSSCorrelation = dubFiller3;
	this->LSSThetaWPrime = dubFiller3;
}

class PointToPointFunc : public NonParametric
{
public:
	PointToPointFunc(std::vector<double>, std::vector<double>);

	void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);
	void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::vector<double> workingChannel, angDist;
	std::vector<int> indices;
	~PointToPointFunc();

};

PointToPointFunc::PointToPointFunc(std::vector<double> workingChannel, std::vector<double> angDist)
{
	this->workingChannel = workingChannel;
	this->angDist = angDist;
}
void PointToPointFunc::muFunc(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &y, std::vector<double> &trueY)
{
	int trueCount = 0, currentIndex, low, high;
	std::vector<int> indicesVec;
	std::vector<bool> badModels;
	std::vector<double> trueYVec, diff;

	diff.resize(y.size(), 0.0);
	badModels.resize(y.size(), true);
	for (int k = 0; k < y.size(); k++)
	{
		low = Tools::max(0, k - 1);
		high = Tools::min(y.size() - 1, k + 1);
		while (low > -1 && !flags[low])
		{
			low--;
		}
		while (high < y.size() && !flags[high])
		{
			high++;
		}
		if (low == -1 || high == y.size() || high - low == 1 || angDist[high] == angDist[low])
		{
			badModels[k] = false;
			diff[k] = 0;
		}
		else
		{
			diff[k] = workingChannel[k] - (workingChannel[low] + ((workingChannel[high] - workingChannel[low]) / (angDist[high] - angDist[low]))*(angDist[k] - angDist[low]));
		}
	}
	for (int i = 0; i < y.size(); i++)
	{
		if (flags[i] && badModels[i])
		{
			trueYVec.push_back(diff[i]);
			indicesVec.push_back(i);
		}
	}
	indices = indicesVec;
	trueY = trueYVec;

}
void PointToPointFunc::muFunc(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &w, std::vector<double> &y, std::vector<double> &trueW, std::vector<double> &trueY)
{
	int trueCount = 0, currentIndex, low, high;
	double xBar, weight;
	std::vector<int> indicesVec;
	std::vector<bool> badModels;
	std::vector<double> trueYVec, trueWVec, diff, weights;

	diff.resize(y.size(), 0.0);
	badModels.resize(y.size(), true);
	weights.resize(y.size(), 0.0);
	for (int k = 0; k < y.size(); k++)
	{
		low = Tools::max(0, k - 1);
		high = Tools::min(y.size() - 1, k + 1);
		while (low > -1 && !flags[low])
		{
			low--;
		}
		while (high < y.size() && !flags[high])
		{
			high++;
		}
		if (low == -1 || high == y.size() || high - low == 1 || angDist[high] == angDist[low])
		{
			badModels[k] = false;
			diff[k] = 0;
			weights[k] = 0;
		}
		else
		{
			diff[k] = workingChannel[k] - (workingChannel[low] + ((workingChannel[high] - workingChannel[low]) / (angDist[high] - angDist[low]))*(angDist[k] - angDist[low]));
			xBar = (w[low] * angDist[low] + w[high] * angDist[high]) / (w[low] + w[high]);
			weight = (1.0 / w[k]);
			weight += pow(((angDist[high] - xBar) * (angDist[high] - angDist[low])), 2.0) / w[low];
			weight += pow(((angDist[low] - xBar) * (angDist[high] - angDist[low])), 2.0) / w[high];
			weight += pow(((angDist[k] - xBar) * (angDist[high] - angDist[low])), 2.0) * ((1.0 / w[high]) + (1.0 / w[low]));
			weights[k] = (1.0 / weight);
		}
	}
	for (int i = 0; i < y.size(); i++)
	{
		if (flags[i] && badModels[i])
		{
			trueYVec.push_back(diff[i]);
			indicesVec.push_back(i);
			trueWVec.push_back(w[i]);
		}
	}
	indices = indicesVec;
	trueY = trueYVec;
	trueW = trueWVec;

}
PointToPointFunc::~PointToPointFunc()
{

}

//operations
void Scan::pointToPointDiff(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &y, std::vector<double> &trueY)
{
	int trueCount = 0, currentIndex, low, high;
	std::vector<int> indicesVec;
	std::vector<bool> badModels;
	std::vector<double> trueYVec, diff;

	diff.resize(y.size(), 0.0);
	badModels.resize(y.size(), true);
	for (int k = 0; k < y.size(); k++)
	{
		low = Tools::max(0, k - 1);
		high = Tools::min(y.size() - 1, k + 1);
		while (low > -1 && !flags[low])
		{
			low--;
		}
		while (high < y.size() && !flags[high])
		{
			high++;
		}
		if (low == -1 || high == y.size() || high - low == 1 || dataProperties.angDist[high] == dataProperties.angDist[low])
		{
			badModels[k] = false;
			diff[k] = 0;
		}
		else
		{
			diff[k] = fluxes.workingChannel[k] - (fluxes.workingChannel[low] + ((fluxes.workingChannel[high] - fluxes.workingChannel[low]) / (dataProperties.angDist[high] - dataProperties.angDist[low]))*(dataProperties.angDist[k] - dataProperties.angDist[low]));
		}
	}
	for (int i = 0; i < y.size(); i++)
	{
		if (flags[i] && badModels[i])
		{
			trueYVec.push_back(diff[i]);
			indicesVec.push_back(i);
		}
	}
	indices = indicesVec;
	trueY = trueYVec;
}
void Scan::calculateScatter()
{
	bool newOne = true;
	std::ofstream OutputFile;
	if (newOne)
	{

		double sum = 0;
		RCR rcr = RCR(SS_MEDIAN_DL);
		PointToPointFunc ptpf(fluxes.workingChannel, dataProperties.angDist);

		rcr.setNonParametricModel(ptpf);
		rcr.performRejection(dataProperties.dataDumps, fluxes.workingChannel);

		// RCR has ~0.2% chance of a runaway rejection which results in a nan sigma.
		// We remove these points (In Processor) before making the model
		this->scanProperties.scatter = .8197*rcr.result.sigma;
		for (int i = 0; i < rcr.result.cleanW.size(); i++)
		{
			sum += rcr.result.cleanW[i];
		}
		this->scanProperties.cleanDumpSum = sum;
	}
	else
	{
		std::vector<double> tempArray, smooth, diff;
		std::vector<bool> checks, checks2;
		int low, high, maxIndex, counter;
		double mu, sigma, max = -999999;
		bool stop = false;
		checks.resize(scanProperties.size, true);
		checks2.resize(scanProperties.size, true);
		diff.resize(scanProperties.size);
		smooth.resize(scanProperties.size);
		while (!stop)
		{
			for (int k = 0; k < scanProperties.size; k++)
			{
				low = Tools::max(0, k - 1);
				high = Tools::min(scanProperties.size - 1, k + 1);
				if (!checks[low] || !checks[high] || high - low == 1 || dataProperties.angDist[high] == dataProperties.angDist[low])
				{
					checks2[k] = false;
				}
				smooth[k] = fluxes.workingChannel[low] + ((fluxes.workingChannel[high] - fluxes.workingChannel[low]) / (dataProperties.angDist[high] - dataProperties.angDist[low]))*(dataProperties.angDist[k] - dataProperties.angDist[low]);
				diff[k] = (fluxes.workingChannel[k] - smooth[k]);
			}

			max = -999999;
			for (int k = 0; k < scanProperties.size; k++)
			{
				if (checks[k] && checks2[k] && std::abs(diff[k]) > max)
				{
					max = std::abs(diff[k]);
					maxIndex = k;
				}
			}
			counter = 0;
			for (int i = 0; i < scanProperties.size; i++)
			{
				if (checks[i] && checks2[i])
				{
					tempArray.push_back(diff[i]);
					counter++;
				}
			}
			mu = Tools::getMedian(tempArray);
			tempArray.clear();
			for (int i = 0; i < scanProperties.size; i++)
			{
				if (checks[i] && checks2[i])
				{
					tempArray.push_back(std::abs(diff[i] - mu));
				}
			}
			sigma = Tools::get68th(tempArray);
			tempArray.clear();
			if (counter*Tools::erfc(std::abs(diff[maxIndex] - mu) / sigma) < .5)
			{
				checks[maxIndex] = false;
			}
			else
			{
				stop = true;
			}
		}
		this->scanProperties.scatter = sigma;
	}
}
void Scan::calculateElevationScatter()
{
	bool newOne = true;
	std::ofstream OutputFile;
	if (newOne)
	{

		double sum = 0;
		RCR rcr = RCR(SS_MEDIAN_DL);
		//PointToPointFunc ptpf(LSSData, elevation);
		PointToPointFunc ptpf(LSSData, dataProperties.elevation);

		rcr.setNonParametricModel(ptpf);
		rcr.performRejection(dataProperties.dataDumps, LSSData);
		this->scanProperties.scatter = .8197*rcr.result.sigma;
		for (int i = 0; i < rcr.result.cleanW.size(); i++)
		{
			sum += rcr.result.cleanW[i];
		}
		this->scanProperties.cleanDumpSum = sum;
	}
	else
	{
		std::vector<double> tempArray, smooth, diff;
		std::vector<bool> checks, checks2;
		int low, high, maxIndex, counter;
		double mu, sigma, max = -999999;
		bool stop = false;
		checks.resize(scanProperties.size, true);
		checks2.resize(scanProperties.size, true);
		diff.resize(scanProperties.size);
		smooth.resize(scanProperties.size);
		while (!stop)
		{
			for (int k = 0; k < scanProperties.size; k++)
			{
				low = Tools::max(0, k - 1);
				high = Tools::min(scanProperties.size - 1, k + 1);
				if (!checks[low] || !checks[high] || high - low == 1 || dataProperties.elevation[high] == dataProperties.elevation[low])
				{
					checks2[k] = false;
				}
				smooth[k] = LSSData[low] + ((LSSData[high] - LSSData[low]) / (dataProperties.elevation[high] - dataProperties.elevation[low]))*(dataProperties.elevation[k] - dataProperties.elevation[low]);
				diff[k] = (LSSData[k] - smooth[k]);
			}

			max = -999999;
			for (int k = 0; k < scanProperties.size; k++)
			{
				if (checks[k] && checks2[k] && std::abs(diff[k]) > max)
				{
					max = std::abs(diff[k]);
					maxIndex = k;
				}
			}
			counter = 0;
			for (int i = 0; i < scanProperties.size; i++)
			{
				if (checks[i] && checks2[i])
				{
					tempArray.push_back(diff[i]);
					counter++;
				}
			}
			mu = Tools::getMedian(tempArray);
			tempArray.clear();
			for (int i = 0; i < scanProperties.size; i++)
			{
				if (checks[i] && checks2[i])
				{
					tempArray.push_back(std::abs(diff[i] - mu));
				}
			}
			sigma = Tools::get68th(tempArray);
			tempArray.clear();
			if (counter*Tools::erfc(std::abs(diff[maxIndex] - mu) / sigma) < .5)
			{
				checks[maxIndex] = false;
			}
			else
			{
				stop = true;
			}
		}
		this->scanProperties.elevationScatter = sigma;
	}
}
void Scan::switchChannels(Channel newChannel)
{
	switch (newChannel)
	{
	case LEFT:
		fluxes.workingChannel = fluxes.lChannel;
		break;
	case RIGHT:
		fluxes.workingChannel = fluxes.rChannel;
		break;
	case COMPOSITE:
		fluxes.workingChannel = fluxes.compChannel;
		break;
	}
}

//transformation
void Scan::updateAngDist()
{
	this->dataProperties.angDist.clear();
	this->dataProperties.angDist.resize(coordinates.workingDec.size(), 0.0);
	for (int i = 1; i < dataProperties.angDist.size(); i++)
	{

		this->dataProperties.angDist[i] = acos(Tools::min(sin(this->coordinates.workingDec[i - 1] * toRad)*sin(this->coordinates.workingDec[i] * toRad) + cos(this->coordinates.workingDec[i - 1] * toRad)*cos(this->coordinates.workingDec[i] * toRad)*cos((this->coordinates.workingRa[i - 1] - this->coordinates.workingRa[i])*toRad), 1.0)) / toRad + dataProperties.angDist[i - 1];
		//this->angDist[i] = Tools::getGCDistance(this->workingDec[i], this->workingRa[i], this->workingDec[i - 1], this->workingRa[i - 1])*toDeg - angDist[i - 1];
		//this->rawDec[i - 1] * toRad;
		//this->rawRa[i - 1] * toRad;
	}
}
void Scan::updateAngDistTemp(double centerDec)
{
	this->dataProperties.angDist.clear();
	this->dataProperties.angDist.resize(coordinates.workingDec.size(), 0.0);
	for (int i = 1; i < dataProperties.angDist.size(); i++)
	{

		//this->dataProperties.angDist[i] = acos(Tools::min(sin(this->workingDec[i - 1] * toRad)*sin(this->workingDec[i] * toRad) + cos(this->workingDec[i - 1] * toRad)*cos(this->workingDec[i] * toRad)*cos((this->workingRa[i - 1] - this->workingRa[i])*toRad), 1.0)) / toRad + dataProperties.angDist[i - 1];
		this->dataProperties.angDist[i] = Tools::getGCDistance(this->coordinates.workingDec[i], this->coordinates.workingRa[i], this->coordinates.workingDec[i - 1], this->coordinates.workingRa[i - 1], centerDec)*toDeg + dataProperties.angDist[i - 1];
		//this->rawDec[i - 1] * toRad;
		//this->rawRa[i - 1] * toRad;
	}
}
void Scan::cosDecTransform(double t_int, double centerDec, double centerRa, double angDistCenterDec)
{
	std::vector<double> decHold, raHold, timeHold;
	std::vector<bool> edgeFlagHold;
	decHold.resize(scanProperties.size);
	raHold.resize(scanProperties.size);
	timeHold.resize(scanProperties.size);
	edgeFlagHold.resize(scanProperties.size, 0);
	int jMin, jMax;
	for (int j = 0; j < scanProperties.size; j++)
	{
		if (t_int >= 0)
		{
			jMin = j;
			while ((dataProperties.time[jMin] >= dataProperties.time[j] - t_int) && jMin > 0)
			{
				jMin = jMin - 1;
			}
			jMax = jMin + 1;
		}
		else
		{
			jMax = j;
			while ((dataProperties.time[jMax] <= dataProperties.time[j] - t_int) && jMax < scanProperties.size - 1)
			{
				jMax = jMax + 1;
			}
			jMin = jMax - 1;
		}

		decHold[j] = (coordinates.workingDec[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingDec[jMax] - coordinates.workingDec[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin])) - centerDec;
		//raHold[j] = (workingRa[jMin] + (time[j] - t_int - time[jMin])*(workingRa[jMax] - workingRa[jMin]) / (time[jMax] - time[jMin]))*cos(centerDec*M_PI / 180.0);
		raHold[j] = ((coordinates.workingRa[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingRa[jMax] - coordinates.workingRa[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin])) - centerRa)*cos(centerDec*M_PI / 180.0);
	}

	this->coordinates.workingDec = decHold;
	this->coordinates.workingRa = raHold;

	//updateAngDistTemp(angDistCenterDec);
}
void Scan::undoCosTransform(double t_int, double centerDec, double centerRa)
{
	std::vector<double> decHold, raHold, timeHold;
	std::vector<bool> edgeFlagHold;
	decHold.resize(scanProperties.size);
	raHold.resize(scanProperties.size);
	timeHold.resize(scanProperties.size);
	edgeFlagHold.resize(scanProperties.size, 0);
	int jMin, jMax;
	for (int j = 0; j < scanProperties.size; j++)
	{
		if (t_int >= 0)
		{
			jMin = j;
			while ((dataProperties.time[jMin] >= dataProperties.time[j] - t_int) && jMin > 0)
			{
				jMin = jMin - 1;
			}
			jMax = jMin + 1;
		}
		else
		{
			jMax = j;
			while ((dataProperties.time[jMax] <= dataProperties.time[j] - t_int) && jMax < scanProperties.size - 1)
			{
				jMax = jMax + 1;
			}
			jMin = jMax - 1;
		}

		decHold[j] = (coordinates.workingDec[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingDec[jMax] - coordinates.workingDec[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin])) + centerDec;
		//raHold[j] = (workingRa[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(workingRa[jMax] - workingRa[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin]))*cos(centerDec*M_PI / 180.0);
		raHold[j] = ((coordinates.workingRa[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingRa[jMax] - coordinates.workingRa[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin]))) / cos(centerDec*M_PI / 180.0) + centerRa;
	}

	this->coordinates.workingDec = decHold;
	this->coordinates.workingRa = raHold;

}
void Scan::dynamicCosDecTransform(double t_int, double CENTERLATI, double CENTERLONG, std::vector<double> LATITUDES)
{
	std::vector<double> bHold, lHold, timeHold;
	std::vector<bool> edgeFlagHold;
	bHold.resize(scanProperties.size);
	lHold.resize(scanProperties.size);
	timeHold.resize(scanProperties.size);
	edgeFlagHold.resize(scanProperties.size, 0);
	int jMin, jMax;
	for (int j = 0; j < scanProperties.size; j++)
	{
		if (t_int >= 0)
		{
			jMin = j;
			while ((dataProperties.time[jMin] >= dataProperties.time[j] - t_int) && jMin > 0)
			{
				jMin = jMin - 1;
			}
			jMax = jMin + 1;
		}
		else
		{
			jMax = j;
			while ((dataProperties.time[jMax] <= dataProperties.time[j] - t_int) && jMax < scanProperties.size - 1)
			{
				jMax = jMax + 1;
			}
			jMin = jMax - 1;
		}

		bHold[j] = (coordinates.workingDec[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingDec[jMax] - coordinates.workingDec[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin])) - CENTERLATI;
		lHold[j] = ((coordinates.workingRa[jMin] + (dataProperties.time[j] - t_int - dataProperties.time[jMin])*(coordinates.workingRa[jMax] - coordinates.workingRa[jMin]) / (dataProperties.time[jMax] - dataProperties.time[jMin])) - CENTERLONG)*cos((LATITUDES[jMax])*M_PI / 180.0);
	}

	this->coordinates.workingDec = bHold;
	this->coordinates.workingRa = lHold;

	//updateAngDistTemp(angDistCenterDec);
}
void Scan::undoDynamicCosDecTransform(double t_int, double CENTERLATI, double CENTERLONG, std::vector<double> LATITUDES)
{
	std::vector<double> decHold, raHold, timeHold;
	std::vector<bool> edgeFlagHold;
	decHold.resize(scanProperties.size);
	raHold.resize(scanProperties.size);
	timeHold.resize(scanProperties.size);
	edgeFlagHold.resize(scanProperties.size, 0);
	int jMin, jMax;
	for (int j = 0; j < scanProperties.size; j++)
	{
		if (t_int >= 0)
		{
			jMin = j;
			while ((dataProperties.time[jMin] >= dataProperties.time[j] - t_int) && jMin > 0)
			{
				jMin = jMin - 1;
			}
			jMax = jMin + 1;
		}
		else
		{
			jMax = j;
			while ((dataProperties.time[jMax] <= dataProperties.time[j] - t_int) && jMax < scanProperties.size - 1)
			{
				jMax = jMax + 1;
			}
			jMin = jMax - 1;
		}

		decHold[j] = coordinates.workingDec[jMax] + CENTERLATI;
		raHold[j] = (coordinates.workingRa[jMax] / cos(LATITUDES[jMax]*M_PI / 180.0)) + CENTERLONG;
	}

	this->coordinates.workingDec = decHold;
	this->coordinates.workingRa = raHold;

}


//point removal
void Scan::removePoints(std::vector<bool> flags)
{
	std::vector<double> wcHold, tgHold, ddHold, raHold, decHold, prHold, lcHold;
	//this->preRFIRa = workingRa;
	//this->preRFIDec = workingDec;
	//this->preRFIThetaGap = thetaGap;



	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			raHold.push_back(this->coordinates.workingRa[i]);
			decHold.push_back(this->coordinates.workingDec[i]);
			wcHold.push_back(this->fluxes.workingChannel[i]);
			tgHold.push_back(this->thetaGap[i]);
			ddHold.push_back(this->dataProperties.dataDumps[i]);
			lcHold.push_back(this->photometryProperties.SSSCorrelation[i]);

			//prHold.push_back(this->postRFIFlux[i]);

		}
	}

	std::ofstream testFile;

	this->coordinates.workingRa = raHold;
	this->coordinates.workingDec = decHold;
	this->fluxes.workingChannel = wcHold;
	this->thetaGap = tgHold;
	this->dataProperties.dataDumps = ddHold;
	this->photometryProperties.SSSCorrelation = lcHold;
	//this->postRFIFlux = prHold;
	//this->preRFISize = this->size;
	this->scanProperties.size = tgHold.size();


	//updateAngDist();
	//cosDecTransform(fileType, 0.0, centerDec);
}
void Scan::removeRFI(std::vector<double> vec)
{
	std::vector<double> wcHold, tgHold, ddHold, raHold, decHold, prHold, lcHold;
	//this->preRFIRa = workingRa;
	//this->preRFIDec = workingDec;
	//this->preRFIThetaGap = thetaGap;


	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i] == vec[i])
		{
			wcHold.push_back(vec[i]);
			raHold.push_back(this->coordinates.workingRa[i]);
			decHold.push_back(this->coordinates.workingDec[i]);
			tgHold.push_back(this->thetaGap[i]);
			ddHold.push_back(this->dataProperties.dataDumps[i]);
			lcHold.push_back(this->photometryProperties.SSSCorrelation[i]);

			//prHold.push_back(this->postRFIFlux[i]);

		}
	}

	std::ofstream testFile;

	this->coordinates.workingRa = raHold;
	this->coordinates.workingDec = decHold;
	this->fluxes.workingChannel = wcHold;
	this->thetaGap = tgHold;
	this->dataProperties.dataDumps = ddHold;
	this->photometryProperties.SSSCorrelation = lcHold;
	//this->postRFIFlux = prHold;
	//this->preRFISize = this->size;
	this->scanProperties.size = tgHold.size();


	//updateAngDist();
	//cosDecTransform(fileType, 0.0, centerDec);
}



void Scan::setLSSStruct()
{
	this->LSSRa = coordinates.workingRa;
	this->LSSDec = coordinates.workingDec;
	this->LSSDataDumps = dataProperties.dataDumps;
	this->scanProperties.LSSSize = LSSRa.size();
}
void Scan::removeLSSPoints(std::vector<bool> flags)

{
	std::vector<double> wcHold, tgHold, ddHold, raHold, decHold, prHold, lcHold;

	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			raHold.push_back(this->LSSRa[i]);
			decHold.push_back(this->LSSDec[i]);
			wcHold.push_back(this->LSSData[i]);
			tgHold.push_back(this->LSSThetaGap[i]);
			ddHold.push_back(this->LSSDataDumps[i]);
			lcHold.push_back(this->LSSLocalModelCorrelation[i]);
		}
	}

	this->LSSRa = raHold;
	this->LSSDec = decHold;
	this->LSSData = wcHold;
	this->LSSThetaGap = tgHold;
	this->LSSDataDumps = ddHold;
	this->LSSLocalModelCorrelation = lcHold;

	this->scanProperties.LSSSize = tgHold.size();

}

void Scan::removeBG(std::vector<double> background)
{
	this->LSSData = background;
	for (int i = 0; i < scanProperties.size; i++)
	{
		this->fluxes.workingChannel[i] = fluxes.workingChannel[i] - background[i];
	}
}
//void Scan::removeRFI(std::vector<double> rfi)
//{
//	//CHECKS IF RFI IS TURNED OFF
//	double sum = 0;
//	double val = 0;
//	for (int i = 0; i < scanProperties.size; i++)
//	{
//		sum += rfi[i];
//	}
//
//	for (int i = 0; i < scanProperties.size; i++)
//	{
//		val = rfi[i];
//		if (sum != 0.0 && val != 0.0)
//		{
//			this->fluxes.workingChannel[i] = rfi[i];
//		}
//	}
//}
void Scan::removeElevation(std::vector<double> elevationBG)
{
	for (int i = 0; i < scanProperties.size; i++)
	{
		this->LSSData[i] = LSSData[i] - elevationBG[i];
	}
}


//getters
int Scan::getSize()
{
	return this->scanProperties.size;
}
int Scan::getRawSize()
{
	return this->scanProperties.rawSize;
}
int Scan::getLSSSize()
{
	return this->scanProperties.LSSSize;
}
int Scan::getCenter()
{
	return scanProperties.center;
}
int Scan::getSurveyNumber()
{
	return scanProperties.surveyNumber;
}
int Scan::getMaxIndex()
{
	return scanProperties.maxIndex;
}
int Scan::getMinIndex()
{
	return scanProperties.minIndex;
}
int Scan::getMaxIndexRaw()
{
	return scanProperties.maxIndexRaw;
}
int Scan::getMinIndexRaw()
{
	return scanProperties.minIndexRaw;
}
int Scan::getEdgePointFlag(int index)
{
	return this->flags.edgePointFlag[index];
}
int Scan::getTurningPointFlag(int index)
{
	return this->flags.turningPointFlag[index];
}
int Scan::getScanNumberInSurvey()
{
	return scanProperties.scanNumberInSurvey;
}
int Scan::getLocalModelInstances(int index)
{
	return photometryProperties.localModelInstances[index];
}
int Scan::getLocalModelRejections(int index)
{
	return photometryProperties.localModelRejections[index];
}
int Scan::getLocalModelCount(int index)
{
	return this->photometryProperties.localModelCount[index];
}
int Scan::getExtraLocalModelCount(int index)
{
	return this->photometryProperties.extraLocalModelCount[index];
}

double Scan::getElevationScatter()
{
	return this->scanProperties.elevationScatter;
}
double Scan::getDumpSum()
{
	return this->scanProperties.cleanDumpSum;
}
double Scan::getIntraScanGap()
{
	return scanProperties.intraScanGap;
}
double Scan::getInterScanGap()
{
	return scanProperties.interScanGap;
}
double Scan::getRawRa(int index)
{
	return this->coordinates.rawRa[index];
}
double Scan::getRawDec(int index)
{
	return this->coordinates.rawDec[index];
}
double Scan::getLSSRa(int index)
{
	return this->LSSRa[index];
}
double Scan::getLSSDec(int index)
{
	return this->LSSDec[index];
}
double Scan::getScatter()
{
	return this->scanProperties.scatter;
}
double Scan::getAngDist(int index)
{
	return this->dataProperties.angDist[index];
}
double Scan::getFlux(int index)
{
	return this->fluxes.workingChannel[index];
}
double Scan::getDec(int index)
{
	return this->coordinates.workingDec[index];
}
double Scan::getRa(int index)
{
	return this->coordinates.workingRa[index];
}
double Scan::getOrigDec(int index)
{
	return this->coordinates.dec[index];
}
double Scan::getOrigRa(int index)
{
	return this->coordinates.ra[index];
}
double Scan::getTSRa(int index)
{
	return this->coordinates.TSRa[index];
}
double Scan::getTSDec(int index)
{
	return this->coordinates.TSDec[index];
}
double Scan::getTime(int index)
{
	return this->dataProperties.time[index];
}
double Scan::getDataDumps(int index)
{
	return this->dataProperties.dataDumps[index];
}
double Scan::getThetaGap(int index)
{
	return this->thetaGap[index];
}
double Scan::getLSSData(int index)
{
	return this->LSSData[index];
}
double Scan::getLSSThetaWPrime(int index)
{
	return LSSThetaWPrime[index];
}
double Scan::getRCRMinThetaGap(int index)
{
	return minRCRThetaGap[index];
}
double Scan::getGMWeight(int index)
{
	return this->photometryProperties.GMWeights[index];
}
double Scan::getSSSCorrelation(int index)
{
	return this->photometryProperties.SSSCorrelation[index];
}

double Scan::getDrop2dValues(int index)
{
	return this->drop2DValues[index];
}
double Scan::getLSSThetaGap(int index)
{
	return this->LSSThetaGap[index];
}

bool Scan::getScanInRa()
{
	return scanProperties.scanInRa;
}
bool Scan::getRFIFlags(int index)
{
	return this->flags.rfiFlags[index];
}
bool Scan::getCentroidFlag(int index)
{
	return this->flags.centroidFlag[index];
}
int Scan::getDropFlag(int index)
{
	return this->flags.dropFlag[index];
}
bool Scan::getLastScan()
{
	return scanProperties.lastScan;
}
int Scan::getIntHolder(int index)
{
	return this->intHolder[index];
}


std::vector<double> Scan::getFlux()
{
	return this->fluxes.workingChannel;
}
std::vector<double> Scan::getDec()
{
	return this->coordinates.workingDec;
}
std::vector<double> Scan::getRawDec()
{
	return this->coordinates.rawDec;
}
std::vector<double> Scan::getLSSDec()
{
	return LSSDec;
}
std::vector<double> Scan::getRa()
{
	return this->coordinates.workingRa;
}
std::vector<double> Scan::getRawRa()
{
	return this->coordinates.rawRa;
}
std::vector<double> Scan::getLSSRa()
{
	return LSSRa;
}


std::vector<double> Scan::getTSDec()
{
	return this->coordinates.TSDec;
}
std::vector<double> Scan::getTSRa()
{
	return this->coordinates.TSRa;
}
std::vector<double> Scan::getAngDist()
{
	return this->dataProperties.angDist;
}
std::vector<double> Scan::getTime()
{
	return this->dataProperties.time;
}
std::vector<double> Scan::getThetaGap()
{
	return this->thetaGap;
}
std::vector<double> Scan::getDataDumps()
{
	return this->dataProperties.dataDumps;
}
std::vector<double> Scan::getLSSData()
{
	return this->LSSData;
}
std::vector<double> Scan::getDrop2dValues()
{
	return this->drop2DValues;
}
std::vector<double> Scan::getElevation()
{
	return this->dataProperties.elevation;
}



//setters

	//vectors
void Scan::setSize(int newSize)
{
	this->scanProperties.size = newSize;
}
void Scan::setRa(std::vector<double> newRa)
{
	this->coordinates.ra.erase(coordinates.ra.begin(), coordinates.ra.begin() + coordinates.ra.size());
	this->coordinates.ra = newRa;
	this->coordinates.workingRa = coordinates.ra;
}
void Scan::setDec(std::vector<double> newDec)
{
	this->coordinates.dec.erase(coordinates.dec.begin(), coordinates.dec.begin() + coordinates.dec.size());
	this->coordinates.dec = newDec;
	this->coordinates.workingDec = coordinates.dec;

}
void Scan::setTSRa(std::vector<double> TSRaHold)
{
	this->coordinates.TSRa = TSRaHold;
}
void Scan::setTSDec(std::vector<double> TSDecHold)
{
	this->coordinates.TSDec = TSDecHold;
}
void Scan::setFlux(int index, double value)
{
	this->fluxes.workingChannel[index] = value;
}
void Scan::setFlux(std::vector<double> fluxVec)
{
	this->fluxes.workingChannel = fluxVec;
}
void Scan::setLSSData(std::vector<double> lssData)
{
	this->LSSData = lssData;
}
void Scan::setThetaGap(std::vector<double> thetaGapVec)
{
	this->thetaGap = thetaGapVec;
}
void Scan::setThetaGap(int index, double value)
{
	this->thetaGap[index] = value;
}
void Scan::setRFIFlags(std::vector<bool> RFIFlag)
{
	this->flags.rfiFlags = RFIFlag;
}
void Scan::setBG(std::vector<double> background)
{
	this->LSSData = background;
}
void Scan::setTurningPointFlag(std::vector<bool> turnPointVec)
{
	this->flags.turningPointFlag = turnPointVec;
}
void Scan::setEdgePointFlag(std::vector<bool> edgePointVec)
{
	this->flags.edgePointFlag = edgePointVec;
}
void Scan::setEdgePointFlag(int index, bool flag)
{
	this->flags.edgePointFlag[index] = flag;
}
void Scan::setCentroidFlag(std::vector<bool> centroidFlags)
{
	this->flags.centroidFlag = centroidFlags;
}
void Scan::setCentroidFlag(int index, bool value)
{
	this->flags.centroidFlag[index] = value;
}
void Scan::setLSSData(int index, double value)
{
	this->LSSData[index] = value;
}
void Scan::setLSSThetaWPrime(int index, double value)
{
	this->LSSThetaWPrime[index] = value;
}
void Scan::setDrop2dValues(int index, double value)
{
	this->drop2DValues[index] = value;
}
void Scan::setLSSThetaGap(std::vector<double> LSSThetaGap)
{
	this->LSSThetaGap = LSSThetaGap;
}
void Scan::setLSSThetaGap(int index, double value)
{
	this->LSSThetaGap[index] = value;
}
void Scan::setDropFlag(int index, int value)
{
	this->flags.dropFlag[index] = value;
}
void Scan::setIntHolder(int index, int value)
{
	this->intHolder[index] = value;
}

	//single values
void Scan::setMaxIndex(int value)
{
	this->scanProperties.maxIndex = value;
}
void Scan::setMinIndex(int value)
{
	this->scanProperties.minIndex = value;
}
void Scan::setMinRCRThetaGap(int index, double value)
{
	this->minRCRThetaGap[index] = value;
}
void Scan::setSurveyNumber(int num)
{
	this->scanProperties.surveyNumber = num;
}
void Scan::setScanNumberInSurvey(int num)
{
	this->scanProperties.scanNumberInSurvey = num;
}
void Scan::setLocalModelInstances(int index, int value)
{
	this->photometryProperties.localModelInstances[index] = value;
}
void Scan::setLocalModelRejections(int index, int value)
{
	this->photometryProperties.localModelRejections[index] = value;
}
void Scan::setGMWeight(int index, double value)
{
	this->photometryProperties.GMWeights[index] = value;
}
void Scan::setSSSCorrelation(int index, double value)
{
	this->photometryProperties.SSSCorrelation[index] = value;
}
void Scan::setLocalModelCount(int index, double value)
{
	this->photometryProperties.localModelCount[index] = value;
}
void Scan::setExtraLocalModelCount(int index, double value)
{
	this->photometryProperties.extraLocalModelCount[index] = value;
}
void Scan::setLastScan(bool val)
{
	this->scanProperties.lastScan = val;
}
void Scan::setScanInRa(bool value)
{
	this->scanProperties.scanInRa = value;
}
void Scan::setCenter(int index)
{
	this->scanProperties.center = index;
}
void Scan::setScatter(double scatter)
{
	this->scanProperties.scatter = scatter;
}
void Scan::setElevationScatter(double scatter)
{
	this->scanProperties.elevationScatter = scatter;
}
void Scan::setIntraScanGap(double value)
{
	this->scanProperties.intraScanGap = value;
}
void Scan::setInterScanGap(double value)
{
	this->scanProperties.interScanGap = value;
}

Scan::~Scan()
{
	//dtor
}