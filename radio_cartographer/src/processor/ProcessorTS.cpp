#include "processor\ProcessorTS.h"
#include "utils\Tools.h"
#include <algorithm>
#include "utils\RCR.h"


ProcessorTS::ProcessorTS()
{
}
ProcessorTS::ProcessorTS(std::vector<Scan> &scans, MapTypes mapType, double psfFWHM, double medianDiffAlongSweeps, bool scansInRa)
{
	this->scans = scans;
	this->mapType = mapType;
	this->psfFWHM = psfFWHM;
	this->medianDiffAlongSweeps = medianDiffAlongSweeps;
	this->scansInRa = scansInRa;
}


std::vector<bool> ProcessorTS::shiftRejection(std::vector<double> shifts)
{
	double minAngle1, maxAngle1, minAngle2, maxAngle2, angleHold, minAngle, maxAngle, angleStep, maxIFFT, hold;
	double deltaI, deltaJ;
	double probHold, prob;
	double probHold1, probHold2;
	double R;
	int j, M, N;
	std::vector<double> dubFiller, angleVec1, angleVec2;
	std::vector<bool> flagsHold;
	flagsHold.resize(shifts.size(), 1);
	N = shifts.size();

	//FIND MINIMUM OF MAXIUMA AND MAXIMUM OF MINIMA FOR ALL SCANS TO CALCULATE R
	for (int i = 0; i < scans.size() - 1; i++)
	{
		minAngle1 = 999999;
		maxAngle1 = -999999;
		minAngle2 = 999999;
		maxAngle2 = -999999;
		if (mapType == DAISY)
		{
			angleVec1 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg, scans[i].getRa(), scans[i].getDec());
			angleVec2 = Tools::daisyAngleBuilder(i + 1, scans[i + 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i + 1].getRa(), scans[i + 1].getDec());
		}
		else
		{
			if (scansInRa)
			{
				angleVec1 = scans[i].getRa();
				angleVec2 = scans[i + 1].getRa();
			}
			else
			{
				angleVec1 = scans[i].getDec();
				angleVec2 = scans[i + 1].getDec();
			}
		}
		for (int j = 0; j < angleVec1.size(); j++)
		{
			angleHold = angleVec1[j];
			if (angleHold < minAngle1)
			{
				minAngle1 = angleHold;
			}
			if (angleHold > maxAngle1)
			{
				maxAngle1 = angleHold;
			}
		}
		for (int j = 0; j < angleVec2.size(); j++)
		{
			angleHold = angleVec2[j];
			if (angleHold < minAngle2)
			{
				minAngle2 = angleHold;
			}
			if (angleHold > maxAngle2)
			{
				maxAngle2 = angleHold;
			}
		}
		minAngle = std::max(minAngle1, minAngle2);
		maxAngle = std::min(maxAngle1, maxAngle2);
	}

	R = maxAngle - minAngle;

	//DETERMINE LIKELINESS OF THE ADJACENT SCANS TO HAVE SIMILAR SHIFTS
	for (int i = 0; i < shifts.size(); i++)
	{

		deltaI = shifts[i];

		//IF ON EDGE SCAN, ONLY USE ONE ADJACENT SCAN'S SHIFT

		//IF BOTH ARE ACCESSIBLE USE BOTH

		if (i == 0)
		{
			j = i + 1;
			M = 1;
		}
		else if (i == shifts.size() - 1)
		{
			j = i - 1;
			M = 1;
		}
		else
		{
			M = 2;
		}

		//CALCULATE LIKELINESS

		if (M == 1)
		{
			deltaJ = shifts[j];
			probHold = (std::abs(deltaI - deltaJ) / R) * (2 - (std::abs(deltaI - deltaJ) / R));
			prob = pow(2, (M - 1))*probHold;
		}
		else
		{
			j = i + 1;
			deltaJ = shifts[j];
			probHold1 = (std::abs(deltaI - deltaJ) / R) * (2 - (std::abs(deltaI - deltaJ) / R));

			j = i - 1;
			deltaJ = shifts[j];
			probHold2 = (std::abs(deltaI - deltaJ) / R) * (2 - (std::abs(deltaI - deltaJ) / R));

			prob = pow(2, (M - 1))*(probHold1*probHold2);
		}

		//REJECT POINTS BASED ON LIKELINESS
		if (prob > (1.0 / (2.0 * N)))
		{
			flagsHold[i] = 0;
		}
	}

	return flagsHold;
}
double ProcessorTS::calculateShift(std::vector<std::vector<double> > &maxIfftInfo)
{
	RCR rcr = RCR(LS_MODE_DL);
	bool stop = false, first = true, alreadySorted = false;
	double muBelow, muAbove, stDevAbove, stDevBelow, belowAboveRatio, maxBelowAboveRatio = -999999;
	double maxBelowAboveIndex, shiftWeightSum = 0.0, weightSum = 0.0;
	double result;
	int counter2 = 0;
	std::vector<bool> checks, shiftRejectionFlags, flagHold;
	std::vector<double> swap, shiftsAbove, shiftsBelow, weightTemps, shifts, weights, tempArray, tempWeights, tempArray2, tempWeights2, weightsBelow, weightsAbove, ratioAbove, ratioBelow;
	std::vector<double> shiftsFinal;

	//CALCULATE THE CORRELATION BETWEEN SCANS
	for (int j = 0; j < maxIfftInfo.size(); j++)
	{
		maxIfftInfo[j][0] = -1.0*pow(-1.0, j) * maxIfftInfo[j][0];
		shifts.push_back(maxIfftInfo[j][0]);
	}

	//REJECT NON-SENSIBLE CORRELATIONS DUE TO RFI
	shiftRejectionFlags.resize(shifts.size(), 1);
	shiftRejectionFlags = shiftRejection(shifts);

	for (int i = 0; i < shifts.size(); i++)
	{
		if (shiftRejectionFlags[i] == 1)
		{
			shiftsFinal.push_back(shifts[i]);
		}
	}


	//CALCULATE AVERAGE SHIFT
	if (shiftsFinal.size() == 0)
	{
		result = 0.0;
	}
	else
	{
		rcr.performBulkRejection(shiftsFinal);
		result = rcr.result.mu;

		flagHold = rcr.result.flags;

		for (int i = 0; i < shiftsFinal.size(); i++)
		{
			if (flagHold[i] == 1)
			{
				counter2++;
			}
		}
		if (counter2 == 0)
		{
			result = 0.0;
		}
	}

	return result;
}
std::vector<std::vector<double>> ProcessorTS::getScanToScanShifts(double sampling)
{
	int requiredPower, maxIFFTBin, h, powerOffset;
	double minAngle1, maxAngle1, minAngle2, maxAngle2, angleHold, minAngle, maxAngle, angleStep, maxIFFT, hold;
	std::vector<double> dubFiller, angleVec1, angleVec2;
	std::vector<std::vector<double> > interpolated, ifftResult, maxIfftInfo;
	for (int i = 0; i < scans.size() - 1; i++)
	{
		minAngle1 = 999999;
		maxAngle1 = -999999;
		minAngle2 = 999999;
		maxAngle2 = -999999;
		if (mapType == DAISY)
		{
			angleVec1 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg, scans[i].getRa(), scans[i].getDec());
			angleVec2 = Tools::daisyAngleBuilder(i + 1, scans[i + 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i + 1].getRa(), scans[i + 1].getDec());
		}
		else
		{
			if (scansInRa)
			{
				angleVec1 = scans[i].getRa();
				angleVec2 = scans[i + 1].getRa();
			}
			else
			{
				angleVec1 = scans[i].getDec();
				angleVec2 = scans[i + 1].getDec();
			}
		}
		for (int j = 0; j < angleVec1.size(); j++)
		{
			angleHold = angleVec1[j];
			if (angleHold < minAngle1)
			{
				minAngle1 = angleHold;
			}
			if (angleHold > maxAngle1)
			{
				maxAngle1 = angleHold;
			}
		}
		for (int j = 0; j < angleVec2.size(); j++)
		{
			angleHold = angleVec2[j];
			if (angleHold < minAngle2)
			{
				minAngle2 = angleHold;
			}
			if (angleHold > maxAngle2)
			{
				maxAngle2 = angleHold;
			}
		}
		minAngle = std::max(minAngle1, minAngle2);
		maxAngle = std::min(maxAngle1, maxAngle2);
		angleStep = std::min((maxAngle1 - minAngle1) / (double(scans[i].getSize())), (maxAngle2 - minAngle2) / (double(scans[i + 1].getSize())));//
		requiredPower = 1;
		while (pow(2, requiredPower) < (maxAngle - minAngle) / angleStep)
		{
			requiredPower++;
		}
		powerOffset = floor(log(15.0 * sampling / psfFWHM) / log(2.0)) + 4;
		requiredPower += powerOffset;
		angleStep = (maxAngle - minAngle) / (pow(2, requiredPower));
		interpolated.resize(3, dubFiller);
		if (angleVec1[angleVec1.size() - 1] > angleVec1[0])
		{
			angleHold = std::max(minAngle1, minAngle2);
			h = 1;
			for (int k = 0; k < pow(2, requiredPower); k++)
			{
				interpolated[0].push_back(angleHold);
				while (h < angleVec1.size() - 1 && angleVec1[h] < angleHold)
				{
					h++;
				}
				if (angleVec1[h] == angleVec1[h - 1])
				{
					interpolated[1].push_back(scans[i].getFlux(h));
				}
				else
				{
					interpolated[1].push_back(scans[i].getFlux(h - 1) + (angleHold - angleVec1[h - 1])*(scans[i].getFlux(h) - scans[i].getFlux(h - 1)) / (angleVec1[h] - angleVec1[h - 1]));
				}
				angleHold += angleStep;
			}
			angleHold = std::max(minAngle1, minAngle2);
			h = scans[i + 1].getSize() - 2;
			for (int k = 0; k < pow(2, requiredPower); k++)
			{
				while (h > 0 && angleVec2[h] < angleHold)
				{
					h--;
				}
				if (angleVec2[h] == angleVec2[h + 1])
				{
					interpolated[2].push_back(scans[i + 1].getFlux(h));
				}
				else
				{
					interpolated[2].push_back(scans[i + 1].getFlux(h + 1) + (angleHold - angleVec2[h + 1])*(scans[i + 1].getFlux(h) - scans[i + 1].getFlux(h + 1)) / (angleVec2[h] - angleVec2[h + 1]));
				}
				angleHold += angleStep;
			}
		}
		else
		{
			angleHold = std::max(minAngle1, minAngle2);
			h = scans[i].getSize() - 2;
			for (int k = 0; k < pow(2, requiredPower); k++)
			{
				interpolated[0].push_back(angleHold);
				while (h > 0 && angleVec1[h] < angleHold)
				{
					h--;
				}
				if (angleVec1[h] == angleVec1[h + 1])
				{
					interpolated[1].push_back(scans[i].getFlux(h));
				}
				else
				{
					interpolated[1].push_back(scans[i].getFlux(h + 1) + (angleHold - angleVec1[h + 1])*(scans[i].getFlux(h) - scans[i].getFlux(h + 1)) / (angleVec1[h] - angleVec1[h + 1]));
				}
				angleHold += angleStep;
			}
			angleHold = std::max(minAngle1, minAngle2);
			h = 1;
			for (int k = 0; k < pow(2, requiredPower); k++)
			{
				while (h < angleVec2.size() - 1 && angleVec2[h] < angleHold)
				{
					h++;
				}
				if (angleVec2[h] == angleVec2[h - 1])
				{
					interpolated[2].push_back(scans[i].getFlux(h));
				}
				else
				{
					interpolated[2].push_back(scans[i + 1].getFlux(h - 1) + (angleHold - angleVec2[h - 1])*(scans[i + 1].getFlux(h) - scans[i + 1].getFlux(h - 1)) / (angleVec2[h] - angleVec2[h - 1]));
				}
				angleHold += angleStep;
			}
		}


		ifftResult = Tools::crossCorrelate(interpolated[1], interpolated[2]);
		maxIFFT = -999999;
		for (int j = 0; j < ifftResult.size(); j++)
		{
			hold = sqrt(ifftResult[j][0] * ifftResult[j][0] + ifftResult[j][1] * ifftResult[j][1]);
			if (hold > maxIFFT)
			{
				maxIFFT = hold;
				maxIFFTBin = j;
			}
		}
		maxIfftInfo.push_back(dubFiller);
		if (maxIFFTBin > pow(2, requiredPower - 1))
		{
			maxIFFTBin -= pow(2, requiredPower);
		}
		maxIfftInfo[i].push_back(maxIFFTBin * angleStep);
		maxIfftInfo[i].push_back(maxIFFT / pow(4.0, requiredPower));
		maxIfftInfo[i].push_back(requiredPower);
		interpolated.clear();
		ifftResult.clear();
		angleVec1.clear();
		angleVec2.clear();
	}

	return maxIfftInfo;
}
double ProcessorTS::find_tInt()
{
	double t_int = 0;
	RCR rcr = RCR(LS_MODE_DL);
	int centerIndex;
	bool stop = false;
	double negation = -1.0;
	double shift, median = 0, sampling;
	std::vector<int> checks;
	std::vector<double> angleHold, angDists, speeds, tempArray, tempWeights, angDistHold, timeHold;
	std::vector<double> centers, speedTemp;
	std::vector<std::vector<double> > maxIfftInfo;
	if (mapType == DAISY)
	{
		for (int i = 0; i < scans.size(); i++)
		{
			centerIndex = scans[i].getCenter();
			speeds.push_back((scans[i].getAngDist(centerIndex + 1) - scans[i].getAngDist(centerIndex - 1)) / (scans[i].getTime(centerIndex + 1) - scans[i].getTime(centerIndex - 1)));
			angDists.push_back((scans[i].getAngDist(centerIndex + 1) - scans[i].getAngDist(centerIndex - 1)) / 2.0);
		}
		rcr.performBulkRejection(angDists);
		sampling = rcr.result.mu;
	}
	else
	{
		for (int i = 0; i < scans.size(); i++)
		{
			if (scansInRa)
			{
				angDistHold = scans[i].getRa();// getAngDist();// STOP CHECK
			}
			else
			{
				angDistHold = scans[i].getDec();
			}

			timeHold = scans[i].getTime();
			for (int j = 1; j < scans[i].getSize(); j++)
			{
				speedTemp.push_back((angDistHold[j] - angDistHold[j - 1]) / (timeHold[j] - timeHold[j - 1]));
			}
			rcr.performBulkRejection(speedTemp);
			speeds.push_back(rcr.result.mu);
			speedTemp.clear();
		}
		sampling = medianDiffAlongSweeps;
	}

	maxIfftInfo = getScanToScanShifts(sampling);
	shift = calculateShift(maxIfftInfo) / 2.0;

	if (mapType == DAISY)
	{
		rcr.performRejection(speeds);
	}
	else
	{
		rcr.performBulkRejection(speeds);
	}


	t_int = shift / rcr.result.mu;


	std::ofstream outputFileAux;
	outputFileAux.open("TS.txt", std::ios_base::app);
	outputFileAux << t_int << "\n";
	outputFileAux.close();

	return t_int;
}

ProcessorTS::~ProcessorTS()
{
}