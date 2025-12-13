#include "ProcessorRFI.h"
#include <iostream>
#include "Tools.h"
#include "math.h"
#include "RCR.h"

static double toRad = M_PI / 180.0;
static double toDeg = 180.0 / M_PI;
std::vector<std::vector<double>> spawnNewRCRThread2(std::vector<double> w, std::vector<double> y)
{
	std::vector<std::vector<double>> rcrOutput;
	rcrOutput.resize(2);
	RCR rcr = RCR(LS_MODE_DL);
	rcr.performBulkRejection(w, y);

	std::vector<double> flags;
	flags.resize(w.size());
	for (int i = 0; i < w.size(); i++)
	{
		flags[i] = (double)rcr.result.flags[i];
	}

	rcrOutput[0].push_back(rcr.result.mu);
	rcrOutput[1] = flags;
	return rcrOutput;
}

ProcessorRFI::ProcessorRFI(std::vector<Scan> &scans, RFIParameters & rfiParams, std::vector<std::vector<double>> & scatter2dSSS, std::vector<std::vector<std::vector<int>>> & classificationsSSS)
{
	this->scans = scans;
	this->scatter2dSSS = scatter2dSSS;
	this->classificationsSSS = classificationsSSS;

	this->psfFWHM = rfiParams.psfFWHM;
	this->rfiScaleBW = rfiParams.rfiScaleBW;
	this->correlatedWeightMap = rfiParams.correlatedWeightMap;
	this->standardGap = rfiParams.standardGap;
	this->partSetProcSSS = rfiParams.partSetProcSSS;
	this->rfiScaleDeg = rfiScaleBW * psfFWHM;

	double decHold, raHold;
	std::vector<bool> boolFiller;
	std::vector<double> dubFiller;
	std::vector<std::vector<bool>> flags;
	std::vector<std::vector<double> > dubFillerFiller;
	
	flags.resize(scans.size(), boolFiller);
	rfiSubtracted.resize(scans.size(), dubFiller);
	rfiValues.resize(scans.size(), dubFillerFiller);
	rfiWeights.resize(scans.size(), dubFillerFiller);
	rfiCounts.resize(scans.size(), dubFillerFiller);
	rfiDists.resize(scans.size(), dubFillerFiller);
	rfiNs.resize(scans.size(), dubFillerFiller);

	for (int i = 0; i < scans.size(); i++)
	{
		rfiValues[i].resize(scans[i].getSize(), dubFiller);
		rfiWeights[i].resize(scans[i].getSize(), dubFiller);
		rfiCounts[i].resize(scans[i].getSize(), dubFiller);
		rfiDists[i].resize(scans[i].getSize(), dubFiller);
		rfiNs[i].resize(scans[i].getSize(), dubFiller);

		rfiSubtracted[i].resize(scans[i].getSize(), 0.0);
		flags[i].resize(scans[i].getSize(), true);
	}

	if (rfiScaleBW != 0.0)
	{
		rfiBuildLocalModels();
		autoCentroid(rfiValues);


		if ((centroidLocations.size() != 0) && ((rfiScaleBW < 1.0) && (rfiScaleBW > pow((standardGap / psfFWHM), 0.5))))
		{
			for (int i = 0; i < centroidLocations.size() / 2; i++)
			{
				decHold = centroidLocations[2 * i];
				raHold = centroidLocations[2 * i + 1];

				//std::cout << "Extra Local Models: \t" << decHold << "\t" << raHold << "\n";
				rfiBuildExtraModels(decHold, raHold);
			}
		}
		rfiSubtracted = rfiBuildGlobalModels();
		setGlobalModelValues(scans, this->scans);
	}
}
ProcessorRFI::~ProcessorRFI()
{
}

void ProcessorRFI::autoCentroid(std::vector<std::vector<std::vector<double> > > rfiValues)
{
	bool firstIteration = true;
	std::vector<bool> sourceFlag;
	std::vector<double> sourceFlux;
	std::vector<int> sourceScan, sourceIndex;

	//FLAGS POINTS WITH FLUX ABOVE 75x LARGER THAN SCATTER AND THAT ARE NOT RFI
	//USED TO BE 100x, NOW 75x. (DAN'S DECISION).
	for (int i = 0; i < scans.size(); i++)
	{
		//std::cout << "scatter:\t" << scans[i].getScatter() << "\n";
		sourceFlag.resize(scans[i].getSize(), true);
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			if ((scans[i].getFlux(j) < (75 * scans[i].getScatter())) || (rfiValues[i][j].size() == 0))
			{
				sourceFlag[j] = false;
			}
			else
			{
				sourceFlag[j] = true;
				sourceFlux.push_back(scans[i].getFlux(j));
				sourceScan.push_back(i);
				sourceIndex.push_back(j);
			}
		}
		scans[i].setCentroidFlag(sourceFlag);
	}

	if (sourceFlux.size() == 0)
	{
		return;
	}


	//FINDS THE LARGEST FLUX
	Tools::sortAll(0, sourceFlux.size() - 1, sourceFlux, sourceScan, sourceIndex);


	//CALCULATES THE WEIGHTED CENTER OF MASS OF FLAGGED FLUXES WITHIN A ONE BEAM WIDTH RADIUS
	bool newSource = false;
	int iHold, jHold, iTemp, jTemp;
	int counter, chs, pointCandidates;
	double wDec, wRa;
	double comDec, comRa;
	double maxFlux, totalFlux;
	double distance;
	double centroidDistance;
	double potentialRa, potentialDec;
	std::vector<double> coordinates, inRange;
	std::vector<double> centroidDistances;
	std::vector<int> possibles;


	while (sourceFlux.size() != 0)
	{
		newSource = false;
		maxFlux = sourceFlux[sourceFlux.size() - 1];
		iHold = sourceScan[sourceScan.size() - 1];
		jHold = sourceIndex[sourceIndex.size() - 1];

		scans[iHold].setCentroidFlag(jHold, false);

		distance = Tools::getGCDistance(scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), partSetProcSSS.centerDecDeg)*toDeg;

		inRange.push_back(scans[iHold].getFlux(jHold));
		inRange.push_back(scans[iHold].getDec(jHold));
		inRange.push_back(scans[iHold].getRa(jHold));
		inRange.push_back(cos((M_PI*distance) / (psfFWHM*2.0)));

		counter = 1;

		pointCandidates = findPossSSS(scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), psfFWHM, possibles);

		for (int i = 0; i < possibles.size() / 2; i++)
		{
			iTemp = possibles[2 * i];
			jTemp = possibles[2 * i + 1];

			distance = Tools::getGCDistance(scans[iTemp].getDec(jTemp), scans[iTemp].getRa(jTemp), scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), partSetProcSSS.centerDecDeg)*toDeg;

			if (distance < psfFWHM)
			{
				inRange.push_back(scans[iTemp].getFlux(jTemp));
				inRange.push_back(scans[iTemp].getDec(jTemp));
				inRange.push_back(scans[iTemp].getRa(jTemp));
				inRange.push_back(cos((M_PI*distance) / (psfFWHM*2.0)));

				counter += 1;

				scans[iTemp].setCentroidFlag(jTemp, false);
			}
		}

		coordinates = determineCenters(counter, inRange, scans[iHold].getThetaGap(jHold), scans[iHold].getDec(jHold), scans[iHold].getRa(jHold));

		if (firstIteration)
		{
			if (counter >= 6)
			{
				centroidLocations.push_back(scans[iHold].getDec(jHold) + coordinates[0]);
				centroidLocations.push_back(scans[iHold].getRa(jHold) + coordinates[1]);
				firstIteration = false;
			}
		}
		else // checks that the same centroid isn't being detected multiple times.
		{
			potentialDec = scans[iHold].getDec(jHold) + coordinates[0];
			potentialRa = scans[iHold].getRa(jHold) + coordinates[1];

			for (int i = 0; i < (centroidLocations.size() / 2); i++)
			{
				centroidDistance = Tools::getGCDistance(centroidLocations[2 * i], centroidLocations[2 * i + 1], potentialDec, potentialRa, partSetProcSSS.centerDecDeg)*toDeg;
				if (centroidDistance > psfFWHM  && counter >= 6)
				{
					newSource = true;
				}
				else
				{
					newSource = false;
					break;
				}
			}

			if (newSource)
			{
				centroidLocations.push_back(scans[iHold].getDec(jHold) + coordinates[0]);
				centroidLocations.push_back(scans[iHold].getRa(jHold) + coordinates[1]);
			}

		}

		inRange.clear();
		sourceFlux.clear();
		sourceScan.clear();
		sourceIndex.clear();

		for (int i = 0; i < scans.size(); i++)
		{
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				if (scans[i].getCentroidFlag(j) == true)
				{
					sourceFlux.push_back(scans[i].getFlux(j));
					sourceScan.push_back(i);
					sourceIndex.push_back(j);
				}
			}
		}
		if (sourceFlux.size() == 0)
		{
			return;
		}
		else
		{
			Tools::sortAll(0, sourceFlux.size() - 1, sourceFlux, sourceScan, sourceIndex);
		}
	}
}
std::vector<double> ProcessorRFI::determineCenters(int pointCount, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
{
	int columnCount = 6, pivotIndex;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower;
	double w = 0, xw = 0, yw = 0, xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, finalAnswer;
	double ww = 0, sumWSquared = 0;

	double weightingFunction = (1.0 / 3.0);
	weightingFunctionTemp = (1.0 / 3.0);//Tools::max(weightingFunction, Tools::min((4.0 / 3.0)*largestGap*toDeg / psfFWHM, 1.0));//maybe weightTempNeeded
	double alpha = (log(0.5)) / (log(cos(M_PI*weightingFunctionTemp / 4.0)));

	std::vector<double> A;
	std::vector<double> b;
	std::vector<double> coef;
	std::vector<double> coordinates;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	for (int i = 0; i < pointCount; i++)
	{

		//CHECK THAT THESE MATCH UP
		zHold = inRange[4 * i];
		xHold = (inRange[4 * i + 1] - i_0);
		yHold = (inRange[4 * i + 2] - j_0);

		wHold = pow(inRange[4 * i + 3], alpha);
		ww += std::abs(pow(wHold, 2));
		w += wHold;
		xw += xHold * wHold;
		yw += yHold * wHold;
		xxw += xHold * xHold*wHold;
		yyw += yHold * yHold*wHold;
		xyw += xHold * yHold*wHold;
		xxxw += xHold * xHold*xHold*wHold;
		xyyw += xHold * yHold*yHold*wHold;
		xxyw += xHold * xHold*yHold*wHold;
		yyyw += yHold * yHold*yHold*wHold;
		xxxxw += xHold * xHold*xHold*xHold*wHold;
		xxyyw += xHold * xHold*yHold*yHold*wHold;
		xxxyw += xHold * yHold*xHold*xHold*wHold;
		yyyyw += yHold * yHold*yHold*yHold*wHold;
		xyyyw += xHold * yHold*yHold*yHold*wHold;
		zw += zHold * wHold;
		xzw += xHold * zHold*wHold;
		yzw += yHold * zHold*wHold;
		xxzw += xHold * xHold*zHold*wHold;
		yyzw += yHold * yHold*zHold*wHold;
		xyzw += xHold * yHold*zHold*wHold;
	}

	sumWSquared = std::abs(pow(w, 2));

	A[0] = w; A[1] = xw; A[2] = yw; A[3] = xxw; A[4] = yyw; A[5] = xyw;
	A[6] = xw;  A[7] = xxw; A[8] = xyw; A[9] = xxxw; A[10] = xyyw; A[11] = xxyw;
	A[12] = yw; A[13] = xyw; A[14] = yyw; A[15] = xxyw; A[16] = yyyw; A[17] = xyyw;
	A[18] = xxw; A[19] = xxxw; A[20] = xxyw; A[21] = xxxxw; A[22] = xxyyw; A[23] = xxxyw;
	A[24] = yyw; A[25] = xyyw; A[26] = yyyw; A[27] = xxyyw; A[28] = yyyyw; A[29] = xyyyw;
	A[30] = xyw; A[31] = xxyw; A[32] = xyyw; A[33] = xxxyw; A[34] = xyyyw; A[35] = xxyyw;
	b[0] = zw; b[1] = xzw; b[2] = yzw; b[3] = xxzw; b[4] = yyzw; b[5] = xyzw;

	coef = Tools::matrixSolver(columnCount, A, b);

	double x = ((coef[2] * coef[5]) - (2 * coef[1] * coef[4])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));
	double y = ((coef[1] * coef[5]) - (2 * coef[2] * coef[3])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));

	coordinates.push_back(x);
	coordinates.push_back(y);

	return coordinates;
}

void ProcessorRFI::setGlobalModelValues(std::vector<Scan> &globalScans, std::vector<Scan> localScans)
{
	for (int i = 0; i < localScans.size(); i++)
	{
		for (int j = 0; j < localScans[i].getSize(); j++)
		{
			globalScans[i].setSSSCorrelation(j, localScans[i].getSSSCorrelation(j));
			globalScans[i].setGMWeight(j, localScans[i].getGMWeight(j));
		}
	}
}

int ProcessorRFI::findPossSSS(double dec, double ra, double limit, std::vector<int> &toRet)
{
	std::vector<int> intFiller;
	toRet.resize(100, 0);
	double distance;
	int pointCounter = 0;
	int decCounter = (int)((dec - partSetProcSSS.minDec) / partSetProcSSS.subDecInc);
	int RACounter = (int)((ra - partSetProcSSS.minRa) / partSetProcSSS.subRaInc);
	int iHold, jHold;


	if (decCounter >= classificationsSSS.size())
	{
		decCounter = classificationsSSS.size() - 1;
	}
	if (RACounter >= classificationsSSS[decCounter].size())
	{
		RACounter = classificationsSSS[decCounter].size() - 1;
	}
	if (decCounter < 0)
	{
		decCounter = 0;
	}
	if (RACounter < 0)
	{
		RACounter = 0;
	}
	for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProcSSS.subDecRes); i++)
	{
		for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProcSSS.subRaRes); j++)
		{
			for (int k = 0; k < classificationsSSS[i][j].size(); k += 2)
			{
				iHold = classificationsSSS[i][j][k];
				jHold = classificationsSSS[i][j][k + 1];

				distance = Tools::getGCDistance(dec, ra, scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), partSetProcSSS.centerDecDeg)*toDeg;

				if (distance < limit)
				{
					if (pointCounter < 50)
					{
						toRet[2 * pointCounter] = iHold;
						toRet[2 * pointCounter + 1] = jHold;
					}
					else
					{
						toRet.push_back(iHold);
						toRet.push_back(jHold);
					}
					pointCounter++;
				}
			}
		}
	}
	return pointCounter;
}
std::vector<double> ProcessorRFI::rfiFit(std::vector<int> &checks, std::vector<double> &inRange)
{
	double sigmaTerm = 0, fluxShapeTerm = 0, shapeSquaredTerm = 0, shapeTerm = 0, fluxTerm = 0;
	double f, g, sigma, sigmaSquare, flux, shape, dump;
	std::vector<double> solution;
	solution.resize(2);
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		if (checks[i])
		{
			flux = inRange[6 * i];
			shape = inRange[6 * i + 1];
			dump = inRange[6 * i + 2];
			sigma = inRange[6 * i + 3]; // sigmaValues[i];
			sigmaSquare = sigma * sigma;

			sigmaTerm += dump / sigmaSquare;
			fluxShapeTerm += (dump * flux * shape) / sigmaSquare;
			shapeSquaredTerm += (dump * shape * shape) / sigmaSquare;
			shapeTerm += dump * shape / sigmaSquare;
			fluxTerm += dump * flux / sigmaSquare;
		}
	}
	f = (2 * sigmaTerm * fluxShapeTerm - fluxTerm * shapeTerm) / (2 * sigmaTerm * shapeSquaredTerm - shapeTerm * shapeTerm);
	g = (fluxTerm - f * shapeTerm) / (2 * sigmaTerm);
	solution[0] = f;
	solution[1] = g;
	return solution;
}
void ProcessorRFI::rfiRemovePoints(double avgSigma, std::vector<int> &checks, std::vector<double> &inRange)
{
	double largestDiff, largest, splusSigmaDiff = 999999, splus = 999999;
	double f, g;
	double shape, weight, modelHold;
	int counter, counterplus, largeIndexI, lastChecked;
	int holdI, holdJ;
	std::vector<double> rfiFitSolution;

	while (splus > avgSigma)
	{
		largestDiff = 0;
		splus = 0.0;
		largest = 0.0;
		counter = 0;
		counterplus = 0;
		largeIndexI = -1;
		rfiFitSolution = rfiFit(checks, inRange);
		f = rfiFitSolution[0];
		g = rfiFitSolution[1];

		for (int i = 0; i < inRange.size() / 6; i++)
		{
			if (checks[i])
			{
				shape = inRange[6 * i];
				weight = inRange[6 * i + 1];
				modelHold = f * weight + g;
				counter++;
				if (f > 0)
				{
					if (shape > modelHold)
					{
						splus = splus + pow(shape - modelHold, 2);
						counterplus++;
					}
					if ((shape - modelHold) > largest)
					{
						largest = shape - modelHold;
						largeIndexI = i;
					}
				}
				else
				{
					splus = splus + pow(shape - modelHold, 2);
					counterplus++;
					if (std::abs(shape - modelHold) > largest)
					{
						largest = shape - modelHold;
						largeIndexI = i;
					}
				}
			}
		}

		if (counter > 1)
		{
			splus = sqrt((splus / (counterplus)));

			if (largeIndexI > -1)
			{
				holdI = inRange[6 * largeIndexI + 4];
				holdJ = inRange[6 * largeIndexI + 5];
			}
			if ((splus > avgSigma) && (largeIndexI > -1))
			{
				checks[largeIndexI] = false;
				lastChecked = largeIndexI;
				splusSigmaDiff = std::abs(avgSigma - splus);
			}
			else
			{
				if (splusSigmaDiff < std::abs(avgSigma - splus))
				{
					checks[lastChecked] = false;
				}
				splus = 0;
			}
		}
		else
		{
			splus = 0;
		}
	}
}
void ProcessorRFI::rfiBuildLocalModels()
{
	std::vector<std::future<std::vector<double>>> futureFiller;
	double decHold, raHold, fluxVal, correlatedScatter;
	for (int i = 0; i < scans.size(); i++)
	{
		futureFiller.resize(scans[i].getSize());
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			correlatedScatter = scatter2dSSS[scans[i].getSurveyNumber()][scans[i].getScanNumberInSurvey()];
			fluxVal = scans[i].getFlux(j);
			decHold = scans[i].getDec(j);
			raHold = scans[i].getRa(j);

			futureFiller[j] = std::async(std::launch::async, &ProcessorRFI::rfiRemovalMulti, this, decHold, raHold, fluxVal, correlatedScatter);
		}

		rfiSetFinalValues(futureFiller);
	}
}
void ProcessorRFI::rfiBuildExtraModels(double decPoint, double raPoint)
{
	int extraSteps, counter = 0;
	double deltaX, deltaY, simpleWidth, distance;
	double decHold, raHold, fluxVal, correlatedScatter;
	double sampleWidth = psfFWHM * (1 - pow(rfiScaleBW, 2));
	std::vector<std::future<std::vector<double>>> futureFiller;

	extraSteps = std::floor((psfFWHM / standardGap)*pow((rfiScaleBW), 2));
	futureFiller.reserve((2*extraSteps) * (2*extraSteps) + 1);
	if (extraSteps % 2 != 0)
	{
		extraSteps += 1;
	}

	for (int i = -extraSteps; i <= extraSteps; i++)
	{
		for (int j = -extraSteps; j <= extraSteps; j++)
		{
			deltaX = i * (standardGap)*(pow((1 / rfiScaleBW), 2) - 1); //In degrees
			deltaY = j * (standardGap)*(pow((1 / rfiScaleBW), 2) - 1);
			decHold = decPoint + deltaX;
			raHold = raPoint + deltaY;
			fluxVal = 999999;//WE ASSUME THESE POINTS ARE HIGHER THAN THE SCATTER
			correlatedScatter = -999999;
			distance = Tools::getGCDistance(decHold, raHold, decPoint, raPoint, partSetProcSSS.centerDecDeg)*toDeg;
		
			if (distance > sampleWidth * 2)
			{
				continue;
			}
			futureFiller.push_back(std::async(std::launch::async, &ProcessorRFI::rfiRemovalMulti, this, decHold, raHold, fluxVal, correlatedScatter));
			counter++;
		}
	}
	rfiSetFinalValues(futureFiller);

}
void ProcessorRFI::rfiSetFinalValues(std::vector<std::future<std::vector<double>>> & futureVec)
{
	int iTemp, jTemp;
	std::vector<double> results;
	for (int j = 0; j < futureVec.size(); j++)
	{
		results = futureVec[j].get();
		for (int k = 0; k < results.size() / 7; k++)
		{
			iTemp = results[7 * k + 0];
			jTemp = results[7 * k + 1];
			rfiWeights[iTemp][jTemp].push_back(results[7 * k + 2]);
			rfiValues[iTemp][jTemp].push_back(results[7 * k + 3]);
			rfiCounts[iTemp][jTemp].push_back(results[7 * k + 4]);
			rfiDists[iTemp][jTemp].push_back(results[7 * k + 5]);
			rfiNs[iTemp][jTemp].push_back(results[7 * k + 6]);
		}
	}
}
std::vector<double> ProcessorRFI::rfiRemovalMulti(double dec, double ra, double fluxVal, double correlatedScatter)
{
	int  holdI, holdJ, pointCandidates = 0;
	double avgSigma = 0, distance, toRad = M_PI / 180.0;
	std::vector<int> checks, possibles;
	std::vector<double> inRange, rfiFitSolution, results;
	pointCandidates = findPossSSS(dec, ra, rfiScaleDeg, possibles);
	inRange.reserve(pointCandidates * 6);
	checks.resize(pointCandidates, true);

	//COLLECT POSSIBLES AND STORE INFORMATION ABOUT PROXIMITY, SCATTER, AND FLUX
	for (int i = 0; i < pointCandidates; i++)
	{
		holdI = possibles[2 * i];
		holdJ = possibles[2 * i + 1];

		distance = Tools::getGCDistance(dec, ra, scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg)*toDeg;

		inRange.push_back(scans[holdI].getFlux(holdJ));
		inRange.push_back(pow(cos((M_PI*distance) / (rfiScaleDeg*2.0)), 2));
		inRange.push_back(scans[holdI].getDataDumps(holdJ));
		inRange.push_back(scatter2dSSS[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()]);
		inRange.push_back(holdI);
		inRange.push_back(holdJ);

		avgSigma += scatter2dSSS[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()];
	}

	avgSigma = (avgSigma / (double)pointCandidates);

	rfiRemovePoints(avgSigma, checks, inRange);	//FIT A LOCAL MODEL AND REJECT MOST DEVIANT POINT UNTIL WITHIN TOLERATED SCATTER
	rfiFitSolution = rfiFit(checks, inRange);	//FIT FINAL LOCAL MODEL
	results = rfiCollectResults(fluxVal, correlatedScatter, rfiFitSolution, checks, inRange);

	return results;
}
std::vector<double> ProcessorRFI::rfiCollectResults(double fluxVal, double correlatedScatter, std::vector<double> &rfiFitSolution, std::vector<int> & checks, std::vector<double> &inRange)
{
	bool countAbove, countBelow, rejectPoints, fluxAbove;
	double dd;
	double distSq, wDistSq, pointCount, wPointCount;
	double LMWeight, LMValue;//local model value
	double f = rfiFitSolution[0];
	double g = rfiFitSolution[1];
	double centerLMFlux = f * 1 + g;
	double jLMFlux;
	//double shapeParam, shapeParamTemp, localModelDumpSum = 0, localModelShapeSum = 0;

	//int normalLocalWeightCounter, normalLocalWeight; WHY IS THE normalLocalWeight AN INTEGER?

	double pw; //proximity weight -- cos(M_PI*dist/2*rfiScaleDeg)^2
	double dw; //dump weight -- dataDump
	double pwSum = 0, dwSum = 0;
	double pwAvg = 0, dwAvg = 0;
	double pwTemp, dwTemp;

	//int normalLocalWeightCounter, normalLocalWeight;
	double pwDelta;
	double pwDeltaCounter;

	double rfiCountHold = 0, pointDistance;
	double iDec, iRa, jDec, jRa, jFlux;
	int checksTrue = 0;
	int iHoldI, iHoldJ, jHoldI, jHoldJ;

	std::vector<double> results;

	for (int i = 0; i < checks.size(); i++)
	{
		if (checks[i])
		{
			pwSum += inRange[6 * i + 1];
			dwSum += inRange[6 * i + 2];
			checksTrue++;
		}
	}

	pwAvg = pwSum / (double)checksTrue;

	results.reserve((inRange.size() / 6) * 7);
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		pwDelta = 0;
		pwDeltaCounter = 0;

		iHoldI = inRange[6 * i + 4];
		iHoldJ = inRange[6 * i + 5];

		iDec = scans[iHoldI].getDec(iHoldJ);
		iRa = scans[iHoldI].getRa(iHoldJ);

		if (checks[i])
		{
			distSq = 0.0;
			wDistSq = 0.0;
			pointCount = 0.0;
			wPointCount = 0.0;
			rfiCountHold = 0.0;

			pw = inRange[6 * i + 1];
			LMValue = f * pw + g;

			for (int j = 0; j < inRange.size() / 6; j++)
			{
				jHoldI = inRange[6 * j + 4];
				jHoldJ = inRange[6 * j + 5];

				jDec = scans[jHoldI].getDec(jHoldJ);
				jRa = scans[jHoldI].getRa(jHoldJ);
				jFlux = scans[jHoldI].getFlux(jHoldJ);

				pwTemp = inRange[6 * j + 1];
				pwDelta += pow(pwTemp - pwAvg, 2);
				pwDeltaCounter++;

				jLMFlux = f * pwTemp + g;

				if (correlatedWeightMap)//We determine the first correction factor and the correlation scale for photometry here.
				{
					countAbove = (centerLMFlux > sqrt(2.0)*correlatedScatter) ? true : false;
					fluxAbove = (jLMFlux > sqrt(2.0)*correlatedScatter) ? true : false;

					dd = scans[jHoldI].getDataDumps(jHoldJ);
					pointDistance = Tools::getGCDistance(iDec, iRa, jDec, jRa, partSetProcSSS.centerDecDeg)*toDeg;//CHECK THAT fwhmRfi*toDeg is right units

					if (checks[j] && (pointDistance < rfiScaleDeg) && ((countAbove && fluxAbove) || (!countAbove && !fluxAbove))) {//Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are (Case 1: above OR Case 2: below) the 2d-scatter. 
						rfiCountHold += 1;
						if ((pointDistance != 0))
						{
							pointCount += 1.0;
							wPointCount += 1.0*dd;
							wDistSq += pow(pointDistance, 2)*dd;
						}
					}
				}
			}

			if (rfiCountHold != rfiCountHold || rfiCountHold == 0)
			{
				rfiCountHold = 1;
			}

			if (inRange.size() / 6 < 3)
			{
				LMWeight = 0;
				for (int k = 0; k < inRange.size() / 6; k++)
				{
					LMWeight += inRange[6 * i + 2];
				}
			}
			else
			{
				LMWeight = (dwSum / (1.0 + 2.0 * pow((pw - pwAvg), 2) / (pwDelta / pwDeltaCounter)));
			}

			results.push_back(iHoldI);
			results.push_back(iHoldJ);
			results.push_back(LMWeight);//w of local model (Apendix C)
			results.push_back(LMValue);//local model value


			results.push_back(rfiCountHold);//number of points that satisfy conditions (1), (2), and (3) in appendix D. (Used to divide local model weights that go into the global model weight)
			results.push_back(wDistSq);
			results.push_back(wPointCount);

		}
	}

	return results;
}
std::vector<std::vector<double> > ProcessorRFI::rfiBuildGlobalModels()
{
	double wSum, dSum, nSum, factor_corr = 3.0;
	std::vector<double> xVec, wVec, dubFiller;
	std::vector<double> flagsTemp;
	std::vector<std::vector<double> > rfiSubtracted, threadVecVec;
	rfiSubtracted.resize(scans.size(), dubFiller);
	std::vector<std::future<std::vector<std::vector<double>>>> futureFiller;
	for (int i = 0; i < scans.size(); i++)
	{
		rfiSubtracted[i].resize(scans[i].getSize(), NAN);
	}

	for (int i = 0; i < rfiValues.size(); i++)
	{
		futureFiller.resize(rfiValues[i].size());//add in

		for (int j = 0; j < rfiValues[i].size(); j++)
		{
			if (rfiValues[i][j].size() > 1)
			{
				xVec = rfiValues[i][j];
				wVec = rfiWeights[i][j];
				futureFiller[j] = std::async(std::launch::async, spawnNewRCRThread2, wVec, xVec);

			}
			else if (rfiValues[i][j].size() == 1)
			{
				rfiSubtracted[i][j] = rfiValues[i][j][0];
				this->scans[i].setGMWeight(j, rfiWeights[i][j][0] / rfiCounts[i][j][0]);
				this->scans[i].setSSSCorrelation(j, 2.0*factor_corr*sqrt(rfiDists[i][j][0] / std::max(rfiNs[i][j][0], 1.0)));
			}
		}

		for (int j = 0; j < rfiValues[i].size(); j++)
		{
			if (rfiValues[i][j].size() > 1)
			{
				threadVecVec = futureFiller[j].get();
				double value = threadVecVec[0][0];
				flagsTemp = threadVecVec[1];
				rfiSubtracted[i][j] = value;

				wSum = 0, dSum = 0.0, nSum = 0.0;
				for (int k = 0; k < rfiWeights[i][j].size(); k++)
				{
					wSum += rfiWeights[i][j][k] / rfiCounts[i][j][k];
					if (flagsTemp[k])
					{
						dSum += rfiDists[i][j][k];
						nSum += rfiNs[i][j][k];
					}
				}
				this->scans[i].setGMWeight(j, wSum);
				this->scans[i].setSSSCorrelation(j, 2.0*factor_corr*sqrt(dSum / std::max(nSum, 1.0)));
			}
		}
	}

	return rfiSubtracted;
}

std::vector<double> ProcessorRFI::getCentroidLocations()
{
	return centroidLocations;
}
std::vector<std::vector<double>> ProcessorRFI::getRFISubtracted()
{
	return rfiSubtracted;
}