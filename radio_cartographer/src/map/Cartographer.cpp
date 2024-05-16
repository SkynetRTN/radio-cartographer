#include "map\Cartographer.h"
#include "utils\Tools.h"
#include "utils\RCR.h"
#include <math.h>
#include <fstream>
#include <future>



static double toRad = M_PI / 180.0;
static double toDeg = 180.0 / M_PI;


//constructors
Cartographer::Cartographer()
{

}
Cartographer::Cartographer(Composite composite, MapParameters params)
{
	this->scans = composite.getScans();
	this->m10PlusProcessing = params.m10PlusProcessing;
	this->partSetVecSSS = composite.getPartSetVecSSS();
	this->partSetVecLSS = composite.getPartSetVecLSS();
	this->psfFWHM = composite.getPSFFWHM();
	this->compPartSetProcSSS = composite.getCompPartSetProcSSS();
	this->compPartSetProcLSS = composite.getCompPartSetProcLSS();
	this->classificationsSSS = composite.getClassificationsSSS();
	this->classificationsLSS = composite.getClassificationsLSS();
	this->processingWeightingFunction = params.processedWeightScale;
	this->resolution = psfFWHM * params.pixelSize;
	this->pCoordinate = composite.getProcessingCoordinate();
	this->mCoordinate = composite.getMappingCoordinate();
	this->standardGaps = composite.getStandardGaps();
	this->rfiScaleDeg = params.rfiScale*psfFWHM;
	this->correlatedWeightMap = params.correlatedWeightMap;
	this->LSSMapping = params.LSSMapping;
	this->SSSMapping = params.SSSMapping;
}

//function calls
Map Cartographer::generateProcMapsMulti()
{
	procMap.reserveGrid(compPartSetProcSSS.maxDec, compPartSetProcSSS.minDec, compPartSetProcSSS.maxRa, compPartSetProcSSS.minRa, resolution, compPartSetProcSSS.centerDecDeg, compPartSetProcSSS.centerRaDeg);

	//INTERPOLATES THE GRID FOR ALL PIXELS

	int totalPixelsDec = (int)(compPartSetProcSSS.maxDec - compPartSetProcSSS.minDec) / (resolution);
	int totalPixelsRa = (int)(compPartSetProcSSS.maxRa - compPartSetProcSSS.minRa) / (resolution);
	int centerPixelRa, centerPixelDec;

	if (totalPixelsRa % 2 == 0)
	{
		centerPixelRa = totalPixelsRa / 2;
	}
	else
	{
		centerPixelRa = (totalPixelsRa - 1) / 2 + 1;
	}
	if (totalPixelsDec % 2 == 0)
	{
		centerPixelDec = totalPixelsDec / 2;
	}
	else
	{
		centerPixelDec = (totalPixelsDec - 1) / 2 + 1;
	}

	// Output file
	std::ofstream outputFile;
	outputFile.open("Output.txt", std::ios_base::app);

	outputFile << "CENTERRA," << compPartSetProcSSS.centerRaDeg << ",Center ra (deg)\n";
	outputFile << "CENTERDE," << compPartSetProcSSS.centerDecDeg << ",Center dec (deg)\n";

	outputFile << "CENTERX," << centerPixelRa << ",Center pixel along the x-axis\n";
	outputFile << "CENTERY," << centerPixelDec << ",Center pixel along the y-axis\n";
	outputFile << "PIXLSIZE," << resolution / psfFWHM << ",Size of single pixel (beams)\n";
	outputFile.close();

	//std::cout << "Beginning Grid Interpolation\n";

	std::vector<std::future<std::vector<double>>> futureVec;
	int xPixel, yPixel, yTemp;
	int liveThreads = 0;
	int completedThreads = 0;
	int maxThreads = 8;
	int counter = 0;

	std::vector<double> results;
	std::vector<int> indices;

	for (double i = compPartSetProcSSS.minDec; i <= compPartSetProcSSS.maxDec; i += resolution)
	{
		xPixel = Tools::determinePixel(i, compPartSetProcSSS.minDec, resolution);
		completedThreads = 0;
		counter = 0;
		futureVec.resize((int)((compPartSetProcSSS.maxRa - compPartSetProcSSS.minRa) / resolution) + 1);
		for (double j = compPartSetProcSSS.minRa; j <= compPartSetProcSSS.maxRa; j += resolution)
		{
			yPixel = Tools::determinePixel(j, compPartSetProcSSS.minRa, resolution);

			futureVec[counter] = std::async(std::launch::async, &Cartographer::processedGridInterpolationMulti, this, i, j);
			liveThreads++;
			counter++;

			if ((liveThreads) >= maxThreads)
			{
				for (int k = completedThreads; k < counter; k++)
				{
					results = futureVec[k].get(); //4
					setMapValues(results[10], results[11], results);
				}
				completedThreads += liveThreads;// 5
				liveThreads = 0;
			}
		}

		for (double j = compPartSetProcSSS.minRa; j <= compPartSetProcSSS.maxRa; j += resolution)
		{
			yPixel = Tools::determinePixel(j, compPartSetProcSSS.minRa, resolution);

			if (yPixel >= completedThreads)
			{
				results = futureVec[yPixel].get();
				setMapValues(xPixel, yPixel, results);
			}
			liveThreads = 0;
		}

	}

	//PROCESSED PATH MAP
	double val;
	double a, b;
	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			a = Tools::min((double)procMap.getSize(0) - 1, round((1.0 / resolution)*(scans[i].getDec(j) - compPartSetProcSSS.minDec)));
			b = Tools::min((double)procMap.getSize(1) - 1, round((1.0 / resolution)*(scans[i].getRa(j) - compPartSetProcSSS.minRa)));

			procMap.setScanNumber(a, b, i);
			procMap.setIndexNumber(a, b, j);

			if (i % 2 == 0)
			{
				procMap.setProcPath(a, b, 1);
			}
			else
			{
				procMap.setProcPath(a, b, 0.5);
			}
		}
	}
	return procMap;
}
std::vector<double> Cartographer::processedGridInterpolationMulti(double i_0, double j_0)
{
	MapTypes mapType;

	double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI;
	bool contained, withinBounds;
	bool m10PlusCriteria = m10PlusProcessing;
	int holdI, holdJ, count = 0, pointCandidates = 0;
	int xPixel, yPixel;
	double distance = 0.0;
	double largestGap, w, wLocalSum;
	double* ret = NULL;
	double scaleMapVal, corrMapVal;
	double factor_w = 3.0, factor_w_sign;
	double largestGapSSS, largestGapLSS;
	double fluxSSS = 0, wSSS = 0, w2SSS = 0, scaleSSS = 0, corrSSS = 0;
	double fluxLSS = 0, wLSS = 0, w2LSS = 0, scaleLSS = 0, corrLSS = 0;

	std::vector<int> possibles;
	std::vector<int> intFiller;
	std::vector<double> dubFiller;
	std::vector<double> fromSolve;
	std::vector<std::vector<int> > scanNumbersSSS, scanNumbersLSS;
	std::vector<double> inRange;
	std::vector<double> results;

	SurfaceType modelType;

	int edge = 1;
	xPixel = Tools::determinePixel(i_0, compPartSetProcSSS.minDec, resolution);
	yPixel = Tools::determinePixel(j_0, compPartSetProcSSS.minRa, resolution);

	fromSolve.resize(4, 0.0);
	mapType = partSetVecSSS[0].mapType;
	results.resize(12, NAN);
	withinBounds = checkEdgeCriteria(j_0, i_0);

	if (withinBounds)
	{
		edge = 10;
	}

	if (SSSMapping)
	{
		inRange = collectPointsSSS(i_0, j_0);
		scanNumbersSSS = determineScanNumbersSSS(inRange);
		modelType = determineSurfaceType(scanNumbersSSS, mapType);
		largestGapSSS = getPixelThetaGapSSS(i_0, j_0, inRange);
		scaleSSS = Tools::max(processingWeightingFunction, Tools::min((4.0 / 3.0)*largestGapSSS / (psfFWHM*toRad), 1.0));
		corrSSS = getCorrMapValueSSS(i_0, j_0, inRange, scaleSSS);
		fromSolve = applySurfaceModel(modelType, false, m10PlusCriteria, inRange, largestGapSSS, i_0, j_0);

		fluxSSS = fromSolve[0];
		wSSS = fromSolve[1];
		w2SSS = fromSolve[2];

		if (modelType == M0)
		{
			if (factor_w < 0.0)
			{
				factor_w_sign = -factor_w;
			}
			else
			{
				factor_w_sign = factor_w;
			}
			corrSSS = factor_w_sign;
		}
	}
	if (LSSMapping)
	{

		inRange = collectPointsLSS(i_0, j_0);
		scanNumbersLSS = determineScanNumbersLSS(inRange);
		modelType = determineSurfaceType(scanNumbersLSS, mapType);
		largestGapLSS = getPixelThetaGapLSS(i_0, j_0, inRange);
		scaleLSS = Tools::max(processingWeightingFunction, Tools::min((4.0 / 3.0)*largestGapLSS / (psfFWHM*toRad), 1.0));
		corrLSS = getCorrMapValueLSS(i_0, j_0, inRange, scaleLSS);

		fromSolve = applySurfaceModel(modelType, true, false, inRange, largestGapLSS, i_0, j_0);
		fluxLSS = fromSolve[0];
		wLSS = fromSolve[1];
		w2LSS = fromSolve[2];
		if (modelType == M0)
		{
			if (factor_w < 0.0)
			{
				factor_w_sign = -factor_w;
			}
			else
			{
				factor_w_sign = factor_w;
			}
			corrLSS = factor_w_sign;
		}
	}

	results[0] = fluxSSS;
	results[1] = wSSS;
	results[2] = w2SSS;
	results[3] = corrSSS;
	results[4] = scaleSSS;

	results[5] = fluxLSS;
	results[6] = wLSS;
	results[7] = w2LSS;
	results[8] = corrLSS;
	results[9] = scaleLSS;

	results[10] = xPixel;
	results[11] = yPixel;

	if (withinBounds == false)
	{
		results[0] = NAN;
		results[1] = NAN;
		results[2] = NAN;

		results[5] = NAN;
		results[6] = NAN;
		results[7] = NAN;
	}

	return results;
}
void Cartographer::setMapValues(int xPixel, int yPixel, std::vector<double>& results)
{

	//results[0] = fluxSSS;
	//results[1] = wSSS;
	//results[2] = w2SSS;
	//results[3] = corrSSS;
	//results[4] = scaleSSS;

	//results[5] = fluxLSS;
	//results[6] = wLSS;
	//results[7] = w2LSS;
	//results[8] = corrLSS;
	//results[9] = scaleLSS;

	//results[10] = xPixel;
	//results[11] = yPixel;

	procMap.setSSSProcFlux(xPixel, yPixel, results[0]);
	procMap.setSSSWeight(xPixel, yPixel, results[1]);
	procMap.setSSSWeight2(xPixel, yPixel, results[2]);
	procMap.setSSSCorrelation(xPixel, yPixel, results[3]);
	procMap.setSSSScale(xPixel, yPixel, results[4]);

	procMap.setLSSProcFlux(xPixel, yPixel, results[5]);
	procMap.setLSSWeight(xPixel, yPixel, results[6]);
	procMap.setLSSWeight2(xPixel, yPixel, results[7]);
	procMap.setLSSCorrelation(xPixel, yPixel, results[8]);
	procMap.setLSSScale(xPixel, yPixel, results[9]);

	//procMap.setDec(xPixel, yPixel, results[10]);
	//procMap.setRa(xPixel, yPixel, results[11]);
}

//misc
int Cartographer::findPossSSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
{
	//TAKES IN INDICIES FOR SPECIFIC POINT
	std::vector<int> intFiller;
	toRet.resize(100, 0.0);

	//STORE DATAPOINT DEC AND RA
	double decTemp = i_0, RATemp = j_0, iHold, jHold;

	//CREATES COUNTERS TO DETERMINE WHAT PIXEL THE DATA POINT IS CURRENTLY ON
	int pointCounter = 0;
	PartitionSet partSetProc = compPartSetProcSSS;

	//BEFORE PRECEEDING UNDERSTAND CLASSIFY() IN SURVEY.CPP FOR DEFINITIONS ON subDecInc AND subDecRes
	//subDecInc IS THE NUMBER OF DEGREES NEEDED TO BE PAST TO BE IN A PARTITION SET RANGE
	int decCounter = (int)((decTemp - partSetProc.minDec) / partSetProc.subDecInc);
	int RACounter = (int)((RATemp - partSetProc.minRa) / partSetProc.subRaInc);

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

	//THE PARTITION SET BELOW CONTAINS THE GRID FOR ALL PIXEL COORDINATES
	//THE LOOPS ITERATES THROUGH THE 4 ADJACENT PARTITIONS IN BOTH RA AND DEC
	//SO LONG AS POINTS FALL WITHIN A CERTAIN RADIUS THEY ARE STORED IN A VECTOR
	double decHold, raHold;
	bool scan0 = false, scan1 = false, scan2 = false, scan3 = false;//
	for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProc.subDecRes); i++)
	{
		for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProc.subRaRes); j++)
		{
			for (int k = 0; k < classificationsSSS[i][j].size(); k += 2)
			{
				iHold = classificationsSSS[i][j][k];
				jHold = classificationsSSS[i][j][k + 1];

				decHold = scans[iHold].getDec(jHold);
				raHold = scans[iHold].getRa(jHold);

				if ((i_0 - limit <= decHold)
					&& (decHold <= i_0 + limit)
					&& (j_0 - limit <= raHold)
					&& (raHold <= j_0 + limit))
					/*
					&& (j_0*cos(i_0*toRad) - limit <= raHold*cos(i_0*toRad))
					&& (raHold*cos(i_0*toRad) <= j_0*cos(i_0*toRad) + limit))
					*/
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

	//RETURNS ALL POINTS WITHIN THE limit RADIUS OF THE POINT
	//STORES IN MEMORY ADDRESS OF toRet THE VALUES STORED ABOVE
	return pointCounter;
}
int Cartographer::findPossLSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
{
	//TAKES IN INDICIES FOR SPECIFIC POINT
	std::vector<int> intFiller;
	toRet.resize(100, 0.0);

	//STORE DATAPOINT DEC AND RA
	double decTemp = i_0, RATemp = j_0, iHold, jHold;

	//CREATES COUNTERS TO DETERMINE WHAT PIXEL THE DATA POINT IS CURRENTLY ON
	int pointCounter = 0;
	PartitionSet partSetProc = compPartSetProcLSS;

	//BEFORE PRECEEDING UNDERSTAND CLASSIFY() IN SURVEY.CPP FOR DEFINITIONS ON subDecInc AND subDecRes
	//subDecInc IS THE NUMBER OF DEGREES NEEDED TO BE PAST TO BE IN A PARTITION SET RANGE
	int decCounter = (int)((decTemp - partSetProc.minDec) / partSetProc.subDecInc);
	int RACounter = (int)((RATemp - partSetProc.minRa) / partSetProc.subRaInc);

	if (decCounter >= classificationsLSS.size())
	{
		decCounter = classificationsLSS.size() - 1;
	}
	if (RACounter >= classificationsLSS[decCounter].size())
	{
		RACounter = classificationsLSS[decCounter].size() - 1;
	}
	if (decCounter < 0)
	{
		decCounter = 0;
	}
	if (RACounter < 0)
	{
		RACounter = 0;
	}

	//THE PARTITION SET BELOW CONTAINS THE GRID FOR ALL PIXEL COORDINATES
	//THE LOOPS ITERATES THROUGH THE 4 ADJACENT PARTITIONS IN BOTH RA AND DEC
	//SO LONG AS POINTS FALL WITHIN A CERTAIN RADIUS THEY ARE STORED IN A VECTOR
	double decHold, raHold;
	bool scan0 = false, scan1 = false, scan2 = false, scan3 = false;//
	for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProc.subDecRes); i++)
	{
		for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProc.subRaRes); j++)
		{
			for (int k = 0; k < classificationsLSS[i][j].size(); k += 2)
			{
				iHold = classificationsLSS[i][j][k];
				jHold = classificationsLSS[i][j][k + 1];

				decHold = scans[iHold].getLSSDec(jHold);
				raHold = scans[iHold].getLSSRa(jHold);

				if ((i_0 - limit <= decHold)
					&& (decHold <= i_0 + limit)
					&& (j_0 - limit <= raHold)
					&& (raHold <= j_0 + limit))
					/*
					&& (j_0*cos(i_0*toRad) - limit <= raHold*cos(i_0*toRad))
					&& (raHold*cos(i_0*toRad) <= j_0*cos(i_0*toRad) + limit))
					*/
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

	//RETURNS ALL POINTS WITHIN THE limit RADIUS OF THE POINT
	//STORES IN MEMORY ADDRESS OF toRet THE VALUES STORED ABOVE
	return pointCounter;
}
double Cartographer::getPixelThetaGapSSS(double i_0, double j_0, std::vector<double> &inRange)
{
	int holdI, holdJ;
	double dataWeightScale, dataWeight, wySum = 0, wSum = 0, exponent;
	double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI, distance;
	for (int i = 0; i < inRange.size()/6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ),compPartSetProcSSS.centerDecDeg);
		dataWeightScale = scans[holdI].getThetaGap(holdJ);

		exponent = -2.329267*log(((dataWeightScale / 2.0) / (psfFWHM*toRad))) - 0.510417;
		dataWeight = pow(-1 * log(distance / (psfFWHM*toRad)), exponent);

		wSum += dataWeight;
		wySum += dataWeightScale * dataWeight;
	}

	return wySum / wSum;
}
double Cartographer::getPixelThetaGapLSS(double i_0, double j_0, std::vector<double> &inRange)
{
	int holdI, holdJ;
	double dataWeightScale, dataWeight, wySum = 0, wSum = 0, exponent;
	double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI, distance;
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getTSDec(holdJ), scans[holdI].getTSRa(holdJ), compPartSetProcSSS.centerDecDeg);
		dataWeightScale = scans[holdI].getLSSThetaGap(holdJ);

		exponent = -2.329267*log(((dataWeightScale / 2.0) / (psfFWHM*toRad))) - 0.510417;
		dataWeight = pow(-1 * log(distance / (psfFWHM*toRad)), exponent);

		wSum += dataWeight;
		wySum += dataWeightScale * dataWeight;
	}

	return wySum / wSum;
}
double Cartographer::getCorrMapValueSSS(double i_0, double j_0, std::vector<double> &inRange, double scaleMapVal)
{
	int holdI, holdJ;
	double factor_w = 3.0;
	double dataWeightScale, dataWeight, wySum = 0, wSum = 0, exponent;
	double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI, distance;
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), compPartSetProcSSS.centerDecDeg);
		dataWeightScale = scans[holdI].getThetaGap(holdJ);

		exponent = -2.329267*log(((dataWeightScale / 2.0) / (psfFWHM*toRad))) - 0.510417;
		dataWeight = pow(-1 * log(distance / (psfFWHM*toRad)), exponent);

		wSum += dataWeight;
		wySum += sqrt(pow(scans[holdI].getSSSCorrelation(holdJ) / psfFWHM, 2) + pow(factor_w*scaleMapVal, 2)) * dataWeight;
	}

	return wySum / wSum;
}
double Cartographer::getCorrMapValueLSS(double i_0, double j_0, std::vector<double> &inRange, double scaleMapVal)
{
	int holdI, holdJ;
	double factor_w = 3.0;
	double dataWeightScale, dataWeight, wySum = 0, wSum = 0, exponent;
	double toRad = M_PI / 180.0, toDeg = 180.0 / M_PI, distance;
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getTSDec(holdJ), scans[holdI].getTSRa(holdJ), compPartSetProcSSS.centerDecDeg);
		dataWeightScale = scans[holdI].getLSSThetaGap(holdJ);

		exponent = -2.329267*log(((dataWeightScale / 2.0) / (psfFWHM*toRad))) - 0.510417;
		dataWeight = pow(-1 * log(distance / (psfFWHM*toRad)), exponent);

		wSum += dataWeight;
		wySum += factor_w*sqrt(pow(scans[holdI].getLSSThetaWPrime(holdJ), 2) + pow(scaleMapVal, 2)) * dataWeight;
	}

	return wySum / wSum;
}
std::vector<double> Cartographer::collectPointsSSS(double i_0, double j_0)
{
	int pointCandidates, holdI, holdJ, count = 0;
	double distance;
	bool contained;
	std::vector<int> possibles, intFiller;
	std::vector<double> inRange;

	pointCandidates = findPossSSS(i_0, j_0, psfFWHM, possibles);
	inRange.resize(pointCandidates * 6);

	//THIS FOR LOOP SAVES ALL POINT CANDIDATES THAT FALL WITHIN ONE BEAMWIDTH
	for (int i = 0; i < pointCandidates; i++)
	{
		holdI = possibles[i * 2]; //SCAN NUMBER
		holdJ = possibles[i * 2 + 1]; //DATA POINT

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), compPartSetProcSSS.centerDecDeg)*toDeg;

		if (distance < psfFWHM)
		{
			//STORES ALL THE INFORMATION OF THE ACCEPTED DATA POINT
			inRange[count * 6] = (scans[holdI].getFlux(holdJ));
			inRange[count * 6 + 1] = cos((M_PI*distance) / (psfFWHM*2.0));
			inRange[count * 6 + 2] = (scans[holdI].getDec(holdJ));
			inRange[count * 6 + 3] = (scans[holdI].getRa(holdJ));
			inRange[count * 6 + 4] = holdI;
			inRange[count * 6 + 5] = holdJ;

			//distance = Tools::getGCDistance(scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), medianDec, medianRa);
			count++;
		}
	}
	inRange.erase(inRange.begin() + 6 * count, inRange.begin() + inRange.size());
	return inRange;
}
std::vector<double> Cartographer::collectPointsLSS(double i_0, double j_0)
{
	int pointCandidates, holdI, holdJ, count = 0;
	double distance;
	bool contained;
	std::vector<int> possibles, intFiller;
	std::vector<double> inRange;

	pointCandidates = findPossLSS(i_0, j_0, psfFWHM, possibles);
	inRange.resize(pointCandidates * 6);

	//THIS FOR LOOP SAVES ALL POINT CANDIDATES THAT FALL WITHIN ONE BEAMWIDTH
	for (int i = 0; i < pointCandidates; i++)
	{
		holdI = possibles[i * 2]; //SCAN NUMBER
		holdJ = possibles[i * 2 + 1]; //DATA POINT

		distance = Tools::getGCDistance(i_0, j_0, scans[holdI].getLSSDec(holdJ), scans[holdI].getLSSRa(holdJ), compPartSetProcSSS.centerDecDeg)*toDeg;

		if (distance < psfFWHM)
		{
			//STORES ALL THE INFORMATION OF THE ACCEPTED DATA POINT
			inRange[count * 6] = (scans[holdI].getLSSData(holdJ));
			inRange[count * 6 + 1] = cos((M_PI*distance) / (psfFWHM*2.0));
			inRange[count * 6 + 2] = (scans[holdI].getLSSDec(holdJ));
			inRange[count * 6 + 3] = (scans[holdI].getLSSRa(holdJ));
			inRange[count * 6 + 4] = holdI;
			inRange[count * 6 + 5] = holdJ;
			count++;
		}
	}
	inRange.erase(inRange.begin() + 6 * count, inRange.begin() + inRange.size());
	return inRange;
}
std::vector<std::vector<int> > Cartographer::determineScanNumbersSSS(std::vector<double> &inRange)
{
	std::vector<int> intFiller;
	std::vector<std::vector<int>> scanNumbersSSS;
	int holdI, holdJ;
	bool contained;
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];
		contained = false;
		//PUSH BACK VALUES OF SCAN NUMBER AND DATA POINT ON THAT SCAN INTO scanNumbers VECTOR
		for (int j = 0; j < scanNumbersSSS.size(); j++)
		{
			if (scanNumbersSSS[j][0] == holdI)
			{
				scanNumbersSSS[j].push_back(holdJ);
				contained = true;
			}
		}
		if (!contained)
		{
			scanNumbersSSS.push_back(intFiller);
			scanNumbersSSS[scanNumbersSSS.size() - 1].push_back(holdI);
			scanNumbersSSS[scanNumbersSSS.size() - 1].push_back(holdJ);
		}
	}
	return scanNumbersSSS;

}
std::vector<std::vector<int> > Cartographer::determineScanNumbersLSS(std::vector<double> &inRange)
{
	std::vector<int> intFiller;
	std::vector<std::vector<int>> scanNumbersLSS;
	int holdI, holdJ;
	bool contained;
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];
		contained = false;
		//PUSH BACK VALUES OF SCAN NUMBER AND DATA POINT ON THAT SCAN INTO scanNumbers VECTOR
		for (int j = 0; j < scanNumbersLSS.size(); j++)
		{
			if (scanNumbersLSS[j][0] == holdI)
			{
				scanNumbersLSS[j].push_back(holdJ);
				contained = true;
			}
		}
		if (!contained)
		{
			scanNumbersLSS.push_back(intFiller);
			scanNumbersLSS[scanNumbersLSS.size() - 1].push_back(holdI);
			scanNumbersLSS[scanNumbersLSS.size() - 1].push_back(holdJ);
		}
	}
	return scanNumbersLSS;
}


//edge calculations
bool Cartographer::checkEdgeCriteria(double ra, double dec)
{
	bool withinBounds = false;
	MapTypes mapType;

	std::vector<double> edgeOneParameters;
	std::vector<double> edgeTwoParameters;
	std::vector<double> edgeThreeParameters;
	std::vector<double> edgeFourParameters;

	std::vector<double> coordinateCheck;
	mapType = partSetVecSSS[0].mapType;
	bool tracking = compPartSetProcSSS.tracking;

	coordinateCheck = determineMappedEdges(ra, dec, tracking, compPartSetProcSSS, pCoordinate, mapType);

	ra = coordinateCheck[0];
	dec = coordinateCheck[1];

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
				withinBounds = true;
			}
		}
		else
		{
			double edgeRadius;
			double distance;
			double medianDec;
			double medianRa;
			double trimSize;
			edgeRadius = partSetVecSSS[i].edgeRadius;
			medianDec = partSetVecSSS[i].medianDec;
			medianRa = partSetVecSSS[i].medianRa;
			trimSize = partSetVecSSS[i].trimSize;

			distance = Tools::getModGCDistance(dec, ra, medianDec, medianRa)*toDeg;

			if ((distance / psfFWHM) < ((edgeRadius - trimSize*psfFWHM) / psfFWHM))
			{
				withinBounds = true;
			}
		}
	}
	return withinBounds;
}
std::vector<double> Cartographer::determineMappedEdges(double LONG, double LATI, bool tracking, PartitionSet compPartSetProc, Coordinates coordinate, MapTypes mapType)
{
	double centerProcLati, centerProcLong, centerMapLati, centerMapLong;
	double unTransLati, unTransLong, unTransLatiVert, unTransLongVert, doAngDeg, undoAngDeg;

	centerProcLati = compPartSetProc.centerLatProcDeg;
	centerProcLong = compPartSetProc.centerLongProcDeg;
	centerMapLati  = compPartSetProc.centerDecDeg;
	centerMapLong  = compPartSetProc.centerRaDeg;

	if (!(mCoordinate == "galactic" && pCoordinate == GALACTIC) && !(mCoordinate == "equatorial" && pCoordinate == EQUATORIAL)) // if mapping coords != processing coords.
	{
		if (pCoordinate == GALACTIC)
		{
			if (tracking)
			{
				centerProcLati = Tools::convertToB(compPartSetProc.medianLongMap, compPartSetProc.medianLatiMap);
				centerProcLong = Tools::convertToL(compPartSetProc.medianLongMap, compPartSetProc.medianLatiMap);
				doAngDeg = compPartSetProc.medianLatiMap;
				undoAngDeg = centerProcLati;
			}
			else
			{
				undoAngDeg = LATI + centerProcLati;
				doAngDeg = compPartSetProc.centerDecDeg;
			}
		}
		else if (pCoordinate == EQUATORIAL)
		{
			if (tracking)
			{
				centerProcLati = Tools::convertToDec(compPartSetProc.medianLongMap, compPartSetProc.medianLatiMap);
				centerProcLong = Tools::convertToRa(compPartSetProc.medianLongMap, compPartSetProc.medianLatiMap);
				doAngDeg = compPartSetProc.medianLatiMap;
				undoAngDeg = centerProcLati;
			}
			else
			{
				undoAngDeg = centerProcLati;
				doAngDeg = LATI;
			}
		}

		unTransLati = LATI + centerProcLati;	// these two lines do the un-transform in the processing coords.
		unTransLong = (LONG / cos(undoAngDeg*toRad)) + centerProcLong;

		if (pCoordinate == GALACTIC)
		{
			unTransLatiVert = Tools::convertToDec(unTransLong, unTransLati);
			unTransLongVert = Tools::convertToRa(unTransLong, unTransLati);
		}
		if (pCoordinate == EQUATORIAL)
		{
			unTransLatiVert = Tools::convertToB(unTransLong, unTransLati);
			unTransLongVert = Tools::convertToL(unTransLong, unTransLati);
		}
		if (mapType == NODDING)
		{
			LATI = unTransLatiVert;
			LONG = unTransLongVert;
		}
		else
		{
			LATI = unTransLatiVert - centerMapLati;	// these two lines re-transform in the mapping coords.
			LONG = (unTransLongVert - centerMapLong)*cos(doAngDeg*toRad);
		}
	}
	else // if mapping coords == processing coords.
	{
		if (mapType == NODDING) // edges for NODDINGs are now checked in un-transformed coordinates.
		{
			LATI = LATI + centerMapLati;
			LONG = (LONG / cos(centerMapLati*toRad)) + centerMapLong;
		}
	}

	std::vector<double> coordinatesToCheck;
	coordinatesToCheck.resize(2);
	coordinatesToCheck[1] = LATI;
	coordinatesToCheck[0] = LONG;

	return coordinatesToCheck;
}

//surface modeling
bool Cartographer::planeApplication(int minScanCount, int minPointsInScan, int minPointsInAp, std::vector<std::vector<int>> scanNumbers, int counter, MapTypes mapType)
{
	bool minApCheck = false, minScanCheck = false, minWScalePointsCheck = false;
	int pointCount;
	int possibles;
	int scanIndex, dataIndex;
	std::vector<int> possiblePoints;
	double wPointsInScan;

	if (mapType == DAISY && scans.size() == 4)
	{
		if (counter == 1)
		{
			minScanCount = 5;
			minPointsInScan = 5;
		}
		else if (counter == 2)
		{
			minScanCount = 4;
			minPointsInScan = 4;
		}
		else
		{
			minScanCount = 2;
			minPointsInScan = 2;
		}
	}

	if (scanNumbers.size() >= minScanCount)
	{
		pointCount = 0;
		for (int i = 0; i < scanNumbers.size(); i++)
		{
			if (scanNumbers[i].size() - 1 >= minPointsInScan)
			{
				minScanCheck = true;
			}
			pointCount += scanNumbers[i].size() - 1;
		}
		if (pointCount >= minPointsInAp)
		{
			minApCheck = true;
		}

		return minApCheck && minScanCheck;
	}
	return false;


	if (minApCheck && minScanCheck)
	{
		return true;
	}
	else
	{
		return false;
	}

}
SurfaceType Cartographer::determineSurfaceType(std::vector<std::vector<int>> & scanNumbers, MapTypes mapType)
{
	SurfaceType modelType;
	if (planeApplication(5, 5, 10, scanNumbers, 1, mapType))
	{
		modelType = M10;
	}
	else if (planeApplication(4, 4, 6, scanNumbers, 2, mapType))
	{
		modelType = M6;
	}
	else if (planeApplication(2, 2, 3, scanNumbers, 3, mapType))
	{
		modelType = M3;
	}
	else
	{
		modelType = M0;
	}
	return modelType;
}
std::vector<double> Cartographer::applySurfaceModel(SurfaceType modelType, bool largeStructure, bool m10PlusCriteria, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
{
	std::vector<double> values;
	switch (modelType)
	{
	case M10:
		values = surfaceFit10(m10PlusCriteria, inRange, largestGap, i_0, j_0);
		break;
	case M6:
		values = surfaceFit6(m10PlusCriteria, inRange, largestGap, i_0, j_0);
		break;
	case M3:
		values = surfaceFit3(m10PlusCriteria, inRange, largestGap, i_0, j_0);
		break;
	case M0:
		double w, wLocalSum;
		w = 0;
		wLocalSum = 0;
		for (int i = 0; i < inRange.size()/6; i++)
		{
			w += inRange[6 * i + 1];
			wLocalSum += inRange[6 * i + 1] * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]);
		}

		values.push_back(NAN);// fromSolve[0];
		values.push_back(w);// fromSolve[1];
		values.push_back(wLocalSum);// fromSolve[2];

		break;
	}
	return values;
}
std::vector<double> Cartographer::surfaceFit3(bool m10PlusCriteria, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
{
	int columnCount = 3, pivotIndex, intersect;
	int scanIndex, dataIndex, iPointCandidates;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower;
	double w = 0,
		xw = 0, yw = 0,
		xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, wLocalW = 0, wLocalW1 = 0, wLocalWInv = 0, finalAnswer;
	double ww = 0, sumWSquared = 0;
	double weightingFunction;

	weightingFunction = processingWeightingFunction;

	weightingFunctionTemp = Tools::max(weightingFunction, Tools::min((4.0 / 3.0)*largestGap / (psfFWHM*toRad), 1.0));
	double alpha = (log(0.5)) / (log(cos(M_PI*weightingFunctionTemp / 4.0)));

	//Photometry variables
	double rfiDistance, correlation = -999999, wcTop = 0.0, wcBottom = 0.0, wCorrelation;

	std::vector<double> A;
	std::vector<double> b;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);

	for (int i = 0; i < inRange.size() / 6; i++)
	{
		intersect = 0;

		zHold = inRange[6 * i];
		xHold = (inRange[6 * i + 2] - i_0);
		yHold = (inRange[6 * i + 3] - j_0);
		
		// CALCULATIONS FOR THE TWO INNER MOST CORRECTION FACTORS
		if (rfiScaleDeg != 0.0 && correlatedWeightMap == true)
		{
			correlation = scans[inRange[6 * i + 4]].getSSSCorrelation(inRange[6 * i + 5]); //in deg?

			for (int j = 0; j < inRange.size() / 6; j++)
			{
				rfiDistance = Tools::getGCDistance(inRange[6 * i + 2], inRange[6 * i + 3], inRange[6 * j + 2], inRange[6 * j + 3], compPartSetProcSSS.centerDecDeg)*toDeg;
				if (rfiDistance < correlation / 2.0) // refining photometry is as simple as tweaking this number.
				{
					intersect += 1;
				}
			}
			if (intersect == 0){
				intersect = 1;
			}
		}

		wHold = pow(inRange[1 + i * 6], alpha);
		ww += std::abs(pow(wHold, 2));
		if (rfiScaleDeg != 0 && correlatedWeightMap == true)
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]) / intersect;
			//wLocalW += pow(scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]), beta);
		}
		else
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5])*scans[inRange[6 * i + 4]].getDataDumps(inRange[6 * i + 5]);
		}
		

		w += wHold;
		xw += xHold*wHold;
		yw += yHold*wHold;
		xxw += xHold*xHold*wHold;
		yyw += yHold*yHold*wHold;
		xyw += xHold*yHold*wHold;
		zw += zHold*wHold;
		xzw += xHold*zHold*wHold;
		yzw += yHold*zHold*wHold;
	}

	A[0] = w; A[1] = xw; A[2] = yw;
	A[3] = xw; A[4] = xxw; A[5] = xyw;
	A[6] = yw; A[7] = xyw; A[8] = yyw;
	b[0] = zw; b[1] = xzw; b[2] = yzw;

	if (m10PlusCriteria == false && m10PlusProcessing == true)
	{
		A[0] = 2 * w;
	}

	finalHold = Tools::pivotSystem(columnCount, A, b);
	finalAnswer = finalHold[1][0] / finalHold[0][0];

	std::vector<double> toRet;
	toRet.resize(4, 0.0);
	if (finalAnswer < 0.0 && m10PlusCriteria == true)
	{
		std::vector<double> m3PlusValues;
		m3PlusValues.resize(4, 0.0);
		m10PlusCriteria = false;
		m3PlusValues= surfaceFit3(m10PlusCriteria, inRange, largestGap, i_0, j_0);

		toRet[0] = m3PlusValues[0];
		toRet[1] = m3PlusValues[1];
		toRet[2] = m3PlusValues[2];
		toRet[3] = m3PlusValues[3];
	}
	else
	{
		toRet[0] = finalAnswer;//CALCULATED FLUX FOR PIXEL
		toRet[1] = w;//WEIGHT FOR PIXEL
		toRet[2] = wLocalW;
		if (rfiScaleDeg != 0.0){
			toRet[3] = correlation;
		}
		else{
			toRet[3] = 0.0;
		}	
	}

	return toRet;
}
std::vector<double> Cartographer::surfaceFit6(bool m10PlusCriteria, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
{
	int columnCount = 6, pivotIndex;
	int scanIndex, dataIndex, iPointCandidates;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower;
	double w = 0, xw = 0, yw = 0, xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, finalAnswer;
	double ww = 0, sumWSquared = 0, wLocalW = 0, wLocalW1 = 0, wLocalWInv = 0;
	double weightingFunction, distanceToCenter;

	weightingFunction = processingWeightingFunction;

	weightingFunctionTemp = Tools::max(weightingFunction, Tools::min((4.0 / 3.0)*largestGap / (psfFWHM*toRad), 1.0));
	double alpha = (log(0.5)) / (log(cos(M_PI*weightingFunctionTemp / 4.0)));
	
	//Photometry variables
	int intersect = 0;
	double rfiDistance, correlation = -999999, wcTop = 0.0, wcBottom = 0.0, wCorrelation;

	std::vector<double> A;
	std::vector<double> b;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		intersect = 0;

		zHold = inRange[6 * i];
		xHold = (inRange[6 * i + 2] - i_0);
		yHold = (inRange[6 * i + 3] - j_0);

		// CALCULATIONS FOR THE TWO INNER MOST CORRECTION FACTORS
		if (rfiScaleDeg != 0.0 && correlatedWeightMap == true)
		{
			correlation = scans[inRange[6 * i + 4]].getSSSCorrelation(inRange[6 * i + 5]);
			for (int j = 0; j < inRange.size() / 6; j++)
			{
				rfiDistance = Tools::getGCDistance(inRange[6 * i + 2], inRange[6 * i + 3], inRange[6 * j + 2], inRange[6 * j + 3], compPartSetProcSSS.centerDecDeg)*toDeg;
				if (rfiDistance < correlation / 2.0) // refining photometry is as simple as tweaking this number.
				{
					intersect += 1;
				}
			}
			if (intersect == 0) {
				intersect = 1;
			}
		}

		wHold = pow(inRange[1 + i * 6], alpha);
		ww += std::abs(pow(wHold, 2));
		if (rfiScaleDeg != 0 && correlatedWeightMap == true)
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]) / intersect;
			//wLocalW += pow(scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]), beta);
		}
		else
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]) * scans[inRange[6 * i + 4]].getDataDumps(inRange[6 * i + 5]);
		}

		w += wHold;
		xw += xHold*wHold;
		yw += yHold*wHold;
		xxw += xHold*xHold*wHold;
		yyw += yHold*yHold*wHold;
		xyw += xHold*yHold*wHold;
		xxxw += xHold*xHold*xHold*wHold;
		xyyw += xHold*yHold*yHold*wHold;
		xxyw += xHold*xHold*yHold*wHold;
		yyyw += yHold*yHold*yHold*wHold;
		xxxxw += xHold*xHold*xHold*xHold*wHold;
		xxyyw += xHold*xHold*yHold*yHold*wHold;
		xxxyw += xHold*yHold*xHold*xHold*wHold;
		yyyyw += yHold*yHold*yHold*yHold*wHold;
		xyyyw += xHold*yHold*yHold*yHold*wHold;
		zw += zHold*wHold;
		xzw += xHold*zHold*wHold;
		yzw += yHold*zHold*wHold;
		xxzw += xHold*xHold*zHold*wHold;
		yyzw += yHold*yHold*zHold*wHold;
		xyzw += xHold*yHold*zHold*wHold;
	}

	sumWSquared = std::abs(pow(w, 2));

	A[0] = w; A[1] = xw; A[2] = yw; A[3] = xxw; A[4] = yyw; A[5] = xyw;
	A[6] = xw;  A[7] = xxw; A[8] = xyw; A[9] = xxxw; A[10] = xyyw; A[11] = xxyw;
	A[12] = yw; A[13] = xyw; A[14] = yyw; A[15] = xxyw; A[16] = yyyw; A[17] = xyyw;
	A[18] = xxw; A[19] = xxxw; A[20] = xxyw; A[21] = xxxxw; A[22] = xxyyw; A[23] = xxxyw;
	A[24] = yyw; A[25] = xyyw; A[26] = yyyw; A[27] = xxyyw; A[28] = yyyyw; A[29] = xyyyw;
	A[30] = xyw; A[31] = xxyw; A[32] = xyyw; A[33] = xxxyw; A[34] = xyyyw; A[35] = xxyyw;
	b[0] = zw; b[1] = xzw; b[2] = yzw; b[3] = xxzw; b[4] = yyzw; b[5] = xyzw;

	if (m10PlusCriteria == false && m10PlusProcessing == true)
	{
		A[0] = 2 * w;
	}

	finalHold = Tools::pivotSystem(columnCount, A, b);
	finalAnswer = finalHold[1][0] / finalHold[0][0];

	std::vector<double> toRet;
	toRet.resize(4, 0.0);
	if (finalAnswer < 0.0 && m10PlusCriteria == true)
	{
		std::vector<double> m6PlusValues;
		m6PlusValues.resize(4, 0.0);
		m10PlusCriteria = false;
		m6PlusValues = surfaceFit6(m10PlusCriteria, inRange, largestGap, i_0, j_0);

		toRet[0] = m6PlusValues[0];
		toRet[1] = m6PlusValues[1];
		toRet[2] = m6PlusValues[2];
		toRet[3] = m6PlusValues[3];
	}
	else
	{
		toRet[0] = finalAnswer;//CALCULATED FLUX FOR PIXEL
		toRet[1] = w;//WEIGHT FOR PIXEL
		toRet[2] = wLocalW;
		if (rfiScaleDeg != 0.0){
			toRet[3] = correlation;
		}
		else{
			toRet[3] = 0.0;
		}
	}

	return toRet;
}
std::vector<double> Cartographer::surfaceFit10(bool m10PlusCriteria, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
{
	int columnCount = 10, pivotIndex;
	int scanIndex, dataIndex, iPointCandidates;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower, alpha, distanceToCenter;
	double ww = 0, sumWSquared = 0, wLocalW = 0, wLocalW1 = 0, wLocalWInv = 0;
	double w = 0, xw = 0, yw = 0, xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		xxxxxw = 0, xxyyyw = 0, xxxxyw = 0, xxxyyw = 0, yyyyyw = 0, xyyyyw = 0,
		xxxxxxw = 0, xxxxxyw = 0, xxxxyyw = 0, xxxyyyw = 0, xxyyyyw = 0, xyyyyyw = 0, yyyyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0, xxxzw = 0, yyyzw = 0, xxyzw = 0, xyyzw = 0,
		xHold, yHold, zHold, wHold, finalAnswer;
	std::vector<double> A;
	std::vector<double> b;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	double weightingFunction;
	weightingFunction = processingWeightingFunction;
	weightingFunctionTemp = Tools::max(weightingFunction, Tools::min((4.0 / 3.0)*largestGap / (psfFWHM*toRad), 1.0));
	alpha = (log(0.5)) / (log(cos(M_PI*weightingFunctionTemp / 4.0)));
	
	//Photometry variables
	int intersect = 0;
	double rfiDistance, correlation = -999999, wcTop = 0.0, wcBottom = 0.0, wCorrelation;

	for (int i = 0; i < inRange.size() / 6; i++)
	{
		intersect = 0;

		zHold = inRange[6 * i];
		xHold = (inRange[6 * i + 2] - i_0);
		yHold = (inRange[6 * i + 3] - j_0);

		// CALCULATIONS FOR THE TWO INNER MOST CORRECTION FACTORS
		if (rfiScaleDeg != 0.0 && correlatedWeightMap == true)
		{
			correlation = scans[inRange[6 * i + 4]].getSSSCorrelation(inRange[6 * i + 5]);

			for (int j = 0; j < inRange.size() / 6; j++)
			{
				rfiDistance = Tools::getGCDistance(inRange[6 * i + 2], inRange[6 * i + 3], inRange[6 * j + 2], inRange[6 * j + 3], compPartSetProcSSS.centerDecDeg)*toDeg;
				if (rfiDistance < correlation / 2.0) // refining photometry is as simple as tweaking this number.
				{
					intersect += 1;
				}
			}
			if (intersect == 0) 
			{
				intersect = 1;
			}
		}

		wHold = pow(inRange[1 + i * 6], alpha);
		ww += std::abs(pow(wHold, 2));
		if (rfiScaleDeg != 0 && correlatedWeightMap == true)
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]) / intersect;
			//wLocalW += pow(scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]), beta);
		}
		else
		{
			wLocalW += wHold * scans[inRange[6 * i + 4]].getGMWeight(inRange[6 * i + 5]) * scans[inRange[6 * i + 4]].getDataDumps(inRange[6 * i + 5]);
		}


		w += wHold;
		xw += xHold*wHold;
		yw += yHold*wHold;
		xxw += xHold*xHold*wHold;
		yyw += yHold*yHold*wHold;
		xyw += xHold*yHold*wHold;
		xxxw += xHold*xHold*xHold*wHold;
		xyyw += xHold*yHold*yHold*wHold;
		xxyw += xHold*xHold*yHold*wHold;
		yyyw += yHold*yHold*yHold*wHold;
		xxxxw += xHold*xHold*xHold*xHold*wHold;
		xxyyw += xHold*xHold*yHold*yHold*wHold;
		xxxyw += xHold*yHold*xHold*xHold*wHold;
		yyyyw += yHold*yHold*yHold*yHold*wHold;
		xyyyw += xHold*yHold*yHold*yHold*wHold;
		xxxxxw += xHold*xHold*xHold*xHold*xHold*wHold;
		xxyyyw += xHold*xHold*yHold*yHold*yHold*wHold;
		xxxxyw += xHold*xHold*xHold*xHold*yHold*wHold;
		xxxyyw += xHold*xHold*xHold*yHold*yHold*wHold;
		yyyyyw += yHold*yHold*yHold*yHold*yHold*wHold;
		xyyyyw += xHold*yHold*yHold*yHold*yHold*wHold;
		xxxxxxw += xHold*xHold*xHold*xHold*xHold*xHold*wHold;
		xxxxxyw += xHold*xHold*xHold*xHold*xHold*yHold*wHold;
		xxxxyyw += xHold*xHold*xHold*xHold*yHold*yHold*wHold;
		xxxyyyw += xHold*xHold*xHold*yHold*yHold*yHold*wHold;
		xxyyyyw += xHold*xHold*yHold*yHold*yHold*yHold*wHold;
		xyyyyyw += xHold*yHold*yHold*yHold*yHold*yHold*wHold;
		yyyyyyw += yHold*yHold*yHold*yHold*yHold*yHold*wHold;
		zw += zHold*wHold;
		xzw += xHold*zHold*wHold;
		yzw += yHold*zHold*wHold;
		xxzw += xHold*xHold*zHold*wHold;
		yyzw += yHold*yHold*zHold*wHold;
		xyzw += xHold*yHold*zHold*wHold;
		xxxzw += xHold*xHold*xHold*zHold*wHold;
		yyyzw += yHold*yHold*yHold*zHold*wHold;
		xxyzw += xHold*xHold*yHold*zHold*wHold;
		xyyzw += xHold*yHold*yHold*zHold*wHold;
	}

	sumWSquared = std::abs(pow(w, 2));

	A[0] = w; A[1] = xw; A[2] = yw; A[3] = xxw; A[4] = yyw; A[5] = xyw; A[6] = xxxw; A[7] = yyyw; A[8] = xxyw; A[9] = xyyw;
	A[10] = xw; A[11] = xxw; A[12] = xyw; A[13] = xxxw; A[14] = xyyw; A[15] = xxyw; A[16] = xxxxw; A[17] = xyyyw; A[18] = xxxyw; A[19] = xxyyw;
	A[20] = yw; A[21] = xyw; A[22] = yyw; A[23] = xxyw; A[24] = yyyw; A[25] = xyyw; A[26] = xxxyw; A[27] = yyyyw; A[28] = xxyyw; A[29] = xyyyw;
	A[30] = xxw; A[31] = xxxw; A[32] = xxyw; A[33] = xxxxw; A[34] = xxyyw; A[35] = xxxyw; A[36] = xxxxxw; A[37] = xxyyyw; A[38] = xxxxyw; A[39] = xxxyyw;
	A[40] = yyw; A[41] = xyyw; A[42] = yyyw; A[43] = xxyyw; A[44] = yyyyw; A[45] = xyyyw; A[46] = xxxyyw; A[47] = yyyyyw; A[48] = xxyyyw; A[49] = xyyyyw;
	A[50] = xyw; A[51] = xxyw; A[52] = xyyw; A[53] = xxxyw; A[54] = xyyyw; A[55] = xxyyw; A[56] = xxxxyw; A[57] = xyyyyw; A[58] = xxxyyw; A[59] = xxyyyw;
	A[60] = xxxw; A[61] = xxxxw; A[62] = xxxyw; A[63] = xxxxxw; A[64] = xxxyyw; A[65] = xxxxyw; A[66] = xxxxxxw; A[67] = xxxyyyw; A[68] = xxxxxyw; A[69] = xxxxyyw;
	A[70] = yyyw; A[71] = xyyyw; A[72] = yyyyw; A[73] = xxyyyw; A[74] = yyyyyw; A[75] = xyyyyw; A[76] = xxxyyyw; A[77] = yyyyyyw; A[78] = xxyyyyw; A[79] = xyyyyyw;
	A[80] = xxyw; A[81] = xxxyw; A[82] = xxyyw; A[83] = xxxxyw; A[84] = xxyyyw; A[85] = xxxyyw; A[86] = xxxxxyw; A[87] = xxyyyyw; A[88] = xxxxyyw; A[89] = xxxyyyw;
	A[90] = xyyw; A[91] = xxyyw; A[92] = xyyyw; A[93] = xxxyyw; A[94] = xyyyyw; A[95] = xxyyyw; A[96] = xxxxyyw; A[97] = xyyyyyw; A[98] = xxxyyyw; A[99] = xxyyyyw;
	b[0] = zw; b[1] = xzw; b[2] = yzw; b[3] = xxzw; b[4] = yyzw; b[5] = xyzw; b[6] = xxxzw; b[7] = yyyzw; b[8] = xxyzw; b[9] = xyyzw;


	if (m10PlusCriteria == false && m10PlusProcessing == true) //ONLY OCCURS IN RECURSIVE CALL
	{
		A[0] = 2 * w;
	}

	finalHold = Tools::pivotSystem(columnCount, A, b);
	finalAnswer = finalHold[1][0] / finalHold[0][0];

	std::vector<double> toRet;
	toRet.resize(4, 0.0);

	if (finalAnswer < 0.0 && m10PlusCriteria == true)
	{
		std::vector<double> m10PlusValues;
		m10PlusValues.resize(4, 0.0);
		m10PlusCriteria = false;
		m10PlusValues = surfaceFit10(m10PlusCriteria, inRange, largestGap, i_0, j_0);

		toRet[0] = m10PlusValues[0];
		toRet[1] = m10PlusValues[1];
		toRet[2] = m10PlusValues[2];
		toRet[3] = m10PlusValues[3];
	}
	else
	{
		toRet[0] = finalAnswer;//CALCULATED FLUX FOR PIXEL
		toRet[1] = w;//WEIGHT FOR PIXEL
		toRet[2] = wLocalW;
		if (rfiScaleDeg != 0.0){
			toRet[3] = correlation;
		}
		else{
			toRet[3] = 0.0;
		}
	}

	return toRet;
}

Cartographer::~Cartographer()
{

}
