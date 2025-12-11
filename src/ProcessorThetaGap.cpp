#include "ProcessorThetaGap.h"
#include "Debugger.h"
#include "Tools.h"
#include "math.h"
#include <iostream>
#include <future>
#include <algorithm>

double toRad = M_PI / 180.0;
double toDeg = 180.0 / M_PI;
ProcessorThetaGap::ProcessorThetaGap()
{
}
ProcessorThetaGap::ProcessorThetaGap(std::vector<Scan> & scans, double psfFWHM, PartitionSet & partSetSSS, PartitionSet & partSetLSS, std::vector<std::vector<std::vector<int> > > & classSSS, std::vector<std::vector<std::vector<int> > > &classLSS)
{
	this->scans = scans;
	this->psfFWHM = psfFWHM;
	this->partSetProcSSS = partSetSSS;
	this->partSetProcLSS = partSetLSS;
	this->classificationsSSS = classSSS;
	this->classificationsLSS = classLSS;
}
int ProcessorThetaGap::findPossSSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
{
	std::vector<int> intFiller;
	toRet.resize(100, 0);
	double decTemp = i_0, RATemp = j_0, toRad = M_PI / 180.0;
	double distance;
	int pointCounter = 0, decCounter = (int)((decTemp - partSetProcSSS.minDec) / partSetProcSSS.subDecInc), RACounter = (int)((RATemp - partSetProcSSS.minRa) / partSetProcSSS.subRaInc), iHold, jHold;
	//int pointCounter = 0;
	//int decCounter = round(((decTemp - partSetProc.minDec) / partSetProc.subDecInc));
	//int RACounter = round((RATemp - partSetProc.minRa) / partSetProc.subRaInc);
	//int iHold, jHold;


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
	//for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProc.subDecRes); i++)
	for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProcSSS.subDecRes); i++)
	{
		//for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProc.subRaRes); j++)
		for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProcSSS.subRaRes); j++)
		{
			for (int k = 0; k < classificationsSSS[i][j].size(); k += 2)
			{
				iHold = classificationsSSS[i][j][k];
				jHold = classificationsSSS[i][j][k + 1];

				distance = Tools::getGCDistance(i_0, j_0, scans[iHold].getDec(jHold), scans[iHold].getRa(jHold), partSetProcSSS.centerDecDeg)*toDeg;

				//if ((i_0 - limit <= scans[iHold].getDec(jHold))
				//	&& (scans[iHold].getDec(jHold) <= i_0 + limit)
				//	&& (j_0 - limit <= scans[iHold].getRa(jHold))//THESE ARE ALREADY COSINE TRANSFORMED??
				//	&& (scans[iHold].getRa(jHold) <= j_0 + limit))

				//&& (j_0*cos(i_0*toRad) - limit <= scans[iHold].getRa(jHold)*cos(i_0*toRad))//THESE ARE ALREADY COSINE TRANSFORMED??
				//&& (scans[iHold].getRa(jHold)*cos(i_0*toRad) <= j_0*cos(i_0*toRad) + limit))

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
int ProcessorThetaGap::findPossLSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
{
	std::vector<int> intFiller;
	toRet.resize(100, 0);
	double decTemp = i_0, RATemp = j_0, toRad = M_PI / 180.0;
	//int pointCounter = 0, decCounter = (int)((decTemp - partSetProc.minDec) / partSetProc.subDecInc), RACounter = (int)((RATemp - partSetProc.minRa) / partSetProc.subRaInc), iHold, jHold;
	int pointCounter = 0;
	int decCounter = round(((decTemp - partSetProcLSS.minDec) / partSetProcLSS.subDecInc));
	int RACounter = round((RATemp - partSetProcLSS.minRa) / partSetProcLSS.subRaInc);
	int iHold, jHold;
	int index = round(limit / partSetProcLSS.subRaInc);


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
	//for (int i = Tools::max(decCounter - 1, 0); i < Tools::min((double)decCounter + 2, partSetProc.subDecRes); i++)
	for (int i = Tools::max(decCounter - index, 0); i < Tools::min((double)decCounter + index + 1, partSetProcLSS.subDecRes); i++)
	{
		//for (int j = Tools::max(RACounter - 1, 0); j < Tools::min((double)RACounter + 2, partSetProc.subRaRes); j++)
		for (int j = Tools::max(RACounter - index, 0); j < Tools::min((double)RACounter + index + 1, partSetProcLSS.subRaRes); j++)
		{
			for (int k = 0; k < classificationsLSS[i][j].size(); k += 2)
			{
				iHold = classificationsLSS[i][j][k];
				jHold = classificationsLSS[i][j][k + 1];
				if ((i_0 - limit <= scans[iHold].getDec(jHold))
					&& (scans[iHold].getDec(jHold) <= i_0 + limit)
					&& (j_0 - limit <= scans[iHold].getRa(jHold))//THESE ARE ALREADY COSINE TRANSFORMED??
					&& (scans[iHold].getRa(jHold) <= j_0 + limit))


					//&& (j_0*cos(i_0*toRad) - limit <= scans[iHold].getRa(jHold)*cos(i_0*toRad))//THESE ARE ALREADY COSINE TRANSFORMED??
					//&& (scans[iHold].getRa(jHold)*cos(i_0*toRad) <= j_0*cos(i_0*toRad) + limit))
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
std::vector<int> ProcessorThetaGap::quadrantSort(double dec, double ra, std::vector<int> &inRange, Quadrant quadrant)
{
	//NOMENCLATURE: 
	//RA AND DEC = "x1,y1"
	//RACHECK AND DECCHECK = "x2,y2"

	//THIS FUNCTION WHICH POINTS OF A GIVEN SET SIT WITHIN A CERTAIN QUADRANT.

	//IT USES RA AND DEC FOR "CARDINAL QUADRANTS" TOP, BOTTOM, LEFT, AND RIGHT, AND SLOPE FOR "DIAGONAL QUADRANTS" 


	double raCheck, decCheck, slope;
	double raNegLine, decNegLine, raPosLine, decPosLine;
	double raLine, raLine2;
	std::vector<int> subInRange;

	for (int k = 0; k < inRange.size() / 2; k++)
	{

		raCheck = scans[inRange[2 * k]].getRa(inRange[2 * k + 1]);
		decCheck = scans[inRange[2 * k]].getDec(inRange[2 * k + 1]);

		raLine = (decCheck - dec)*(-1) + ra;
		raLine2 = (decCheck - dec)*(1) + ra;

		if (quadrant == TOP)
		{
			if ((decCheck > dec) && (raCheck < raLine2) && (raCheck > raLine))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == BOTTOM)
		{
			if ((decCheck < dec) && (raCheck > raLine2) && (raCheck < raLine))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == LEFT_QUAD)
		{
			if ((raCheck > raLine2) && (raCheck > raLine))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == RIGHT_QUAD)
		{
			if ((raCheck < raLine2) && (raCheck < raLine))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == DIAG_TOP_RIGHT)
		{
			if ((decCheck > dec) && (raCheck < ra))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == DIAG_BOTTOM_RIGHT)
		{
			if ((decCheck < dec) && (raCheck < ra))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == DIAG_TOP_LEFT)
		{
			if ((decCheck > dec) && (raCheck > ra))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
		else if (quadrant == DIAG_BOTTOM_LEFT)
		{
			if ((decCheck < dec) && (raCheck > ra))
			{
				subInRange.push_back(inRange[2 * k]);
				subInRange.push_back(inRange[2 * k + 1]);
			}
		}
	}
	return subInRange;
}
void ProcessorThetaGap::determineCircleParams(double dec, double ra, Quadrant quadrant, std::vector<int> &inRange, std::vector<double> &circleParams)
{
	int scanIndex, dataIndex;

	double raCheck, decCheck;
	double radiusSquared, radius, diameter, distance;
	double x0, y0, x1, x2, y1, y2;
	double kDiag1, kDiag2;
	double decMin, decMax, raMin, raMax;

	kDiag1 = -1.0*1.0;
	kDiag2 = -1.0*-1.0;
	//THE FIRST -1.0 HERE HIGHLIGHTS THAT RA INCREASES RIGHT TO LEFT.  THEREFORE WHAT IS NORMALLY SEEN AS A NEGATIVE SLOPE, IS ACTUALLY A "POSITIVE" SLOPE IN THIS REFLECTION. 

	circleParams.resize(6 * inRange.size() / 2, 0.0);

	x1 = ra;
	y1 = dec;

	for (int k = 0; k < inRange.size() / 2; k++)
	{
		scanIndex = inRange[2 * k];
		dataIndex = inRange[2 * k + 1];

		x2 = scans[scanIndex].getRa(dataIndex);
		y2 = scans[scanIndex].getDec(dataIndex);

		raCheck = x2;
		decCheck = y2;

		if (quadrant == TOP || quadrant == BOTTOM)
		{
			y0 = (-1 * pow(x1, 2) + 2 * x1*x2 - pow(x2, 2) + pow(y1, 2) - pow(y2, 2)) / (2 * (y1 - y2));
			x0 = x1;
		}
		else if (quadrant == LEFT_QUAD || quadrant == RIGHT_QUAD)
		{
			x0 = (-1 * pow(x1, 2) + pow(x2, 2) + pow(y1, 2) - 2 * y1*y2 + pow(y2, 2)) / (-2 * x1 + 2 * x2);
			y0 = y1;
		}
		else if (quadrant == DIAG_TOP_RIGHT || quadrant == DIAG_BOTTOM_LEFT)
		{
			x0 = -1 * ((-1 * pow(x1, 2) + pow(x2, 2) + 2 * x1*(y1 - y2) + pow((y1 - y2), 2)) / (2 * (x1 - x2 - y1 + y2)));
			y0 = (pow(x1, 2) - 2 * x1*x2 + pow(x2, 2) + 2 * x1*y1 - 2 * x2*y1 - pow(y1, 2) + pow(y2, 2)) / (2 * (x1 - x2 - y1 + y2));
		}
		else if (quadrant == DIAG_BOTTOM_RIGHT || quadrant == DIAG_TOP_LEFT)
		{
			x0 = -1 * ((-1 * pow(x1, 2) + pow(x2, 2) - 2 * x1*(y1 - y2) + pow((y1 - y2), 2)) / (2 * (x1 - x2 + y1 - y2)));
			y0 = -1 * ((pow(x1, 2) + pow(x2, 2) + 2 * x2*y1 - pow(y1, 2) - 2 * x1*(x2 + y1) + pow(y2, 2))) / (2 * (x1 - x2 + y1 - y2));
		}

		radiusSquared = pow((x1 - x0), 2) + pow((y1 - y0), 2);
		if (radiusSquared != radiusSquared)
		{
			Debugger::print("Warn", "Non-sense radius calculated");
			radiusSquared = 0;
		}

		radius = std::abs(sqrt(radiusSquared));
		diameter = radius * 2;
		distance = Tools::getPythDistance(dec, ra, decCheck, raCheck);

		//CircleParameters: scanIndex, dataIndex, distance, radius, circleCenterRa, circleCenterDec

		circleParams[k * 6] = inRange[2 * k];
		circleParams[k * 6 + 1] = inRange[2 * k + 1];
		circleParams[k * 6 + 2] = distance;
		circleParams[k * 6 + 3] = diameter;
		circleParams[k * 6 + 4] = x0;
		circleParams[k * 6 + 5] = y0;
	}
}
double ProcessorThetaGap::maxGapQuadrant(std::vector<double> &circleParams)
{
	std::vector<int> scanIndices, dataIndices, innerPoints;
	std::vector<double> diameter, circleCenterRa, circleCenterDec, distance;
	int iIndex = 0, jIndex = 0, interiorCounter = 0;
	double raHold, decHold, distToCenter, thetaGap;
	bool thetaGapBool = false, largestInteriorPoint = false, samePoint = false;

	if (circleParams.size() == 0)
	{
		thetaGap = 999999;
		return thetaGap;
	}

	//CircleParameters: scanIndex, dataIndex, distance, radius, circleCenterRa, circleCenterDec

	scanIndices.resize(circleParams.size() / 6);
	dataIndices.resize(circleParams.size() / 6);
	distance.resize(circleParams.size() / 6);
	diameter.resize(circleParams.size() / 6);
	circleCenterRa.resize(circleParams.size() / 6);
	circleCenterDec.resize(circleParams.size() / 6);

	for (int i = 0; i < circleParams.size() / 6; i++)
	{
		scanIndices[i] = circleParams[6 * i];
		dataIndices[i] = circleParams[6 * i + 1];
		distance[i] = circleParams[6 * i + 2];
		diameter[i] = circleParams[6 * i + 3];
		circleCenterRa[i] = circleParams[6 * i + 4];
		circleCenterDec[i] = circleParams[6 * i + 5];
	}

	Tools::sortAll(0, diameter.size() - 1, diameter, circleCenterRa, circleCenterDec);
	Tools::sortAll(0, distance.size() - 1, distance, scanIndices, dataIndices);


	//ITERATE UNTIL POINT IS WITHIN A DIAMETER OF CIRCLE

	if (diameter.size() == 1)
	{
		thetaGap = diameter[0];
		return thetaGap;
	}

	while (distance[iIndex] > diameter[jIndex])
	{
		jIndex++;
		if (jIndex == diameter.size())
		{
			jIndex--;
			break;
		}
	}

	while (thetaGapBool == false)
	{
		//COLLECT INDICES OF ALL OF THE POINTS WITHIN THE DIAMETER
		for (int k = 0; k < distance.size(); k++)
		{
			if (diameter[jIndex] > distance[k])
			{
				innerPoints.push_back(k);
			}
		}

		//ONCE A POINT IS FOUND INSIDE OF A CIRCLE
		int counter = 0;
		double raPoint, decPoint;
		while (counter < innerPoints.size())
		{
			iIndex = innerPoints[counter];
			//CALCULATE THE DISTANCE FROM THE CENTER OF THE CIRCLE TO THE POINT

			raPoint = scans[scanIndices[iIndex]].getRa(dataIndices[iIndex]);
			decPoint = scans[scanIndices[iIndex]].getDec(dataIndices[iIndex]);
			distToCenter = Tools::getPythDistance(circleCenterDec[jIndex], circleCenterRa[jIndex], decPoint, raPoint) + 0.000002;//THIS IS TO HELP FOR WHEN THE DISTANCE TO A POINT IS EXACTLY THE RADIUS

																																 //IF IT FALLS INSIDE OF THE RADIUS, USE THE PREVIOUS DIAMETER
			if (distToCenter < diameter[jIndex] / 2)
			{
				thetaGap = diameter[jIndex - 1];
				thetaGapBool = true;
				break;
			}
			counter++;
		}
		if (thetaGapBool == false)
		{
			jIndex++;
			if (jIndex == diameter.size())
			{
				thetaGap = 999999;
				return thetaGap;
			}
		}
		innerPoints.clear();
	}
	return thetaGap;
}
double ProcessorThetaGap::circleThetaGap(int i_0, int j_0, bool LSSOn)
{
	std::vector<int> inRangeTemp, inRange, topRange, bottomRange, leftRange, rightRange, topRightDiagRange, topLeftDiagRange, bottomRightDiagRange, bottomLeftDiagRange;
	std::vector<int> scanIndices, dataIndices;
	std::vector<double> circlecircleParams, indicesVec;
	std::vector<double> circleParamsTop, circleParamsBottom, circleParamsLeft, circleParamsRight;
	std::vector<double> circleParamsTopRightDiag, circleParamsTopLeftDiag, circleParamsBottomRightDiag, circleParamsBottomLeftDiag, circleParamsAll;
	std::vector<double> diameter, circleCenterRa, circleCenterDec, distance;

	bool sameNumber;
	int sizeTop, sizeBottom, sizeLeft, sizeRight;
	int scanNum, dataNum, pointCandidates;
	double raHold, decHold, distToCenter;
	double topMax, bottomMax, leftMax, rightMax;
	double topRightDiagMax, topLeftDiagMax, bottomRightDiagMax, bottomLeftDiagMax;
	double topDiagsMax, bottomDiagsMax, diagsMax;
	double maxTopBottom, maxLeftRight, quadsMax, maxStandard, maxAll;
	double thetaGap, medianStandardGap;

	std::ofstream RangeFile;

	if (LSSOn)
	{
		raHold = scans[i_0].getLSSRa(j_0);
		decHold = scans[i_0].getLSSDec(j_0);
		pointCandidates = findPossLSS(decHold, raHold, psfFWHM, inRangeTemp);
	}
	else
	{
		raHold = scans[i_0].getRa(j_0);
		decHold = scans[i_0].getDec(j_0);
		pointCandidates = findPossSSS(decHold, raHold, psfFWHM, inRangeTemp);
	}



	//for (int i = 0; i < inRangeTemp.size() / 2; i++)
	//{
	//	sameNumber = false;
	//	if (inRangeTemp.size() == 0)
	//	{
	//		inRange.push_back(inRangeTemp[2 * i]);
	//		inRange.push_back(inRangeTemp[2 * i + 1]);
	//	}
	//	for (int j = 0; j < inRange.size() / 2; j++)
	//	{
	//		if (inRangeTemp[2 * i] == inRange[2 * j] && inRangeTemp[2 * i + 1] == inRange[2 * j + 1])
	//		{
	//			sameNumber = true;
	//		}
	//	}
	//	if (sameNumber == false)
	//	{
	//		inRange.push_back(inRangeTemp[2 * i]);
	//		inRange.push_back(inRangeTemp[2 * i + 1]);
	//	}
	//}


	double distanceCheck;
	for (int i = 0; i < inRangeTemp.size() / 2; i++)
	{
		distanceCheck = Tools::getGCDistance(decHold, raHold, scans[inRangeTemp[2 * i]].getDec(inRangeTemp[2 * i + 1]), scans[inRangeTemp[2 * i]].getRa(inRangeTemp[2 * i + 1]), partSetProcSSS.centerDecDeg)*toDeg;//WE NEED TO THINK ABOUT IF WE WANT TO TRANSFORM LSS GRID DIFFERENTLY
		sameNumber = false;
		if (inRangeTemp.size() == 0 && distanceCheck > scans[i_0].getRCRMinThetaGap(j_0))
		{
			inRange.push_back(inRangeTemp[2 * i]);
			inRange.push_back(inRangeTemp[2 * i + 1]);
		}
		for (int j = 0; j < inRange.size() / 2; j++)
		{
			if (inRangeTemp[2 * i] == inRange[2 * j] && inRangeTemp[2 * i + 1] == inRange[2 * j + 1])
			{
				sameNumber = true;
			}
		}
		if (sameNumber == false && distanceCheck > scans[i_0].getRCRMinThetaGap(j_0))
		{
			inRange.push_back(inRangeTemp[2 * i]);
			inRange.push_back(inRangeTemp[2 * i + 1]);
		}
	}


	//CARDINAL QUADRANTS, N,S,E,W
	topRange = quadrantSort(decHold, raHold, inRange, TOP);
	bottomRange = quadrantSort(decHold, raHold, inRange, BOTTOM);
	leftRange = quadrantSort(decHold, raHold, inRange, LEFT_QUAD);
	rightRange = quadrantSort(decHold, raHold, inRange, RIGHT_QUAD);

	determineCircleParams(decHold, raHold, TOP, topRange, circleParamsTop);
	determineCircleParams(decHold, raHold, BOTTOM, bottomRange, circleParamsBottom);
	determineCircleParams(decHold, raHold, LEFT_QUAD, leftRange, circleParamsLeft);
	determineCircleParams(decHold, raHold, RIGHT_QUAD, rightRange, circleParamsRight);


	////DIAGONAL QUADRANTS NE,SE,SW,NW
	topRightDiagRange = quadrantSort(decHold, raHold, inRange, DIAG_TOP_RIGHT);
	topLeftDiagRange = quadrantSort(decHold, raHold, inRange, DIAG_TOP_LEFT);
	bottomRightDiagRange = quadrantSort(decHold, raHold, inRange, DIAG_BOTTOM_RIGHT);
	bottomLeftDiagRange = quadrantSort(decHold, raHold, inRange, DIAG_BOTTOM_LEFT);

	determineCircleParams(decHold, raHold, DIAG_TOP_RIGHT, topRightDiagRange, circleParamsTopRightDiag);
	determineCircleParams(decHold, raHold, DIAG_TOP_LEFT, topLeftDiagRange, circleParamsTopLeftDiag);
	determineCircleParams(decHold, raHold, DIAG_BOTTOM_RIGHT, bottomRightDiagRange, circleParamsBottomRightDiag);
	determineCircleParams(decHold, raHold, DIAG_BOTTOM_LEFT, bottomLeftDiagRange, circleParamsBottomLeftDiag);


	//CALCULATE THE MAXIMUM THETA GAP
	topMax = maxGapQuadrant(circleParamsTop);
	bottomMax = maxGapQuadrant(circleParamsBottom);
	leftMax = maxGapQuadrant(circleParamsLeft);
	rightMax = maxGapQuadrant(circleParamsRight);

	topRightDiagMax = maxGapQuadrant(circleParamsTopRightDiag);
	topLeftDiagMax = maxGapQuadrant(circleParamsTopLeftDiag);
	bottomRightDiagMax = maxGapQuadrant(circleParamsBottomRightDiag);
	bottomLeftDiagMax = maxGapQuadrant(circleParamsBottomLeftDiag);

	//DETERMINE THE LARGEST THETA GAP
	topDiagsMax = Tools::max(topRightDiagMax, topLeftDiagMax);
	bottomDiagsMax = Tools::max(bottomRightDiagMax, bottomLeftDiagMax);
	diagsMax = Tools::max(topDiagsMax, bottomDiagsMax);

	maxTopBottom = Tools::max(topMax, bottomMax);
	maxLeftRight = Tools::max(leftMax, rightMax);
	quadsMax = Tools::max(maxTopBottom, maxLeftRight);

	//maxAll = quadsMax;
	maxAll = Tools::max(quadsMax, diagsMax);

	thetaGap = maxAll * toRad;

	//THE 4.0/3.0 CHECK IS APPLIED NOW IN THE THETA GAP SETTER 

	return thetaGap; //IN RADIANS
}

std::vector<std::vector<double>> ProcessorThetaGap::calculateThetaGapSSS()
{
	std::vector<std::vector<double>> thetaGapGrid;
	std::vector<double> thetaGapVec;
	std::vector<std::future<double>> futureVec;
	thetaGapGrid.resize(scans.size());
	double value;
	for (int i = 0; i < scans.size(); i++)
	{
		futureVec.resize(scans[i].getSize());
		thetaGapVec.resize(scans[i].getSize(), -999999);

		for (int j = 0; j < scans[i].getSize(); j++)
		{
			futureVec[j] = std::async(std::launch::async, &ProcessorThetaGap::circleThetaGap, this, i, j, false);
		}
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			value = futureVec[j].get();
			if (value > psfFWHM*toRad*(0.75))
			{
				value = psfFWHM * toRad*(0.75);
			}
			thetaGapVec[j] = value;
		}
		thetaGapGrid[i] = thetaGapVec;
		thetaGapVec.clear();
	}
	return thetaGapGrid;
}
std::vector<std::vector<double>> ProcessorThetaGap::calculateThetaGapLSS()
{
	std::vector<std::vector<double>> thetaGapGrid;
	std::vector<double> thetaGapVec;
	std::vector<std::future<double>> futureVec;
	thetaGapGrid.resize(scans.size());
	double value;
	for (int i = 0; i < scans.size(); i++)
	{
		futureVec.resize(scans[i].getLSSSize());
		thetaGapVec.resize(scans[i].getLSSSize(), -999999);
		for (int j = 0; j < scans[i].getLSSSize(); j++)
		{
			futureVec[j] = std::async(std::launch::async, &ProcessorThetaGap::circleThetaGap, this, i, j, true);
		}
		for (int j = 0; j < scans[i].getLSSSize(); j++)
		{
			value = futureVec[j].get();
			if (value > psfFWHM*toRad*(0.75))
			{
				value = psfFWHM * toRad*(0.75);
			}
			thetaGapVec[j] = value;
		}
		thetaGapGrid[i] = thetaGapVec;
		thetaGapVec.clear();
	}
	return thetaGapGrid;
}
ProcessorThetaGap::~ProcessorThetaGap()
{
}
/*
void ProcessorThetaGap::safeguardEdgeThetaGaps(Composite &composite, double weightScale)
{
	int counter = 0;
	int possibles;
	double decHold;
	double maxLocalThetaGap = -999999;
	std::vector<int> inRangeTemp;
	std::vector<double> dubFiller;
	std::vector<std::vector<double> > thetaGapFill;

	MapTypes scanMapType;

	for (int i = 0; i < scans.size(); i++)
	{
		thetaGapFill.push_back(dubFiller);
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			if (scans[i].getThetaGap(j) > psfFWHM*toRad*(0.75))
			{
				scans[i].setThetaGap(j, psfFWHM*toRad*(0.75));
			}
			thetaGapFill[i].push_back(scans[i].getThetaGap(j));
		}
	}

	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			scanMapType = composite.getMapTypes(scans[i].getSurveyNumber());
			if (scanMapType != DAISY)
			{
				if (((i == 0 || i == scans.size()) || ((j == 0 || j == scans[i].getSize() - 1)) ||
					(scans[i].getScanNumberInSurvey() == 0 || scans[i].getLastScan())))//CHECKS IF FIRST OR LAST SCAN IN SURVEY
				{
					possibles = findPossSSS(scans[i].getDec(j), scans[i].getRa(j), (4.0 / 3.0)*std::max(scans[i].getIntraScanGap(), scans[i].getInterScanGap()), inRangeTemp);//SCAN GAPS IN DEGREES, THETA GAPS IN RADIANS
					maxLocalThetaGap = -999999;
					for (int k = 0; k < possibles; k++)
					{
						if (maxLocalThetaGap < scans[inRangeTemp[2 * k]].getThetaGap(inRangeTemp[2 * k + 1]))
						{
							maxLocalThetaGap = scans[inRangeTemp[2 * k]].getThetaGap(inRangeTemp[2 * k + 1]);
						}
					}
					for (int k = 0; k < possibles; k++)
					{
						thetaGapFill[inRangeTemp[2 * k]][inRangeTemp[2 * k + 1]] = maxLocalThetaGap;
					}
				}
			}
			else
			{
				if ((j == 0 || j == scans[i].getSize() - 1))
				{
					possibles = findPossSSS(scans[i].getDec(j), scans[i].getRa(j), (4.0 / 3.0)*std::max(scans[i].getIntraScanGap(), scans[i].getInterScanGap()), inRangeTemp);//SCAN GAPS IN DEGREES, THETA GAPS IN RADIANS
					maxLocalThetaGap = -999999;
					for (int k = 0; k < possibles; k++)
					{
						if (maxLocalThetaGap < scans[inRangeTemp[2 * k]].getThetaGap(inRangeTemp[2 * k + 1]))
						{
							maxLocalThetaGap = scans[inRangeTemp[2 * k]].getThetaGap(inRangeTemp[2 * k + 1]);
						}
					}
					for (int k = 0; k < possibles; k++)
					{
						thetaGapFill[inRangeTemp[2 * k]][inRangeTemp[2 * k + 1]] = maxLocalThetaGap;
					}
				}
			}

		}
	}

	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			scans[i].setThetaGap(j, thetaGapFill[i][j]);
		}
	}

	composite.setScans(scans);
	this->scans = scans;
}
*/