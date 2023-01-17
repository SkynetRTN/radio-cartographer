#include "Processor.h"
#include "OutputFile.h"
#include "Debugger.h"
#include "Scan.h"
#include "Tools.h"
#include "RCR.h"
#include <math.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "Map.h"
#include "BackgroundCUDA.h"
#include "ProcessorRFI.h"
#include "ProcessorThetaGap.h"
#include "ProcessorTS.h"


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

std::vector<std::vector<double>> spawnNewRCRThread(std::vector<double> w, std::vector<double> y)
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


Processor::Processor()
{
	globalMaxDec = -999999;
	globalMinDec = 999999;
	globalMaxRa = -999999;
	globalMinRa = 999999;
}
Processor::Processor(Survey &survey)
{
	mapType = survey.getMapType();
	scans = survey.getScans();
	scansInRa = survey.getScanDirection();
	this->pCoordinate = survey.getProcessingCoordinate();
	this->mCoordinate = survey.getMappingCoordinate();
}
Processor::Processor(ProcessorParameters & procParams)
{
	globalMaxDec = -999999;
	globalMinDec = 999999;
	globalMaxRa = -999999;
	globalMinRa = 999999;
	bgScaleBW = procParams.bgScaleBW;
	rfiScaleBW = procParams.rfiScaleBW;
	wScaleBW = procParams.wScaleBW;
	timeShift = procParams.timeShift;
	timeShiftValue = procParams.timeShiftValue;
	wCorrMap = procParams.wCorrMap;
	lssProc = procParams.lssProc;

}

//processing commands
void Processor::performBGSubtractionMulti(Survey &survey)
{
	int counter = 0, completedThreads = 0, liveThreads = 0;
	int maxThreads = 8;
	std::vector<double> results;
	Output output;
	output.printBgSubtractionInfo(bgScaleBW);
	characterizeSurvey(survey);
	BackgroundCUDA bgCuda;
	std::vector<std::future<std::vector<double>>> futureVec;
	futureVec.resize(scans.size());

	Debugger::clockTime(true);

	if (bgScaleBW != 0.0)
	{
		for (int i = 0; i < scans.size(); i++)
		{
			Debugger::print("Info", i);
			bgCuda = BackgroundCUDA(scans[i], false);
			futureVec[i] = std::async(std::launch::async, &BackgroundCUDA::calculateBGMulti, bgCuda, bgScaleBW*psfFWHM);
			counter++;
			liveThreads++;

			if (liveThreads >= maxThreads)
			{
				for (int i = completedThreads; i < counter; i++)
				{
					results = futureVec[i].get();
					scans[i].removeBG(results);
				}
				completedThreads += liveThreads;
				liveThreads = 0;
			}
		}
		for (int i = completedThreads; i < scans.size(); i++)
		{
			results = futureVec[i].get();
			scans[i].removeBG(results);
		}
	}


	Debugger::clockTime(false);
	survey.setScans(scans);
}
void Processor::performTimeShifting(Survey &survey)
{
	characterizeSurvey(survey);

	double t_int = 0;
	if (timeShift == true)
	{
		if (timeShiftValue == 0.0) {
			ProcessorTS procTS = ProcessorTS(scans, mapType, psfFWHM, medianDiffAlongSweeps, scansInRa);
			t_int = procTS.find_tInt();
		}
		else {
			t_int = timeShiftValue // user has specified their own time shift value
		}
		survey.setTimeShift(t_int);
	}

	Output output;
	output.printTimeShiftInfo(timeShift, t_int);

	std::vector<double> holdDec, holdRa;
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].cosDecTransform(t_int, 0, 0, partSetProcSSS.centerDecDeg);
		scans[i].updateAngDistTemp(partSetProcSSS.centerDecDeg);
		holdDec = scans[i].getDec();
		holdRa = scans[i].getRa();
		scans[i].setTSDec(holdDec);
		scans[i].setTSRa(holdRa);
	}

	survey.setScans(scans);

	PartitionSet partSet = determineProcSurveyDimensions(survey, false);
	survey.setPartSetProcSSS(partSet);
	classifySSS(survey);

}
void Processor::processLargeScaleStructure(Survey &survey)
{
	scans = survey.getScans();
	partSetProcLSS = survey.getPartSetProcSSS();
	survey.setPartSetProcLSS(partSetProcLSS);
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setLSSStruct();
	}
	survey.setScans(scans);
	classifyLSS(survey);
	lssThetaGapCalc(survey);

	if (lssProc)
	{
		lssDropCorrection(survey);
		lssSet2DScatter(survey);
		lssPolynomialRFISubtraction(survey, bgScaleBW, wScaleBW);
		//lssElevationSubtraction(survey);
	}
}
void Processor::set2DScatter(Survey& survey)
{
	scans = survey.getScans();
	RCR rcr = RCR(LS_MODE_DL);
	int jBefore = 0, jAfter = 0;
	std::vector<double> regRet;
	double scatter2DSum, iterCount = 0, scatter2DAve;
	double minAngle1, maxAngle1, minAngle2, maxAngle2, minAngle3, maxAngle3, angleHold, minAngle, maxAngle, toPush, weight, weightLow, weightHigh, xBar, angleHigh, angleLow, angleMid, weightMid;
	std::vector<double> angleVec1, angleVec2, angleVec3, delta, weights, scatter2DAverages;
	std::vector<double> scatter2dTemp;
	scatter2dTemp.reserve(scans.size());
	std::ofstream scatter2DFile;
	scatter2DSum = 0;

	for (int i = 1; i < scans.size() - 1; i++)
	{
		minAngle1 = 999999;
		maxAngle1 = -999999;
		minAngle2 = 999999;
		maxAngle2 = -999999;
		minAngle3 = 999999;
		maxAngle3 = -999999;
		if (mapType == DAISY)
		{
			angleVec1 = Tools::daisyAngleBuilder(i - 1, scans[i - 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i - 1].getRa(), scans[i - 1].getDec());
			angleVec2 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg, scans[i].getRa(), scans[i].getDec());
			angleVec3 = Tools::daisyAngleBuilder(i + 1, scans[i + 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i + 1].getRa(), scans[i + 1].getDec());
		}
		else
		{
			if (scansInRa)
			{
				angleVec1 = scans[i - 1].getRa();
				angleVec2 = scans[i].getRa();
				angleVec3 = scans[i + 1].getRa();
			}
			else
			{
				angleVec1 = scans[i - 1].getDec();
				angleVec2 = scans[i].getDec();
				angleVec3 = scans[i + 1].getDec();
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
		for (int j = 0; j < angleVec3.size(); j++)
		{
			angleHold = angleVec3[j];
			if (angleHold < minAngle3)
			{
				minAngle3 = angleHold;
			}
			if (angleHold > maxAngle3)
			{
				maxAngle3 = angleHold;
			}
		}
		minAngle = minAngle1;
		if (minAngle2 > minAngle)
		{
			minAngle = minAngle2;
		}
		if (minAngle3 > minAngle)
		{
			minAngle = minAngle3;
		}
		maxAngle = maxAngle1;
		if (maxAngle2 < maxAngle)
		{
			maxAngle = maxAngle2;
		}
		if (maxAngle3 < maxAngle)
		{
			maxAngle = maxAngle3;
		}
		if (angleVec2[angleVec2.size() - 1] > angleVec2[0])
		{
			int j = 0, k = 0;
			jBefore = angleVec1.size() - 1, jAfter = angleVec3.size() - 1;
			while (minAngle > angleVec2[j])
			{
				j++;
			}
			while (maxAngle > angleVec2[k])
			{
				k++;
			}
			for (; j < k; j++)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore--;
				}
				if (jBefore < angleVec1.size() - 1 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore + 1] - angleVec2[j]))
				{
					jBefore++;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter--;
				}
				if (jAfter < angleVec3.size() - 1 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter + 1] - angleVec2[j]))
				{
					jAfter++;
				}
				if (scansInRa)
				{
					toPush = scans[i].getFlux(j) - (scans[i - 1].getFlux(jBefore) + (scans[i + 1].getFlux(jAfter) - scans[i - 1].getFlux(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getDec(jAfter);
					angleMid = scans[i].getDec(j);
					angleLow = scans[i - 1].getDec(jBefore);
				}
				else
				{
					toPush = scans[i].getFlux(j) - (scans[i - 1].getFlux(jBefore) + (scans[i + 1].getFlux(jAfter) - scans[i - 1].getFlux(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getRa(jAfter);
					angleMid = scans[i].getRa(j);
					angleLow = scans[i - 1].getRa(jBefore);
				}

				weightLow = scans[i - 1].getDataDumps(jBefore);
				weightHigh = scans[i + 1].getDataDumps(jAfter);
				xBar = (weightLow * angleLow + weightHigh * angleHigh) / (weightLow + weightHigh);
				weight = (1.0 / scans[i].getDataDumps(j));
				weight += pow(((angleHigh - xBar) * (angleHigh - angleLow)), 2.0) / weightLow;
				weight += pow(((angleLow - xBar) * (angleHigh - angleLow)), 2.0) / weightHigh;
				weight += pow(((angleMid - xBar) * (angleHigh - angleLow)), 2.0) * ((1.0 / weightHigh) + (1.0 / weightLow));
				weights.push_back(1.0 / weight);

			}
			rcr.performBulkRejection(delta);
			scatter2dTemp.push_back(std::sqrt(rcr.result.mu * rcr.result.mu + rcr.result.sigma * rcr.result.sigma));
			delta.clear();
			weights.clear();
		}
		else
		{
			int j = angleVec2.size() - 1, k = angleVec2.size() - 1;
			jBefore = 0, jAfter = 0;
			while (maxAngle > angleVec2[k])
			{
				k--;
			}
			while (minAngle > angleVec2[j])
			{
				j--;
			}
			for (; j > k; j--)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore++;
				}
				if (jBefore > 0 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore - 1] - angleVec2[j]))
				{
					jBefore--;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter++;
				}
				if (jAfter > 0 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter - 1] - angleVec2[j]))
				{
					jAfter--;
				}
				if (scansInRa)
				{
					toPush = scans[i].getFlux(j) - (scans[i - 1].getFlux(jBefore) + (scans[i + 1].getFlux(jAfter) - scans[i - 1].getFlux(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getDec(jAfter);
					angleMid = scans[i].getDec(j);
					angleLow = scans[i - 1].getDec(jBefore);
				}
				else
				{
					toPush = scans[i].getFlux(j) - (scans[i - 1].getFlux(jBefore) + (scans[i + 1].getFlux(jAfter) - scans[i - 1].getFlux(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getRa(jAfter);
					angleMid = scans[i].getRa(j);
					angleLow = scans[i - 1].getRa(jBefore);
				}

				weightLow = scans[i - 1].getDataDumps(jBefore);
				weightHigh = scans[i + 1].getDataDumps(jAfter);
				xBar = (weightLow * angleLow + weightHigh * angleHigh) / (weightLow + weightHigh);
				weight = (1.0 / scans[i].getDataDumps(j));
				weight += pow(((angleHigh - xBar) * (angleHigh - angleLow)), 2.0) / weightLow;
				weight += pow(((angleLow - xBar) * (angleHigh - angleLow)), 2.0) / weightHigh;
				weight += pow(((angleMid - xBar) * (angleHigh - angleLow)), 2.0) * ((1.0 / weightHigh) + (1.0 / weightLow));
				weights.push_back(1.0 / weight);

			}
			rcr.performBulkRejection(delta);
			scatter2dTemp.push_back(std::sqrt(rcr.result.mu * rcr.result.mu + rcr.result.sigma * rcr.result.sigma));
			delta.clear();
			weights.clear();
		}
	}

	std::vector<double> indices, dumpSums;

	for (int i = 1; i < scans.size() - 1; i++)
	{
		dumpSums.push_back(scans[i].getDumpSum());
		indices.push_back(i);

	}
	for (int i = 0; i < scatter2dTemp.size(); i++)
	{
		scatter2dTemp[i] *= .8137;
	}
	RCR rcr2 = RCR(LS_MODE_68);
	LinearModel model = LinearModel(dumpSums, indices, scatter2dTemp);
	rcr2.setParametricModel(model);
	rcr2.performBulkRejection(dumpSums, scatter2dTemp);

	bool useAverage = false;
	double average2dScatter;

	if (scatter2dTemp.size() == 2)
	{
		useAverage = true;
		average2dScatter = (scatter2dTemp[0] + scatter2dTemp[1]) / 2.0;
	}

	double m, b;
	m = model.m;
	b = model.b;
	xBar = model.xBar;
	scatter2dTemp.clear();
	scatter2dTemp.resize(scans.size(), 0.0);

	for (int i = 0; i < scans.size(); i++)
	{
		if (useAverage)
		{
			scatter2dTemp[i] = average2dScatter;
		}
		else
		{
			scatter2dTemp[i] = m * (i - xBar) + b;
		}
		scatter2dTemp[i] = std::max(scans[i].getScatter(), scatter2dTemp[i]);
	}

	survey.set2DScatterVec(scatter2dTemp);
}
void Processor::performRFIRejectionMulti(Composite &composite, std::vector<Survey> &surveys)
{
	Debugger::print("Info", "Starting RFI");
	PartitionSet partSet;
	partSet = determineProcCompositeDimensions(composite, false);
	composite.setCompPartSetProcSSS(partSet);
	partSet = determineProcCompositeDimensions(composite, true);
	composite.setCompPartSetProcLSS(partSet);
	classifySSS(composite);
	classifyLSS(composite);

	scans = composite.getScans();
	partSetProcSSS = composite.getCompPartSetProcSSS();
	classificationsSSS = composite.getClassificationsSSS();
	scatter2dSSS = composite.getScatter2D();
	standardGaps = composite.getStandardGaps();

	double rfiScaleDeg = rfiScaleBW * psfFWHM;
	double rfiScaleRad = rfiScaleDeg * M_PI / 180.0;
	
	std::vector<MapTypes> mapHolder;
	mapHolder.resize(surveys.size());
	for (int i = 0; i < surveys.size(); i++)
	{
		mapHolder[i] = surveys[i].getMapType();
	}

	double standardGap = 999999;
	double standardGapTemp;
	for (int i = 0; i < standardGaps.size(); i++)
	{
		standardGapTemp = standardGaps[i];
		if (mapHolder[i] == DAISY)
		{
			standardGapTemp = (standardGapTemp / 2)*scans.size();
		}
		if (standardGapTemp < standardGap)
		{
			standardGap = standardGapTemp;
		}
	}

	Output output;
	output.printRfiInfo(rfiScaleBW);
	std::vector<std::vector<double>> rfiSubtracted;

	if (rfiScaleBW != 0.0)
	{
		RFIParameters rfiParams;
		rfiParams.correlatedWeightMap = wCorrMap;
		rfiParams.psfFWHM = psfFWHM;
		rfiParams.rfiScaleBW = rfiScaleBW;
		rfiParams.standardGap = standardGap;
		rfiParams.partSetProcSSS = partSetProcSSS;

		ProcessorRFI processorRFI(scans, rfiParams, scatter2dSSS, classificationsSSS);
		rfiSubtracted = processorRFI.getRFISubtracted();
		centroidLocations = processorRFI.getCentroidLocations();

		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].removeRFI(rfiSubtracted[i]); //resets the working channel
												  //scans[i].removePoints(flags[i]); // removes flagged points (this is also used in trimming edges, so we might be able to just rewrite this as a clip RFI Data Function that takes in the RFISubtractedVector
			scans[i].updateAngDistTemp(partSetProcSSS.centerDecDeg);
		}
	}

	composite.setScans(scans);
	repartitionScans(composite, surveys);
	composite.setScans(scans);

	partSet = determineProcCompositeDimensions(composite, false);
	composite.setCompPartSetProcSSS(partSet);
	classifySSS(composite);
}
void Processor::calculateProcThetaGapMulti(Composite & composite)
{
	Debugger::print("Info", "Starting Theta Gap");
	std::vector<double> thetaGapVec;
	std::vector<std::vector<double>> thetaGapGridSSS, thetaGapGridLSS;
	std::vector<std::future<double>> futureVec;
	double value;
	scans = composite.getScans();
	classificationsSSS = composite.getClassificationsSSS();
	classificationsLSS = composite.getClassificationsLSS();
	
	ProcessorThetaGap processorTG(scans, psfFWHM, partSetProcSSS, partSetProcLSS, classificationsSSS, classificationsLSS);
	
	if (LSSProcessing)
	{
		thetaGapGridSSS = processorTG.calculateThetaGapSSS();
		thetaGapGridLSS = processorTG.calculateThetaGapLSS();
		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setThetaGap(thetaGapGridSSS[i]);
			scans[i].setLSSThetaGap(thetaGapGridLSS[i]);
		}
	}
	else
	{
		thetaGapGridSSS = processorTG.calculateThetaGapSSS();
		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setThetaGap(thetaGapGridSSS[i]);
		}
	}

	/*
	for (int i = 0; i < scans.size(); i++)
	{
		futureVec.resize(scans[i].getSize());
		thetaGapVec.resize(scans[i].getSize(), -999999);

		for (int j = 0; j < scans[i].getSize(); j++)
		{
			futureVec[j] = std::async(std::launch::async, &Processor::circleThetaGap, this, i, j, false);
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
		scans[i].setThetaGap(thetaGapVec);
		thetaGapVec.clear();
	}
	if (LSSProcessing)
	{
		for (int i = 0; i < scans.size(); i++)
		{
			futureVec.resize(scans[i].getLSSSize());
			thetaGapVec.resize(scans[i].getLSSSize(), -999999);
			for (int j = 0; j < scans[i].getLSSSize(); j++)
			{
				futureVec[j] = std::async(std::launch::async, &Processor::circleThetaGap, this, i, j, LSSProcessing);
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
			scans[i].setLSSThetaGap(thetaGapVec);
			thetaGapVec.clear();
		}
	}
	*/
	composite.setScans(scans);
}

//data pre-processing
void Processor::characterizeSurvey(Survey &survey)
{
	mapType = survey.getMapType();
	scans = survey.getScans();
	psfFWHM = survey.getPSFFWHM();
	medianDiffAlongSweeps = survey.getMedianDiffAlongSweeps();
	scansInRa = survey.getScanDirection();
}
void Processor::determineDataBoundaries(Survey &survey)
{
	scans = survey.getScans();
	mapType = survey.getMapType();
	bool scansInRaTemp = survey.getScanDirection();

	RCR rcr = RCR(SS_MEDIAN_DL);

	double scanDecMin = 999999, scanDecMax = -999999;
	double scanRaMin = 999999, scanRaMax = -999999;
	double decVal, raVal;
	std::vector<double> minDecVals, maxDecVals;
	std::vector<double> minRaVals, maxRaVals;
	std::vector<double> decHold, decAll;
	std::vector<double> raHold, raAll;

	for (int i = 0; i < scans.size(); i++)
	{
		decHold = scans[i].getDec();
		raHold = scans[i].getRa();
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			decAll.push_back(decHold[j]);
			raAll.push_back(raHold[j]);
		}
	}

	//IF OVERLAYING IMAGES, EACH SURVEY SHOULD BE CENTERED ABOUT 0 WHEN DEC TRANSFORMED
	rcr.performBulkRejection(decAll);
	survey.setMedianDec(rcr.result.mu);

	rcr.performBulkRejection(raAll);
	survey.setMedianRa(rcr.result.mu);

	std::cout << "scan: " << scansInRaTemp << std::endl;

	if (mapType != DAISY)
	{
		if (scansInRaTemp == true)
		{
			minDecVals = scans[0].getDec();
			maxDecVals = scans[scans.size() - 1].getDec();
			for (int i = 0; i < scans.size(); i++)
			{
				scanRaMin = 999999;
				scanRaMax = -999999;
				for (int j = 0; j < scans[i].getSize(); j++)
				{
					raVal = scans[i].getRa(j);
					if (raVal > scanRaMax)
					{
						scanRaMax = raVal;
					}
					if (raVal < scanRaMin)
					{
						scanRaMin = raVal;
					}
				}
				minRaVals.push_back(scanRaMin);
				maxRaVals.push_back(scanRaMax);
			}

		}
		else
		{
			minRaVals = scans[0].getRa();
			maxRaVals = scans[scans.size() - 1].getRa();
			for (int i = 0; i < scans.size(); i++)
			{
				scanDecMin = 999999;
				scanDecMax = -999999;
				for (int j = 0; j < scans[i].getSize(); j++)
				{
					decVal = scans[i].getDec(j);
					if (decVal > scanDecMax)
					{
						scanDecMax = decVal;
					}
					if (decVal < scanDecMin)
					{
						scanDecMin = decVal;
					}
				}
				minDecVals.push_back(scanDecMin);
				maxDecVals.push_back(scanDecMax);
			}

		}

		std::cout << minDecVals.size() << std::endl;
		std::cout << maxDecVals.size() << std::endl;
		std::cout << minRaVals.size() << std::endl;
		std::cout << maxRaVals.size() << std::endl;

		rcr.performBulkRejection(minDecVals);
		if (globalMinDec > rcr.result.mu)
		{
			globalMinDec = rcr.result.mu;
		}

		rcr.performBulkRejection(maxDecVals);
		if (globalMaxDec < rcr.result.mu)
		{
			globalMaxDec = rcr.result.mu;
		}

		rcr.performBulkRejection(minRaVals);
		if (globalMinRa > rcr.result.mu)
		{
			globalMinRa = rcr.result.mu;
		}

		rcr.performBulkRejection(maxRaVals);
		if (globalMaxRa < rcr.result.mu)
		{
			globalMaxRa = rcr.result.mu;
		}
	}
	else
	{
		std::vector<double> radii;
		for (int i = 0; i < scans.size(); i++)
		{
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				//scans[i].cosDecTransform(0, survey.getMedianDec(), survey.getMedianRa(), survey.getMedianDec());
				radii.push_back(Tools::getGCDistance(scans[i].getDec(j), scans[i].getRa(j)*cos(survey.getMedianDec()), survey.getMedianDec(), survey.getMedianRa(), survey.getMedianDec())*toDeg);//THIS IS WRONG: DAISYS ARE DAISIES WHEN IN RA*C0S(DEC) SPACE, BUT WHEN STILL RAW RA AND DEC, THEY ARE ELONGATED ALONG THE RA AXIS, SO WE ARE TRYING TO UNSTRETCH THEM TWICE...
				//scans[i].undoCosTransform(survey.getMedianDec(), survey.getMedianDec(), survey.getMedianRa());
			}
		}
		RCR rcr2 = RCR(LS_MODE_DL);
		rcr2.performBulkRejection(radii);
		survey.setEdgeRadius(rcr2.result.mu);
		if (survey.getMedianDec() + rcr2.result.mu > globalMaxDec)
		{
			globalMaxDec = survey.getMedianDec() + rcr2.result.mu;
		}

		if (survey.getMedianDec() - rcr2.result.mu < globalMinDec)
		{
			globalMinDec = survey.getMedianDec() - rcr2.result.mu;
		}

		if (survey.getMedianRa() + rcr2.result.mu > globalMaxRa)
		{
			globalMaxRa = survey.getMedianRa() + rcr2.result.mu;
		}

		if (survey.getMedianRa() - rcr2.result.mu < globalMinRa)
		{
			globalMinRa = survey.getMedianRa() - rcr2.result.mu;
		}

	}

	globalCenterDec = 0.5*(globalMaxDec - globalMinDec) + globalMinDec;
	globalCenterRa = 0.5*(globalMaxRa - globalMinRa) + globalMinRa;

	partSetProcSSS.centerDecDeg = globalCenterDec;
	partSetProcSSS.centerRaDeg = globalCenterRa;
	
	survey.setPartSetProcSSS(partSetProcSSS);//WE ONLY NEED SSS BECAUSE LSS CHANGES ONLY ONCE WE'VE REJECTED POINTS FROM SSS
}
void Processor::characterizeData(Survey &survey)
{
	scans = survey.getScans();
	mapType = survey.getMapType();

	if (mapType != DAISY)
	{
		scansInRa = survey.getScanDirection();
		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setScanInRa(scansInRa);
		}
	}
		
	for (int i = 0; i < scans.size(); i++)
	{
		//if (survey.getTracking())
		//{
		//	scans[i].cosDecTransform(survey.getTimeShift(), survey.getMedianDec(), survey.getMedianRa(), survey.getMedianDec());
		//	scans[i].updateAngDistTemp(survey.getMedianDec());
		//}
		//else
		//{
		//	scans[i].cosDecTransform(survey.getTimeShift(), globalCenterDec, globalCenterRa, globalCenterDec);
		//	scans[i].updateAngDistTemp(globalCenterDec);
		//}

		 if (survey.getTracking())
		 {
			if (survey.getProcessingCoordinate() == EQUATORIAL)
			{
				scans[i].dynamicCosDecTransform(survey.getTimeShift(), survey.getMedianDec(), survey.getMedianRa(), scans[i].getDec());
				scans[i].updateAngDistTemp(survey.getMedianDec());
			}
			else
			{
				scans[i].cosDecTransform(survey.getTimeShift(), survey.getMedianDec(), survey.getMedianRa(), survey.getMedianDec());
				scans[i].updateAngDistTemp(survey.getMedianDec());
			}
		}
		else
		{
			if (survey.getProcessingCoordinate() == EQUATORIAL)
			{
				scans[i].dynamicCosDecTransform(survey.getTimeShift(), globalCenterDec, globalCenterRa, scans[i].getDec());
				scans[i].updateAngDistTemp(globalCenterDec);
			}
			else
			{
				scans[i].cosDecTransform(survey.getTimeShift(), globalCenterDec, globalCenterRa, globalCenterDec);
				scans[i].updateAngDistTemp(globalCenterDec);
			}
		}
	}

	partSetProcSSS = survey.getPartSetProcSSS();

	if (survey.getTracking())
	{
		medianLatiMapped = partSetProcSSS.medianDec;
		medianLongMapped = partSetProcSSS.medianRa;
		partSetProcSSS.medianDec = 0;
		partSetProcSSS.medianRa = 0;
		trackingForEdges = true;
	}
	else
	{
		partSetProcSSS.medianDec = survey.getMedianDec() - globalCenterDec;
		partSetProcSSS.medianRa = (survey.getMedianRa() - globalCenterRa)*cos(globalCenterDec);
	}

	survey.setPartSetProcSSS(partSetProcSSS);


	switchChannels(survey.getChannel(), survey);
	calculateScatter();

	survey.setScans(scans);
	survey.setStandardThetaGap();

	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			scans[i].setMinRCRThetaGap(j, survey.getMinGapThreshold());
		}
	}
	survey.setScans(scans);

	Debugger::print("Info", "Scatters loaded");
}
void Processor::switchChannels(Channel newChannel, Survey &survey)
{
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].switchChannels(newChannel);
	}
}
void Processor::calculateScatter()
{
	double m, b, xBar;
	std::vector<double> indices, scatters, dumpSums, ret;

	for (int j = 0; j < scans.size(); j++)
	{
		scans[j].calculateScatter();
		scatters.push_back(scans[j].getScatter());
		dumpSums.push_back(scans[j].getDumpSum());
		indices.push_back(j);
	}

	RCR rcr = RCR(LS_MODE_68);
	LinearModel model = LinearModel(dumpSums, indices, scatters);
	rcr.setParametricModel(model);

	rcr.performBulkRejection(dumpSums, scatters);

	m = model.m;
	b = model.b;
	xBar = model.xBar;

	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setScatter(m * (i - xBar) + b);
		if (scans[i].getScatter() < 0.0)
		{
			scans[i].setScatter(scatters[i]);
		}
	}
}

//tools
int Processor::findPossSSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
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

				if(distance < limit)
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
int Processor::findPossLSS(double i_0, double j_0, double limit, std::vector<int> &toRet)
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
void Processor::classifySSS(Survey &survey)
{
	scans = survey.getScans();
	partSetProcSSS = survey.getPartSetProcSSS();
	//std::vector<std::vector<std::vector<int> > > classificationsSSS;
	classificationsSSS = survey.getClassificationsSSS();

	std::vector<int> intFiller;
	std::vector<std::vector<int> > intFillerFiller;
	std::vector<double> raHold, decHold;
	int decCounter, raCounter;
	double decTemp, raTemp;
	classificationsSSS.clear();

	for (int i = 0; i < partSetProcSSS.subDecRes; i++)
	{
		classificationsSSS.push_back(intFillerFiller);
		for (int j = 0; j < partSetProcSSS.subRaRes; j++)
		{
			classificationsSSS[i].push_back(intFiller);
		}
	}
	for (int i = 0; i < scans.size(); i++)
	{
		raHold = scans[i].getRa();
		decHold = scans[i].getDec();
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			raTemp = raHold[j];
			decTemp = decHold[j];
			decCounter = (int)((decTemp - partSetProcSSS.minDec) / partSetProcSSS.subDecInc);
			raCounter = (int)((raTemp - partSetProcSSS.minRa) / partSetProcSSS.subRaInc);
			if (decCounter <= 0)
			{
				decCounter = 0;
			}
			if (raCounter <= 0)
			{
				raCounter = 0;
			}
			if (decCounter >= classificationsSSS.size())
			{
				decCounter = classificationsSSS.size() - 1;
			}
			if (raCounter >= classificationsSSS[decCounter].size())
			{
				raCounter = classificationsSSS[decCounter].size() - 1;
			}
			classificationsSSS[decCounter][raCounter].push_back(i);
			classificationsSSS[decCounter][raCounter].push_back(j);
		}
	}
	survey.setClassificationsSSS(classificationsSSS);
}
void Processor::classifyLSS(Survey &survey)
{
	scans = survey.getScans();
	partSetProcLSS = survey.getPartSetProcLSS();
	classificationsLSS = survey.getClassificationsLSS();

	std::vector<int> intFiller;
	std::vector<std::vector<int> > intFillerFiller;
	std::vector<double> raHold, decHold;
	int decCounter, raCounter;
	double decTemp, raTemp;
	classificationsLSS.clear();

	for (int i = 0; i < partSetProcLSS.subDecRes; i++)
	{
		classificationsLSS.push_back(intFillerFiller);
		for (int j = 0; j < partSetProcLSS.subRaRes; j++)
		{
			classificationsLSS[i].push_back(intFiller);
		}
	}
	for (int i = 0; i < scans.size(); i++)
	{
		raHold = scans[i].getLSSRa();
		decHold = scans[i].getLSSDec();
		for (int j = 0; j < scans[i].getLSSSize(); j++)
		{
			raTemp = raHold[j];
			decTemp = decHold[j];
			decCounter = (int)((decTemp - partSetProcLSS.minDec) / partSetProcLSS.subDecInc);
			raCounter = (int)((raTemp - partSetProcLSS.minRa) / partSetProcLSS.subRaInc);
			if (decCounter <= 0)
			{
				decCounter = 0;
			}
			if (raCounter <= 0)
			{
				raCounter = 0;
			}
			if (decCounter >= classificationsLSS.size())
			{
				decCounter = classificationsLSS.size() - 1;
			}
			if (raCounter >= classificationsLSS[decCounter].size())
			{
				raCounter = classificationsLSS[decCounter].size() - 1;
			}
			classificationsLSS[decCounter][raCounter].push_back(i);
			classificationsLSS[decCounter][raCounter].push_back(j);
		}
	}
	//survey.setPartSetProcLSS(partSetProcLSS);
	survey.setClassificationsLSS(classificationsLSS);
}
void Processor::classifySSS(Composite & composite)
{
	std::vector<int> intFiller;
	std::vector<std::vector<int> > intFillerFiller;
	std::vector<double> raHold, decHold;
	int decCounter, raCounter;
	double decTemp, raTemp;
	//std::vector<std::vector<std::vector<int> > > classificationsSSS;

	scans = composite.getScans();
	partSetProcSSS = composite.getCompPartSetProcSSS();
	classificationsSSS = composite.getClassificationsSSS();
	classificationsSSS.clear();

	for (int i = 0; i < partSetProcSSS.subDecRes; i++)
	{
		classificationsSSS.push_back(intFillerFiller);
		for (int j = 0; j < partSetProcSSS.subRaRes; j++)
		{
			classificationsSSS[i].push_back(intFiller);
		}
	}
	for (int i = 0; i < scans.size(); i++)
	{
		raHold = scans[i].getRa();
		decHold = scans[i].getDec();
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			raTemp = raHold[j];
			decTemp = decHold[j];
			decCounter = (int)((decTemp - partSetProcSSS.minDec) / partSetProcSSS.subDecInc);
			raCounter = (int)((raTemp - partSetProcSSS.minRa) / partSetProcSSS.subRaInc);
			if (decCounter <= 0)
			{
				decCounter = 0;
			}
			if (raCounter <= 0)
			{
				raCounter = 0;
			}
			if (decCounter >= classificationsSSS.size())
			{
				decCounter = classificationsSSS.size() - 1;
			}
			if (raCounter >= classificationsSSS[decCounter].size())
			{
				raCounter = classificationsSSS[decCounter].size() - 1;
			}
			classificationsSSS[decCounter][raCounter].push_back(i);
			classificationsSSS[decCounter][raCounter].push_back(j);
		}
	}
	//composite.setCompPartSetProcSSS(classificationsSSS);
	composite.setClassificationsSSS(classificationsSSS);
}
void Processor::classifyLSS(Composite & composite)
{
	scans = composite.getScans();
	partSetProcLSS = composite.getCompPartSetProcLSS();
	classificationsLSS = composite.getClassificationsLSS();

	std::vector<int> intFiller;
	std::vector<std::vector<int> > intFillerFiller;
	std::vector<double> raHold, decHold;
	int decCounter, raCounter;
	double decTemp, raTemp;
	classificationsLSS.clear();

	for (int i = 0; i < partSetProcLSS.subDecRes; i++)
	{
		classificationsLSS.push_back(intFillerFiller);
		for (int j = 0; j < partSetProcLSS.subRaRes; j++)
		{
			classificationsLSS[i].push_back(intFiller);
		}
	}
	for (int i = 0; i < scans.size(); i++)
	{
		raHold = scans[i].getLSSRa();
		decHold = scans[i].getLSSDec();
		for (int j = 0; j < scans[i].getLSSSize(); j++)
		{
			raTemp = raHold[j];
			decTemp = decHold[j];
			decCounter = (int)((decTemp - partSetProcLSS.minDec) / partSetProcLSS.subDecInc);
			raCounter = (int)((raTemp - partSetProcLSS.minRa) / partSetProcLSS.subRaInc);
			if (decCounter <= 0)
			{
				decCounter = 0;
			}
			if (raCounter <= 0)
			{
				raCounter = 0;
			}
			if (decCounter >= classificationsLSS.size())
			{
				decCounter = classificationsLSS.size() - 1;
			}
			if (raCounter >= classificationsLSS[decCounter].size())
			{
				raCounter = classificationsLSS[decCounter].size() - 1;
			}
			classificationsLSS[decCounter][raCounter].push_back(i);
			classificationsLSS[decCounter][raCounter].push_back(j);
		}
	}
	//composite.setCompPartSetProcLSS(partSetProcLSS);
	composite.setClassificationsLSS(classificationsLSS);
}

PartitionSet Processor::determineProcSurveyDimensions(Survey &survey, bool LSS)
{
	//THIS FUNCTION DETERMINES SURVEY BOUNDARIES, BOTH PROCESSED AND RAW
	PartitionSet partSetProc;
	scans = survey.getScans();
	if (LSS)
	{
		partSetProc = survey.getPartSetProcLSS();
	}
	else
	{
		partSetProc = survey.getPartSetProcSSS();
	}
	
	RCR rcr = RCR(LS_MODE_DL);

	double decTemp, raTemp;
	partSetProc.minDec = 999999, partSetProc.minRa = 999999, partSetProc.maxDec = -999999, partSetProc.maxRa = -999999;
	partSetProc.mapType = survey.getMapType();

	int size;
	//USED TO FIND RCR MIN AND MAX DEC AND RA
	for (int i = 0; i < scans.size(); i++)
	{
		if (LSS)
		{
			size = scans[i].getLSSSize();
		}
		else
		{
			size = scans[i].getSize();
		}

		for (int j = 0; j < size; j++)
		{
			if (LSS)
			{
				raTemp = scans[i].getLSSRa(j);
				decTemp = scans[i].getLSSDec(j);
			}
			else
			{
				raTemp = scans[i].getRa(j);
				decTemp = scans[i].getDec(j);
			}


			if (raTemp > partSetProc.maxRa)
			{
				partSetProc.maxRa = raTemp;
			}
			if (raTemp < partSetProc.minRa)
			{
				partSetProc.minRa = raTemp;
			}
			if (decTemp > partSetProc.maxDec)
			{
				partSetProc.maxDec = decTemp;
			}
			if (decTemp < partSetProc.minDec)
			{
				partSetProc.minDec = decTemp;
			}
		}
	}

	if (partSetProc.mapType == DAISY)
	{
		RCR rcr = RCR(LS_MODE_DL);
		int centerIndex;
		double distance;
		std::vector<double> radiiHold;

		for (int i = 0; i < scans.size(); i++)
		{
			centerIndex = scans[i].getCenter();

			if (LSS)
			{
				size = scans[i].getLSSSize();
			}
			else
			{
				size = scans[i].getSize();
			}

			for (int j = 0; j < size; j++)
			{
				//distance = Tools::getGCDistance(scans[i].getDec(j), scans[i].getRa(j), partSetProc.medianDec, partSetProc.medianRa, partSetProc.centerDecDeg)*toDeg;//LAST TERM SHOULD PROBABLY BE 0
				//ONCE DAISY PATHS ARE CORRECTED THIS SHOULD NO LONGER BE MOD
				if (LSS)
				{
					distance = Tools::getModGCDistance(scans[i].getLSSDec(j), scans[i].getLSSRa(j), partSetProc.medianDec, partSetProc.medianRa)*toDeg;
				}
				else
				{
					distance = Tools::getModGCDistance(scans[i].getDec(j), scans[i].getRa(j), partSetProc.medianDec, partSetProc.medianRa)*toDeg;
				}
				
				radiiHold.push_back(distance);
			}
		}
		rcr.performBulkRejection(radiiHold);
		partSetProc.edgeRadius = rcr.result.mu; //DEGREES
		Debugger::print("Info", "Edge Radius in BW", partSetProc.edgeRadius / psfFWHM);
	}

	// If mapped in equatorial but want to process in galactic
	if (survey.getMappingCoordinate() == "equatorial" && survey.getProcessingCoordinate() == GALACTIC)
	{
		globalCenterLatProc  = Tools::convertToB(globalCenterRa, globalCenterDec);
		globalCenterLongProc = Tools::convertToL(globalCenterRa, globalCenterDec);
	}
	// If mapped in galactic but want to process in equatorial
	else if (survey.getMappingCoordinate() == "galactic" && survey.getProcessingCoordinate() == EQUATORIAL)
	{
		globalCenterLatProc  = Tools::convertToDec(globalCenterRa, globalCenterDec);
		globalCenterLongProc = Tools::convertToRa(globalCenterRa, globalCenterDec);
	}

	partSetProc.centerLatProcDeg = globalCenterLatProc;   // never used if coordinate == MAPPED.
	partSetProc.centerLongProcDeg = globalCenterLongProc;

	partSetProc.centerDecDeg = globalCenterDec;
	partSetProc.centerRaDeg = globalCenterRa;

	partSetProc.subDecRes = Tools::max(1, (int)((partSetProc.maxDec - partSetProc.minDec) / Tools::max(psfFWHM*rfiScaleBW, psfFWHM)));
	partSetProc.subRaRes = Tools::max(1, (int)(((partSetProc.maxRa - partSetProc.minRa) / Tools::max(psfFWHM*rfiScaleBW, psfFWHM))*Tools::min(cos(partSetProc.maxDec / toDeg), cos(partSetProc.minDec / toDeg))));

	partSetProc.subDecInc = (partSetProc.maxDec - partSetProc.minDec) / partSetProc.subDecRes;
	partSetProc.subRaInc = (partSetProc.maxRa - partSetProc.minRa) / partSetProc.subRaRes;

	partSetProc.medianLatiMap = medianLatiMapped;  // never used if tracking == false.
	partSetProc.medianLongMap = medianLongMapped;

	return partSetProc;
}
PartitionSet Processor::determineProcCompositeDimensions(Composite &composite, bool LSS)
{
	PartitionSet partSetProc;
	scans = composite.getScans();
	partSetProc = composite.getCompPartSetProcSSS();

	if (LSS)
	{
		partSetProc = composite.getCompPartSetProcLSS();
	}
	else
	{
		partSetProc = composite.getCompPartSetProcSSS();
	}

	double decTemp, raTemp;
	partSetProc.minDec = 999999, partSetProc.minRa = 999999, partSetProc.maxDec = -999999, partSetProc.maxRa = -999999;

	std::vector<double> raHold, decHold;

	int size;
	//USED TO FIND RCR MIN AND MAX DEC AND RA
	for (int i = 0; i < scans.size(); i++)
	{
		if (LSS)
		{
			size = scans[i].getLSSSize();
		}
		else
		{
			size = scans[i].getSize();
		}

		for (int j = 0; j < size; j++)
		{
			if (LSS)
			{
				raTemp = scans[i].getLSSRa(j);
				decTemp = scans[i].getLSSDec(j);
			}
			else
			{
				raTemp = scans[i].getRa(j);
				decTemp = scans[i].getDec(j);
			}


			if (raTemp > partSetProc.maxRa)
			{
				partSetProc.maxRa = raTemp;
			}
			if (raTemp < partSetProc.minRa)
			{
				partSetProc.minRa = raTemp;
			}
			if (decTemp > partSetProc.maxDec)
			{
				partSetProc.maxDec = decTemp;
			}
			if (decTemp < partSetProc.minDec)
			{
				partSetProc.minDec = decTemp;
			}
		}
	}

	if (composite.getMappingCoordinate() == "equatorial" && composite.getProcessingCoordinate() == GALACTIC)
	{
		globalCenterLatProc  = Tools::convertToB(globalCenterRa, globalCenterDec);
		globalCenterLongProc = Tools::convertToL(globalCenterRa, globalCenterDec);
	}
	else if (composite.getMappingCoordinate() == "galactic" && composite.getProcessingCoordinate() == EQUATORIAL)
	{
		globalCenterLatProc = Tools::convertToDec(globalCenterRa, globalCenterDec);
		globalCenterLongProc = Tools::convertToRa(globalCenterRa, globalCenterDec);
	}

	partSetProc.centerLatProcDeg = globalCenterLatProc;    // never used if coordinate == MAPPED.
	partSetProc.centerLongProcDeg = globalCenterLongProc;

	partSetProc.centerDecDeg = globalCenterDec;
	partSetProc.centerRaDeg = globalCenterRa;

	partSetProc.tracking = trackingForEdges;

	partSetProc.subDecRes = Tools::max(1, (int)((partSetProc.maxDec - partSetProc.minDec) / Tools::max(psfFWHM*rfiScaleBW, psfFWHM)));
	partSetProc.subRaRes = Tools::max(1, (int)(((partSetProc.maxRa - partSetProc.minRa) / Tools::max(psfFWHM*rfiScaleBW, psfFWHM))*Tools::min(cos(partSetProc.maxDec / toDeg), cos(partSetProc.minDec / toDeg))));

	partSetProc.subDecInc = (partSetProc.maxDec - partSetProc.minDec) / partSetProc.subDecRes;
	partSetProc.subRaInc = (partSetProc.maxRa - partSetProc.minRa) / partSetProc.subRaRes;

	partSetProc.medianLatiMap = medianLatiMapped;   // never used if tracking == false.
	partSetProc.medianLongMap = medianLongMapped;

	//composite.setCompPartSetProcSSS(partSetProc);
	//composite.setCompPartSetProcLSS(partSetProc);
	return partSetProc;
}

void Processor::forceThetaGap(Composite & composite, double value)
{
	scans = composite.getScans();
	double forcedVal = value * 0.75*psfFWHM*toRad;
	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			scans[i].setThetaGap(j, forcedVal);
		}
	}
	composite.setScans(scans);
}
void Processor::repartitionScans(Composite & composite, std::vector<Survey> & surveys)
{
	std::vector<Scan> scansHold;
	PartitionSet partSet;
	for (int i = 0; i < surveys.size(); i++)
	{
		for (int j = 0; j < scans.size(); j++)
		{
			if (surveys[i].getSurveyNumber() == scans[j].getSurveyNumber())
			{
				scansHold.push_back(scans[j]);
			}
		}
		surveys[i].setScans(scansHold);
		
		partSet = determineProcSurveyDimensions(surveys[i], false);
		surveys[i].setPartSetProcSSS(partSet);
		partSet = determineProcSurveyDimensions(surveys[i], true);
		surveys[i].setPartSetProcLSS(partSet);

		classifySSS(surveys[i]);
		classifyLSS(surveys[i]);

		surveys[i].setCentroidLocations(centroidLocations);
		scansHold.clear();
		scans = composite.getScans();
	}
}


//SMALL SCALE STRUCTURE PROCESSING

//shift
std::vector<bool> Processor::shiftRejection(std::vector<double> shifts)
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
double Processor::calculateShift(std::vector<std::vector<double> > &maxIfftInfo)
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

		Debugger::print("Info", rcr.result.mu, rcr.result.sigma);

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
std::vector<std::vector<double>> Processor::getScanToScanShifts(double sampling)
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
			angleVec1 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg,  scans[i].getRa(), scans[i].getDec());
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
		Debugger::print("Info", i);
	}

	return maxIfftInfo;
}



//LARGE SCALE STRUCTURE PROCESSING
void Processor::lssDropCorrection(Survey &survey)
{

	//stop pushing back the pre and post indices and instead just use them to call the continued iterations of 2D Scatter
	RCR rcr = RCR(SS_MEDIAN_DL);
	scans = survey.getScans();

	int indexTemp;
	double max = 0, percentRejected = 0, boostValue = 0;

	std::vector<int> preIndices, proIndices;
	std::vector<double> delta, maxDetails;
	
	while (true)
	{
		delta.clear();

		preIndices.clear();
		proIndices.clear();

		delta = setDropDeltas();
		setDropFlags(delta);

		std::vector<double> holdAng1, holdAng2, holdAng3;
		//COLLECT THE PRE AND POST INDICES FOR EACH POINT TO BE TESTED
		for (int i = 0; i < scans.size(); i++)
		{
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				if (scans[i].getDropFlag(j) == 1 || scans[i].getDropFlag(j) == -1)
				{
					if (scansInRa)
					{
						holdAng1 = scans[i].getRa();
						if (i != 0)
						{
							holdAng2 = scans[i - 1].getRa();
						}
						if (i != scans.size() - 1)
						{
							holdAng3 = scans[i + 1].getRa();
						}
					}
					else
					{
						holdAng1 = scans[i].getDec();
						if (i != 0)
						{
							holdAng2 = scans[i - 1].getDec();
						}
						if (i != scans.size() - 1)
						{
							holdAng3 = scans[i + 1].getDec();
						}
					}

					if (i == 0)
					{
						preIndices.push_back(-1);
						//proIndices.push_back(Tools::determineNearestIndex(j, holdAng1, holdAng3));

						indexTemp = Tools::determineNearestIndex(j, holdAng1, holdAng3);
						//indexTemp = correctIndex(j, indexTemp, i, i + 1);
						indexTemp = correctIndex(i, j, i + 1, indexTemp);// j, indexTemp, i, i + 1);
						proIndices.push_back(indexTemp);
					}
					else if (i == scans.size() - 1)
					{
						indexTemp = Tools::determineNearestIndex(j, holdAng1, holdAng2);
						indexTemp = correctIndex(i, j, i - 1, indexTemp);//correctIndex(j, indexTemp, i, i - 1);
						preIndices.push_back(indexTemp);

						//preIndices.push_back(Tools::determineNearestIndex(j, holdAng1, holdAng2)); 
						proIndices.push_back(999999);
					}
					else
					{
						indexTemp = Tools::determineNearestIndex(j, holdAng1, holdAng2);
						indexTemp = correctIndex(i, j, i - 1, indexTemp);//correctIndex(j, indexTemp, i, i - 1);
						preIndices.push_back(indexTemp); // CHANGE TO BE DEC OR RA

						indexTemp = Tools::determineNearestIndex(j, holdAng1, holdAng3);
						indexTemp = correctIndex(i, j, i + 1, indexTemp);// correctIndex(j, indexTemp, i, i + 1);
						proIndices.push_back(indexTemp);

						//preIndices.push_back(Tools::determineNearestIndex(j, holdAng1, holdAng2)); // CHANGE TO BE DEC OR RA
						//proIndices.push_back(Tools::determineNearestIndex(j, holdAng1, holdAng3));
					}
				}
			}
		}

		maxDetails = determinePercentDrop(preIndices, proIndices);

		if (max < 0.5)
		{
			Debugger::print("Info", "No Drops");
			return;
		}
		else
		{
			correctDrop(maxDetails);
		}

		survey.setScans(scans);
	}
}
void Processor::lssThetaGapCalc(Survey &survey)
{
	double value;
	std::vector<std::vector<double>> thetaGapLSS;
	scans = survey.getScans();
	ProcessorThetaGap procTG = ProcessorThetaGap(scans, psfFWHM, partSetProcSSS, partSetProcLSS, classificationsSSS, classificationsLSS);
	thetaGapLSS = procTG.calculateThetaGapLSS();
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setLSSThetaGap(thetaGapLSS[i]);
	}
	survey.setScans(scans);
}
void Processor::lssSet2DScatter(Survey &survey)
{
	scans = survey.getScans();
	RCR rcr = RCR(LS_MODE_DL);
	int jBefore = 0, jAfter = 0;
	std::vector<double> regRet;
	double scatter2DSum, iterCount = 0, scatter2DAve;
	double minAngle1, maxAngle1, minAngle2, maxAngle2, minAngle3, maxAngle3, angleHold, minAngle, maxAngle, toPush, weight, weightLow, weightHigh, xBar, angleHigh, angleLow, angleMid, weightMid;
	std::vector<double> angleVec1, angleVec2, angleVec3, delta, weights, scatter2DAverages;
	std::vector<double> scatter2dTemp;
	scatter2dTemp.reserve(scans.size());
	std::ofstream scatter2DFile;
	scatter2DSum = 0;

	for (int i = 1; i < scans.size() - 1; i++)
	{
		minAngle1 = 999999;
		maxAngle1 = -999999;
		minAngle2 = 999999;
		maxAngle2 = -999999;
		minAngle3 = 999999;
		maxAngle3 = -999999;
		if (mapType == DAISY)
		{
			angleVec1 = Tools::daisyAngleBuilder(i - 1, scans[i - 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i - 1].getRa(), scans[i - 1].getDec());
			angleVec2 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg, scans[i].getRa(), scans[i].getDec());
			angleVec3 = Tools::daisyAngleBuilder(i + 1, scans[i + 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i + 1].getRa(), scans[i + 1].getDec());
		}
		else
		{
			if (scansInRa)
			{
				angleVec1 = scans[i - 1].getRa();
				angleVec2 = scans[i].getRa();
				angleVec3 = scans[i + 1].getRa();
			}
			else
			{
				angleVec1 = scans[i - 1].getDec();
				angleVec2 = scans[i].getDec();
				angleVec3 = scans[i + 1].getDec();
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
		for (int j = 0; j < angleVec3.size(); j++)
		{
			angleHold = angleVec3[j];
			if (angleHold < minAngle3)
			{
				minAngle3 = angleHold;
			}
			if (angleHold > maxAngle3)
			{
				maxAngle3 = angleHold;
			}
		}
		minAngle = minAngle1;
		if (minAngle2 > minAngle)
		{
			minAngle = minAngle2;
		}
		if (minAngle3 > minAngle)
		{
			minAngle = minAngle3;
		}
		maxAngle = maxAngle1;
		if (maxAngle2 < maxAngle)
		{
			maxAngle = maxAngle2;
		}
		if (maxAngle3 < maxAngle)
		{
			maxAngle = maxAngle3;
		}
		if (angleVec2[angleVec2.size() - 1] > angleVec2[0])
		{
			int j = 0, k = 0;
			jBefore = angleVec1.size() - 1, jAfter = angleVec3.size() - 1;
			while (minAngle > angleVec2[j])
			{
				j++;
			}
			while (maxAngle > angleVec2[k])
			{
				k++;
			}
			for (; j < k; j++)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore--;
				}
				if (jBefore < angleVec1.size() - 1 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore + 1] - angleVec2[j]))
				{
					jBefore++;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter--;
				}
				if (jAfter < angleVec3.size() - 1 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter + 1] - angleVec2[j]))
				{
					jAfter++;
				}
				if (scansInRa)
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getDec(jAfter);
					angleMid = scans[i].getDec(j);
					angleLow = scans[i - 1].getDec(jBefore);
				}
				else
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getRa(jAfter);
					angleMid = scans[i].getRa(j);
					angleLow = scans[i - 1].getRa(jBefore);
				}

				weightLow = scans[i - 1].getDataDumps(jBefore);
				weightHigh = scans[i + 1].getDataDumps(jAfter);
				xBar = (weightLow * angleLow + weightHigh * angleHigh) / (weightLow + weightHigh);
				weight = (1.0 / scans[i].getDataDumps(j));
				weight += pow(((angleHigh - xBar) * (angleHigh - angleLow)), 2.0) / weightLow;
				weight += pow(((angleLow - xBar) * (angleHigh - angleLow)), 2.0) / weightHigh;
				weight += pow(((angleMid - xBar) * (angleHigh - angleLow)), 2.0) * ((1.0 / weightHigh) + (1.0 / weightLow));
				weights.push_back(1.0 / weight);

			}
			rcr.performBulkRejection(delta);
			scatter2dTemp.push_back(std::sqrt(rcr.result.mu * rcr.result.mu + rcr.result.sigma * rcr.result.sigma));
			delta.clear();
			weights.clear();
		}
		else
		{
			int j = angleVec2.size() - 1, k = angleVec2.size() - 1;
			jBefore = 0, jAfter = 0;
			while (maxAngle > angleVec2[k])
			{
				k--;
			}
			while (minAngle > angleVec2[j])
			{
				j--;
			}
			for (; j > k; j--)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore++;
				}
				if (jBefore > 0 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore - 1] - angleVec2[j]))
				{
					jBefore--;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter++;
				}
				if (jAfter > 0 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter - 1] - angleVec2[j]))
				{
					jAfter--;
				}
				if (scansInRa)
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getDec(jAfter);
					angleMid = scans[i].getDec(j);
					angleLow = scans[i - 1].getDec(jBefore);
				}
				else
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					angleHigh = scans[i + 1].getRa(jAfter);
					angleMid = scans[i].getRa(j);
					angleLow = scans[i - 1].getRa(jBefore);
				}

				weightLow = scans[i - 1].getDataDumps(jBefore);
				weightHigh = scans[i + 1].getDataDumps(jAfter);
				xBar = (weightLow * angleLow + weightHigh * angleHigh) / (weightLow + weightHigh);
				weight = (1.0 / scans[i].getDataDumps(j));
				weight += pow(((angleHigh - xBar) * (angleHigh - angleLow)), 2.0) / weightLow;
				weight += pow(((angleLow - xBar) * (angleHigh - angleLow)), 2.0) / weightHigh;
				weight += pow(((angleMid - xBar) * (angleHigh - angleLow)), 2.0) * ((1.0 / weightHigh) + (1.0 / weightLow));
				weights.push_back(1.0 / weight);

			}
			rcr.performBulkRejection(delta);
			scatter2dTemp.push_back(std::sqrt(rcr.result.mu * rcr.result.mu + rcr.result.sigma * rcr.result.sigma));
			delta.clear();
			weights.clear();
		}
	}

	std::vector<double> indices, dumpSums;

	for (int i = 1; i < scans.size() - 1; i++)
	{
		dumpSums.push_back(scans[i].getDumpSum());
		indices.push_back(i);

	}
	for (int i = 0; i < scatter2dTemp.size(); i++)
	{
		scatter2dTemp[i] *= .8137;
	}
	RCR rcr2 = RCR(LS_MODE_68);
	LinearModel model = LinearModel(dumpSums, indices, scatter2dTemp);
	rcr2.setParametricModel(model);
	rcr2.performBulkRejection(dumpSums, scatter2dTemp);

	bool useAverage = false;
	double average2dScatter;

	if (scatter2dTemp.size() == 2)
	{
		useAverage = true;
		average2dScatter = (scatter2dTemp[0] + scatter2dTemp[1]) / 2.0;
	}

	double m, b;
	m = model.m;
	b = model.b;
	xBar = model.xBar;
	scatter2dTemp.clear();
	scatter2dTemp.resize(scans.size(), 0.0);

	for (int i = 0; i < scans.size(); i++)
	{
		if (useAverage)
		{
			scatter2dTemp[i] = average2dScatter;
		}
		else
		{
			scatter2dTemp[i] = m*(i - xBar) + b;
		}
		scatter2dTemp[i] = std::max(scans[i].getScatter(), scatter2dTemp[i]);
	}

	survey.setLSS2DScatterVec(scatter2dTemp);

	std::ofstream scatterFile;
	scatterFile.open("LSS_Scatter.txt");
	for (int i = 0; i < scans.size(); i++)
	{
		scatterFile << scatter2dTemp[i] << "\n";
	}
	scatterFile.close();
}
void Processor::lssPolynomialRFISubtraction(Survey &survey, double bgScale, double wScale)
{
	bool multi = true;
	std::vector<double> dubFiller;
	std::vector<std::vector<double> > cleanedLSSGrid;

	for (int i = 0; i < scans.size(); i++)
	{
		dubFiller.resize(scans[i].getSize(), 0.0);
		cleanedLSSGrid.push_back(dubFiller);
	}
	
	//cleanedLSSGrid[58][9] = lssRecursivePolyFit(survey, 120, 40, bgScale, wScale);
	//cleanedLSSGrid[58][9] = lssRecursivePolyFit(survey, 59, 10, bgScale, wScale);
	
	if (multi)
	{
		std::vector<std::future<double>> futureVec;
		for (int i = 0; i < scans.size(); i++)
		{
			futureVec.resize(scans[i].getSize());
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				futureVec[j] = std::async(std::launch::async, &Processor::lssRecursivePolyFit, this, std::ref(survey), i, j, bgScale, wScale);
			}
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				cleanedLSSGrid[i][j] = futureVec[j].get();
			}
		}
	}
	else
	{
		for (int i = 0; i < scans.size(); i++)
		{
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				cleanedLSSGrid[i][j] = lssRecursivePolyFit(survey, i, j, bgScale, wScale);
			}
		}
	}


	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setLSSData(cleanedLSSGrid[i]);
	}
	survey.setScans(scans);
}
void Processor::lssElevationSubtraction(Survey &survey)
{
	int counter = 0;
	int liveThreads = 0;
	int maxThreads = 8;
	int completedThreads = 0;
	std::vector<double> results;
	BackgroundCUDA bgCuda;
	std::vector<std::future<std::vector<double>>> futureVec;
	futureVec.resize(scans.size());
	scans = survey.getScans();
	calculateSurveyElevationScatter(); //LOOKS LIKE THIS IS 2D BACKGROUND AGAIN INSTEAD OF THIS CALCULATE SCATTER V2

	for (int i = 0; i < scans.size(); i++)
	{
		Debugger::print("Info", i);
		bgCuda = BackgroundCUDA(scans[i], true);
		futureVec[i] = std::async(std::launch::async, &BackgroundCUDA::calculateBGMulti, bgCuda, bgScaleBW*psfFWHM);
		counter++;
		liveThreads++;

		if (liveThreads >= maxThreads)
		{
			for (int i = completedThreads; i < counter; i++)
			{
				results = futureVec[i].get();
				scans[i].removeBG(results);
			}
			completedThreads += liveThreads;
			liveThreads = 0;
		}
	}
	for (int i = completedThreads; i < scans.size(); i++)
	{
		results = futureVec[i].get();
		scans[i].removeElevation(results);
	}

	survey.setScans(scans);
}


//lss drop
std::vector<double> Processor::determinePercentDrop(std::vector<int> & preIndices, std::vector<int> & proIndices)
{
	int indexCounter = 0, indexCounter2 = 0;
	int rejectedPoints = 0, totalPoints = 0;
	int maxScan = 0, maxData = 0, maxPre = 0, maxPro = 0;
	int maxProTemp;
	int postsPre = -1, presPost = -1;
	int indexTemp;
	double max = 0, percentRejected = 0, boostValue = 0;
	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			if (scans[i].getDropFlag(j) == 1 || scans[i].getDropFlag(j) == -1)
			{
				rejectedPoints = 0;
				totalPoints = 0;

				if (i == scans.size() - 1)
				{
					maxProTemp = 0;
				}
				else
				{
					maxProTemp = scans[i + 1].getSize();
				}

				if (i != 0)
				{
					if (preIndices[indexCounter2] == scans[i - 1].getSize())
					{
						preIndices[indexCounter2] = scans[i - 1].getSize() - 1;
					}
				}

				if (proIndices[indexCounter2] == -1)
				{
					proIndices[indexCounter2] = 0;
				}

				//BEFORE AND AFTER POINT IN SAME SCAN
				for (int k = 0; k < scans[i].getSize(); k++)
				{
					if (scans[i].getDropFlag(k) == 1)
					{
						rejectedPoints++;
						totalPoints++;
					}
					else if (scans[i].getDropFlag(k) == 0)
					{
						totalPoints++;
					}
				}

				//PREVIOUS SCAN UP UNTIL PRE INDEX
				for (int k = 0; k < preIndices[indexCounter2]; k++)
				{
					if (scans[i - 1].getDropFlag(k) == 1)
					{
						rejectedPoints++;
						totalPoints++;
					}
					else if (scans[i - 1].getDropFlag(k) == 0)
					{
						totalPoints++;
					}
				}


				//PROCEEDING SCAN UP UNTIL PRO INDEX
				for (int k = proIndices[indexCounter2]; k < maxProTemp; k++) //WE NEED A maxProTemp BECAUSE WE CAN'T CALL scans[i+1] SO WE JUST PROVIDE A FILLER
				{

					if (scans[i + 1].getDropFlag(k) == 1)
					{
						rejectedPoints++;
						totalPoints++;
					}
					else if (scans[i + 1].getDropFlag(k) == 0)
					{
						totalPoints++;
					}
				}




				/*
				if (i == 0 || i == scans.size() - 1)
				{
				totalPoints++;
				}
				*/

				if (totalPoints == 0)
				{
					percentRejected = 0.0;
				}
				else
				{
					percentRejected = (double)rejectedPoints / (double)totalPoints;
				}

				if (percentRejected > max)// && (i != 0 && i != scans.size()-1))// did not solve the problem
				{
					max = percentRejected;
					//std::cout << "max\t" << max << "\n";
					maxScan = i;
					maxData = j;
					maxPre = preIndices[indexCounter2];
					maxPro = proIndices[indexCounter2];
				}

				indexCounter2++;
			}
		}
	}
	std::vector<double> maxDetails;
	maxDetails.resize(5);
	maxDetails[0] = max;
	maxDetails[1] = maxScan;
	maxDetails[2] = maxData;
	maxDetails[3] = maxPro;
	maxDetails[4] = maxPre;

	return maxDetails;
}
std::vector<double> Processor::setDropDeltas()
{
	int jBefore = 0, jAfter = 0;
	double minAngle1, maxAngle1, minAngle2, maxAngle2, minAngle3, maxAngle3, angleHold, minAngle, maxAngle, toPush, weight, weightLow, weightHigh, xBar, angleHigh, angleLow, angleMid, weightMid;
	std::vector<double> angleVec1, angleVec2, angleVec3, delta, scatter2DAverages;
	for (int i = 1; i < scans.size() - 1; i++)
	{
		minAngle1 = 999999;
		maxAngle1 = -999999;
		minAngle2 = 999999;
		maxAngle2 = -999999;
		minAngle3 = 999999;
		maxAngle3 = -999999;
		if (mapType == DAISY)
		{
			angleVec1 = Tools::daisyAngleBuilder(i - 1, scans[i - 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i - 1].getRa(), scans[i - 1].getDec());
			angleVec2 = Tools::daisyAngleBuilder(i, scans[i].getCenter(), partSetProcSSS.centerDecDeg, scans[i].getRa(), scans[i].getDec());
			angleVec3 = Tools::daisyAngleBuilder(i + 1, scans[i + 1].getCenter(), partSetProcSSS.centerDecDeg, scans[i + 1].getRa(), scans[i + 1].getDec());
		}
		else
		{
			if (scansInRa)
			{
				angleVec1 = scans[i - 1].getRa();
				angleVec2 = scans[i].getRa();
				angleVec3 = scans[i + 1].getRa();
			}
			else
			{
				angleVec1 = scans[i - 1].getDec();
				angleVec2 = scans[i].getDec();
				angleVec3 = scans[i + 1].getDec();
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
		for (int j = 0; j < angleVec3.size(); j++)
		{
			angleHold = angleVec3[j];
			if (angleHold < minAngle3)
			{
				minAngle3 = angleHold;
			}
			if (angleHold > maxAngle3)
			{
				maxAngle3 = angleHold;
			}
		}
		minAngle = minAngle1;
		if (minAngle2 > minAngle)
		{
			minAngle = minAngle2;
		}
		if (minAngle3 > minAngle)
		{
			minAngle = minAngle3;
		}
		maxAngle = maxAngle1;
		if (maxAngle2 < maxAngle)
		{
			maxAngle = maxAngle2;
		}
		if (maxAngle3 < maxAngle)
		{
			maxAngle = maxAngle3;
		}
		if (angleVec2[angleVec2.size() - 1] > angleVec2[0])
		{
			int j = 0, k = 0;
			jBefore = angleVec1.size() - 1, jAfter = angleVec3.size() - 1;
			while (minAngle > angleVec2[j])
			{
				j++;
			}
			while (maxAngle > angleVec2[k])
			{
				k++;
			}
			for (; j < k; j++)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore--;
				}
				if (jBefore < angleVec1.size() - 1 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore + 1] - angleVec2[j]))
				{
					jBefore++;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter--;
				}
				if (jAfter < angleVec3.size() - 1 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter + 1] - angleVec2[j]))
				{
					jAfter++;
				}
				if (scansInRa)
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					scans[i].setDrop2dValues(j, toPush);

				}
				else
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					scans[i].setDrop2dValues(j, toPush);
				}
			}
		}
		else
		{
			int j = angleVec2.size() - 1, k = angleVec2.size() - 1;
			jBefore = 0, jAfter = 0;
			while (maxAngle > angleVec2[k])
			{
				k--;
			}
			while (minAngle > angleVec2[j])
			{
				j--;
			}
			for (; j > k; j--)
			{
				while (angleVec1[jBefore] < angleVec2[j])
				{
					jBefore++;
				}
				if (jBefore > 0 && std::abs(angleVec1[jBefore] - angleVec2[j]) > std::abs(angleVec1[jBefore - 1] - angleVec2[j]))
				{
					jBefore--;
				}
				while (angleVec3[jAfter] < angleVec2[j])
				{
					jAfter++;

				}
				if (jAfter > 0 && std::abs(angleVec3[jAfter] - angleVec2[j]) > std::abs(angleVec3[jAfter - 1] - angleVec2[j]))
				{
					jAfter--;
				}
				if (scansInRa)
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getDec(jAfter) - scans[i - 1].getDec(jBefore))*(scans[i].getDec(j) - scans[i - 1].getDec(jBefore)));
					delta.push_back(toPush);
					scans[i].setDrop2dValues(j, toPush);
				}
				else
				{
					toPush = scans[i].getLSSData(j) - (scans[i - 1].getLSSData(jBefore) + (scans[i + 1].getLSSData(jAfter) - scans[i - 1].getLSSData(jBefore)) / (scans[i + 1].getRa(jAfter) - scans[i - 1].getRa(jBefore))*(scans[i].getRa(j) - scans[i - 1].getRa(jBefore)));
					delta.push_back(toPush);
					scans[i].setDrop2dValues(j, toPush);
				}
			}
		}
	}
	return delta;
}
void Processor::setDropFlags(std::vector<double> & delta)
{
	std::vector<double> nonRejected;
	//FIND THE AVERAGE DEVIATION FROM NOISE LEVEL AND THE LARGEST AND SMALLEST NON-REJECTED POINTS
	RCR rcr = RCR(SS_MEDIAN_DL);
	rcr.performBulkRejection(delta);
	std::vector<bool> flags = rcr.result.flags;
	nonRejected = rcr.result.cleanY;
	double maxVal = -999999;// Tools::max(nonRejected);
	double minVal = 999999;// Tools::min(nonRejected);

	for (int i = 0; i < nonRejected.size(); i++)
	{
		if (minVal > nonRejected[i])
		{
			minVal = nonRejected[i];
		}
		if (maxVal < nonRejected[i])
		{
			maxVal = nonRejected[i];
		}
	}

	//IDENTIFY REJECTED POINTS (1), THOSE WITH DIFFERENCES BUT NOT REJECTED (0), AND THOSE THAT HAVE NO DIFFERENCES (-1)
	for (int i = 0; i < scans.size(); i++)
	{
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			if (scans[i].getDrop2dValues(j) != -999999)
			{
				if ((scans[i].getDrop2dValues(j) > maxVal) || (scans[i].getDrop2dValues(j) < minVal))
				{
					scans[i].setDropFlag(j, 1);
				}
				else
				{
					scans[i].setDropFlag(j, 0);
				}
			}
		}
	}
}
void Processor::correctDrop(std::vector<double> maxDetails)
{
	double max = maxDetails[0];
	int maxScan = (int)maxDetails[1];
	int maxData = (int)maxDetails[2];
	int maxPro = (int)maxDetails[3];
	int maxPre = (int)maxDetails[4];
	std::vector<double> angTemp1, angTemp2, data;
	std::vector<bool> flagHolder;
	double boostValue;
	int maxProTemp, postsPre, presPost;

	if (maxScan == scans.size() - 1)
	{
		maxProTemp = 0;
	}
	else
	{
		maxProTemp = scans[maxScan + 1].getSize();
	}

	if (maxScan == 0)
	{
		if (scansInRa)
		{
			angTemp1 = scans[maxScan + 1].getRa();
			angTemp2 = scans[maxScan].getRa();
		}
		else
		{
			angTemp1 = scans[maxScan + 1].getDec();
			angTemp2 = scans[maxScan].getDec();
		}

		postsPre = Tools::determineNearestIndex(maxPro, angTemp1, angTemp2);
		maxData = postsPre;
	}
	else if (maxScan == scans.size() - 1)
	{
		if (scansInRa)
		{
			angTemp1 = scans[maxScan - 1].getRa();
			angTemp2 = scans[maxScan].getRa();
		}
		else
		{
			angTemp1 = scans[maxScan - 1].getDec();
			angTemp2 = scans[maxScan].getDec();
		}

		presPost = Tools::determineNearestIndex(maxPro, angTemp1, angTemp2);
		maxData = presPost;
	}

	for (int j = 0; j <= maxData; j++)
	{
		if (scans[maxScan].getDrop2dValues(j) != -999999)
		{
			data.push_back(-1 * scans[maxScan].getDrop2dValues(j));
		}
	}
	for (int j = maxData + 1; j < scans[maxScan].getSize(); j++)
	{
		if (scans[maxScan].getDrop2dValues(j) != -999999)
		{
			data.push_back(scans[maxScan].getDrop2dValues(j));
		}
	}
	for (int j = 0; j < maxPre; j++)
	{
		if (scans[maxScan - 1].getDrop2dValues(j) != -999999)
		{
			data.push_back(-1 * scans[maxScan - 1].getDrop2dValues(j));
		}
	}
	for (int j = maxPro; j < maxProTemp; j++)
	{
		if (scans[maxScan + 1].getDrop2dValues(j) != -999999)
		{
			data.push_back(scans[maxScan + 1].getDrop2dValues(j));
		}
	}

	RCR rcr2 = RCR(LS_MODE_DL);
	rcr2.performBulkRejection(data);

	std::ofstream outFile;
	outFile.open("Rejections.txt");
	outFile << "Indices:\t" << maxScan << "\t" << maxData << "\t" << scans[maxScan].getSize() << "\n";
	//outFile << "Boost Value:\t" << boostValue << "\n";
	flagHolder = rcr2.result.flags;
	for (int q = 0; q < data.size(); q++)
	{
		outFile << data[q] << "\t" << flagHolder[q] << "\n";
	}
	outFile.close();


	boostValue = 2 * rcr2.result.mu;

	//std::cout << "Indices:\t" << maxScan << "\t" << maxData << "\t" << scans[maxScan].getSize() << "\n";
	Debugger::print("Info", "Boost Value", boostValue);

	for (int i = maxScan; i < scans.size(); i++)
	{
		if (i == maxScan)
		{
			for (int j = maxData + 1; j < scans[i].getSize(); j++)
			{
				scans[i].setLSSData(j, scans[i].getLSSData(j) - boostValue);
			}
		}
		else
		{
			for (int j = 0; j < scans[i].getSize(); j++)
			{
				scans[i].setLSSData(j, scans[i].getLSSData(j) - boostValue);
			}
		}
	}

}
int Processor::correctIndex(int i1, int j1, int i2, int j2)//int j1, int j2, int i1, int i2)
{
	bool incWithIndex;
	bool proceeding = false;
	bool preceeding = false;
	bool j2Greater;
	std::vector<double> vec1, vec2;

	if (scansInRa)
	{
		vec1 = scans[i1].getRa();
		vec2 = scans[i2].getRa();
	}
	else
	{
		vec1 = scans[i1].getDec();
		vec2 = scans[i2].getDec();
	}

	//ENSURE STANDARD CASE
	if (j2 == -1)
	{
		return -1;
	}

	if (j2 == vec2.size())
	{
		return vec2.size();
	}

	//SCAN BEFORE OR AFTER
	if (i2 > i1)
	{
		proceeding = true;
	}
	else
	{
		proceeding = false;
	}

	if (vec1[vec1.size() - 1] > vec1[0])
	{
		incWithIndex = true;
	}
	else
	{
		incWithIndex = false;
	}


	while (true)
	{
		if (vec2[j2] > vec1[j1])
		{
			j2Greater = true;
		}
		else
		{
			j2Greater = false;
		}

		if (proceeding && incWithIndex && j2Greater)
		{
			j2++;
		}
		else if (proceeding && !incWithIndex && !j2Greater)
		{
			j2++;
		}
		else if (!proceeding && incWithIndex && !j2Greater)
		{
			j2--;
		}
		else if (!proceeding && !incWithIndex && j2Greater)
		{
			j2--;
		}
		else
		{
			return j2;
		}
	}

	/*
	if (i2 > i1)
	{
	proceeding = true;
	}
	else
	{
	preceeding = true;
	}

	if (vec1[vec1.size() - 1] > vec1[0])
	{
	posScanDirection = true;
	}
	else
	{
	posScanDirection = false;
	}


	if (posScanDirection)
	{
	if (proceeding)
	{
	while (vec2[j2] < vec1[j1])
	{
	j2++;
	}
	}

	if (preceeding)
	{
	while (vec2[j2] > vec1[j1])
	{
	j2--;
	}
	}
	}
	else
	{

	if (proceeding)
	{
	while (vec2[j2] > vec1[j1])
	{
	j2++;
	}
	}

	if (preceeding)
	{
	while (vec2[j2] < vec1[j1])
	{
	j2--;
	}
	}
	}
	*/
	return j2;
}

//elevation
void Processor::calculateSurveyElevationScatter()
{
	double m, b, xBar;
	std::vector<double> indices, scatters, dumpSums, ret;

	for (int j = 0; j < scans.size(); j++)
	{
		scans[j].calculateElevationScatter();
		scatters.push_back(scans[j].getElevationScatter());
		dumpSums.push_back(scans[j].getDumpSum());
		indices.push_back(j);
	}

	RCR rcr = RCR(LS_MODE_68);
	LinearModel model = LinearModel(dumpSums, indices, scatters);
	rcr.setParametricModel(model);

	rcr.performBulkRejection(dumpSums, scatters);

	m = model.m;
	b = model.b;
	xBar = model.xBar;

	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setElevationScatter(m * (i - xBar) + b);
		if (scans[i].getScatter() < 0.0)
		{
			scans[i].setElevationScatter(scatters[i]);
		}
	}
}
void Processor::calculateElevationBG(int scanIndex, double baseline)
{
	Debugger::print("Info", scanIndex);
	int size = scans[scanIndex].getSize();
	double elevationScatter = scans[scanIndex].getElevationScatter();
	std::vector<int> intFiller;
	std::vector<double> dubFiller, data, weights, bg, elevation, LSSData, dumps;
	std::vector<double> ra, dec;
	std::vector<std::vector<double> > bgData, bgWeights;

	elevation = scans[scanIndex].getElevation();
	LSSData = scans[scanIndex].getLSSData();
	dumps = scans[scanIndex].getDataDumps();
	ra = scans[scanIndex].getRa();
	dec = scans[scanIndex].getDec();

	bg.resize(size, 999999);
	bgData.resize(size, dubFiller);
	bgWeights.resize(size, dubFiller);

	int end = 0, start = 0, counter = 0, low, high, trueCount;
	double mean, maxDev = 0;

	RCR rcr = RCR(LS_MODE_DL);

	//FIND THE MINIMUM AND MAXIMUM INDEX VALUES FOR A BASELINE LENGTH STARTING AT INDEX start
	for (; start < size - 2; start++)
	{
		while (end < size - 1 && elevation[end] - elevation[start] < baseline)
		{
			end++;
		}
		end--;
		if (end == size - 2 && elevation[end + 1] - elevation[start] < baseline)
		{
			end++;
		}
		//CALCULATE THE LOCAL MODEL VALUES
		if (end != start)
		{
			//calculateBaselines(true, start, end, size, elevationScatter, elevation, LSSData, dumps, bgData, bgWeights, ra, dec);
		}

	}
	end = size - 1;
	start = size - 1;
	for (; end > 1; end--)
	{
		while (elevation[end] - elevation[start] < baseline && start > 0)
		{
			start--;
		}
		start++;
		if (start == 1 && elevation[end] - elevation[start - 1] < baseline)
		{
			start--;
		}
		if (start != end)
		{
			//calculateBaselines(false, start, end, size, elevationScatter, elevation, LSSData, dumps, bgData, bgWeights, ra, dec);
		}
	}

	for (int i = 0; i < bgData.size(); i++)
	{
		trueCount = bgData[i].size();
		if (trueCount == 0)
		{
			mean = 999999;
		}
		else
		{
			for (int j = 0; j < bgData[i].size(); j++)
			{
				if (bgData[i][j] == NAN || bgData[i][j] != bgData[i][j] || bgWeights[i][j] == NAN || bgWeights[i][j] != bgWeights[i][j])
				{
					trueCount--;
				}
				else
				{
					data.push_back(bgData[i][j]);
					weights.push_back(bgWeights[i][j]);
				}
			}
			if (trueCount == 0)
			{
				mean = 999999;
			}
			else if (trueCount == 1)
			{
				mean = bgData[i][0];
			}
			else
			{
				rcr.performBulkRejection(weights, data);
				mean = rcr.result.mu;
			}
			data.clear();
			weights.clear();
		}
		if (trueCount > 0)
		{
			bg[i] = mean;
		}
		else
		{
			bg[i] = 999999;
		}
	}

	for (int k = 0; k < size; k++)
	{
		if (bg[k] == 999999)
		{
			low = Tools::max(0, k - 1);
			high = Tools::min(k + 1, size - 1);
			while (bg[low] == 999999 && low > 1)
			{
				low--;
			}
			while (bg[high] == 999999 && high < size - 1)
			{
				high++;
			}
			if (bg[low] != 999999 && bg[high] != 999999)
			{
				bg[k] = bg[low] + (bg[high] - bg[low]) / (high - low)*(k - low);
			}
			else if (bg[low] == 999999)
			{
				bg[low] = bg[high];
			}
			else if (bg[high] == 999999)
			{
				bg[high] = bg[low];
			}
		}
	}
	scans[scanIndex].removeElevation(bg);
}

//lss rfi
bool Processor::verifyOutlierState(SurfaceType modelType, int i_0, int j_0, double thetaFWHM, std::vector<double> &coef, std::vector<double> &inRange)
{

	double thetaPerp;
	std::vector<double> weights, quadCoefs, ra, dec;
	std::vector<bool> checks;
	checks.resize(scans[i_0].getSize(), 1);
	weights.resize(scans[i_0].getSize(), 1);
	ra = scans[i_0].getRa();
	dec = scans[i_0].getDec();

	double secondDerivative = calculateSecondDerivitive(modelType, i_0, j_0, thetaFWHM, coef);

	if (scansInRa)
	{
		quadCoefs = Tools::quadraticRegressionPivot(checks, weights, ra, dec);
	}
	else
	{
		quadCoefs = Tools::quadraticRegressionPivot(checks, weights, dec, ra);
	}

	double b = quadCoefs[1];
	double c = quadCoefs[0];

	if (scansInRa)
	{
		thetaPerp = 1 / Tools::signedArcTan(b + 2 * c * scans[i_0].getRa(j_0), 1);
	}
	else
	{
		thetaPerp = 1 / Tools::signedArcTan(b + 2 * c * scans[i_0].getDec(j_0), 1);
	}

	if (thetaPerp * secondDerivative > 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}
bool Processor::planeApplication(int minScanCount, int minPointsInScan, int minPointsInAp, std::vector<std::vector<int>> scanNumbers, int counter, MapTypes mapType)
{
	bool minApCheck = false, minScanCheck = false, minWScalePointsCheck = false;
	double wPointsInScan;
	int pointCount;
	int possibles;
	int scanIndex, dataIndex;
	std::vector<int> possiblePoints;

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
void Processor::rejectModelPoints(int i_0, int j_0, double thetaMinPrime, double thetaFWHM, double avgScatter, std::vector<double> &inRange, std::vector<bool> &checks)
{
	double largestGap, thetaWPrime, alpha;
	double secondDeriv, thetaPerp, signedDelta;
	double sigma = 999999, delta, max = -999999, xPrime, yPrime, zPrime, wPrime, thetaPrime;
	double w = 0, ww = 0, wzz = 0, wDelta;
	double wVarience;
	int count, iTemp, jTemp, maxI;
	std::vector<double> coef, deltas;
	std::vector<int> deltasIndex;
	std::vector<std::vector<int> > scanNumbers;
	std::ofstream SurPoints;
	SurfaceType modelType;

	deltas.reserve(inRange.size() / 6);
	deltasIndex.reserve(inRange.size() / 6);

	largestGap = scans[i_0].getLSSThetaGap(j_0); //circleThetaGap(i_0, j_0, true);

	while (sigma > avgScatter)
	{
		//largestGap = lssCircleThetaGap(i_0, j_0, inRange, checks);// in radians
		thetaWPrime = Tools::max(thetaMinPrime, Tools::min(((4.0 / 3.0)*largestGap*toDeg / psfFWHM), thetaFWHM));
		alpha = -1.0*std::log(2) / (std::log(std::cos(M_PI*thetaWPrime / (4.0*thetaFWHM))));

		count = 0;
		sigma = 0;
		maxI = -1;

		scanNumbers = calculateScanNumbers(count, inRange, checks);
		modelType = determineSurfaceType(scanNumbers, mapType);
		coef = determineModelCoef(modelType, count, inRange, checks, alpha, scans[i_0].getDec(j_0), scans[i_0].getRa(j_0));

		//CALCULATE DELTAS OFF OF MODEL
		for (int i = 0; i < inRange.size() / 6; i++)
		{
			if (checks[i] == false)
			{
				continue;
			}

			iTemp = inRange[4 + i * 6];
			jTemp = inRange[5 + i * 6];

			xPrime = (scans[iTemp].getDec(jTemp) - scans[i_0].getDec(j_0));
			yPrime = (scans[iTemp].getRa(jTemp) - scans[i_0].getRa(j_0));
			zPrime = applySurfaceModel(xPrime, yPrime, inRange[i * 6], coef); //	WE NEED TO CHANGE MODELS IF IT GIVES INFINITY
			wPrime = std::pow(inRange[1 + i * 6], alpha);

			if (zPrime != zPrime)
			{
				Debugger::print("Info", "NAN");
			}

			signedDelta = inRange[6 * i] - zPrime;
			delta = std::abs(inRange[6 * i] - zPrime);

			thetaPerp = calculateThetaPerp(iTemp, jTemp);
			secondDeriv = calculateSecondDerivitive(modelType, iTemp, jTemp, thetaPerp, coef); //possibly not being calculated for m6

			w += wPrime;
			ww += wPrime * wPrime;
			wzz += wPrime * delta*delta;
			wVarience = wPrime * delta * delta;


			if (wVarience > max)
			{
				max = wVarience;
				maxI = i;
			}

			//if ((signedDelta > 0 && secondDeriv < 0) || (signedDelta < 0 && secondDeriv > 0))
			//{
			//	if (wVarience > max)
			//	{
			//		max = wVarience;
			//		maxI = i;
			//	}
			//}
		}

		double deltaTemp;
		switch (modelType) {
		case(M10):
			deltaTemp = 10;
			break;
		case(M6):
			deltaTemp = 6;
			break;
		case(M3):
			deltaTemp = 3;
			break;
		case(M0):
			deltaTemp = 0;
			break;
		}

		sigma = std::sqrt(wzz / (w - (deltaTemp)*(ww / w)));

		//std::cout << sigma << "\t" << sigma/count << "\t" << avgScatter << "\n";

		if ((sigma) < avgScatter)
		{
			break;
		}

		if (count > 1 && maxI != -1)
		{
			if (i_0 == 59 && j_0 == 10)
			{
				SurPoints.open("SurPoints.txt");
				SurPoints << sigma << "\t" << avgScatter << "\t" << inRange[6 * maxI + 4] << "\t" << inRange[6 * maxI + 5] << "\t" << "rejected\n";
				SurPoints.close();
			}

			checks[maxI] = false;
			max = 0;
		}
		else
		{
			sigma = -999999;
		}

	}
}
double Processor::lssRecursivePolyFit(Survey &survey, int i_0, int j_0, double bgScale, double wScale)
{
	int pointCandidates, count = 0, holdI, holdJ;
	double distance;
	double avgScatter = 0;
	double xPrime, yPrime, zPrime;
	SurfaceType modelType;
	std::vector<int> possibles;
	std::vector<bool> checks;
	std::vector<double> inRange, coef;
	std::vector<std::vector<int> > scanNumbers;

	double preScanDist = Tools::getGCDistance(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), scans[i_0].getDec(0), scans[i_0].getRa(0), partSetProcSSS.centerDecDeg)*toDeg;
	double proScanDist = Tools::getGCDistance(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), scans[i_0].getDec(scans[i_0].getSize() - 1), scans[i_0].getRa(scans[i_0].getSize() - 1), partSetProcSSS.centerDecDeg)*toDeg;
	double thetaBG = bgScale;// UNITS OF BEAMS
	double thetaEnd = std::max(proScanDist, preScanDist) / psfFWHM;// UNITS OF BEAMS
	double thetaFWHM = std::max(std::min(thetaBG / 3.0, 2.0*thetaEnd / 3.0), 1.0);// UNITS OF BEAMS

	double largestGap = scans[i_0].getLSSThetaGap(j_0);
	double thetaMinPrime = wScale * thetaFWHM;// 2 * thetaFWHM / 3.0;
	double thetaWPrime = Tools::max(thetaMinPrime, Tools::min(((4.0 / 3.0)*largestGap*toDeg / psfFWHM), thetaFWHM));
	double alpha = -1.0*std::log(2) / (std::log(std::cos(M_PI*thetaWPrime / (4.0*thetaFWHM))));

	std::vector<double> LSS2DScatter = survey.getLSS2DScatter();
	
	pointCandidates = findPossLSS(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), thetaFWHM*psfFWHM, possibles);
	inRange.reserve(pointCandidates * 6); //THIS MAY BE CHANGED TO RESERVE INSTEAD OF RESIZE
	checks.resize(pointCandidates, true);

	//CALCULATE EXPECTED SCATTER OVER GIVEN RADIUS
	for (int k = 0; k < pointCandidates; k++)
	{
		holdI = possibles[k * 2];
		holdJ = possibles[k * 2 + 1];

		distance = Tools::getGCDistance(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg)*toDeg;

		if (distance < thetaFWHM*psfFWHM)
		{
			//STORES ALL THE INFORMATION OF THE ACCEPTED DATA POINT

			inRange.push_back(scans[holdI].getLSSData(holdJ));
			inRange.push_back(cos((M_PI*distance) / (thetaFWHM*psfFWHM*2.0)));
			inRange.push_back(scans[holdI].getDataDumps(holdJ));
			inRange.push_back(LSS2DScatter[scans[holdI].getScanNumberInSurvey()]);
			inRange.push_back(holdI);
			inRange.push_back(holdJ);

			count++;
			avgScatter += LSS2DScatter[scans[holdI].getScanNumberInSurvey()];
		}
	}
	avgScatter = avgScatter / count;//; // (double)count;
	inRange.erase(inRange.begin() + 6 * count, inRange.begin() + inRange.size());
	checks.erase(checks.begin() + 1 * count, checks.begin() + checks.size());

	//BEGIN REJECTING POINTS UNTIL NOISE WITHIN SCATTER
	rejectModelPoints(i_0, j_0, thetaMinPrime, thetaFWHM, avgScatter, inRange, checks);

	//FINAL FIT
	count = 0;
	for (int i = 0; i < scanNumbers.size(); i++)
	{
		scanNumbers[i].clear();
	}
	scanNumbers.clear();
	scanNumbers = calculateScanNumbers(count, inRange, checks);
	modelType = determineSurfaceType(scanNumbers, mapType);

	//largestGap = lssCircleThetaGap(i_0, j_0, inRange, checks);
	thetaWPrime = Tools::max(thetaMinPrime, Tools::min(((4.0 / 3.0)*largestGap*toDeg / psfFWHM), thetaFWHM));
	alpha = -1.0*std::log(2) / (std::log(std::cos(M_PI*thetaWPrime / (4.0*thetaFWHM))));

	coef = determineModelCoef(modelType, count, inRange, checks, alpha, scans[i_0].getDec(j_0), scans[i_0].getRa(j_0));

	xPrime = 0.0;
	yPrime = 0.0;
	zPrime = applySurfaceModel(xPrime, yPrime, scans[i_0].getLSSData(j_0), coef);

	std::vector<double> results = LSSCollectResults(i_0, j_0, thetaFWHM, LSS2DScatter, coef, checks, inRange);
	//LSSRFICollectResults
	//scans[i_0].setLSSThetaWPrime(j_0, thetaWPrime);
	//scans[i_0].setLSSCorrelation(thetaWPrime); //in beams
	//scans[i_0].setLSSSummedWeights(wDumpSum/count);
	return zPrime;
}
double Processor::applySurfaceModel(double xPrime, double yPrime, double zOrig, std::vector<double> &coef)
{
	double zPrime;
	switch (coef.size())
	{
	case 10:
		zPrime = coef[0] + coef[1] * xPrime + coef[2] * yPrime + coef[3] * xPrime*xPrime + coef[4] * yPrime*yPrime + coef[5] * xPrime*yPrime + coef[6] * xPrime*xPrime*xPrime + coef[7] * yPrime*yPrime*yPrime + coef[8] * xPrime*xPrime*yPrime + coef[9] * xPrime*yPrime*yPrime;
		break;
	case 6:
		zPrime = coef[0] + coef[1] * xPrime + coef[2] * yPrime + coef[3] * xPrime*xPrime + coef[4] * yPrime*yPrime + coef[5] * xPrime*yPrime;
		break;
	case 3:
		zPrime = coef[0] + coef[1] * xPrime + coef[2] * yPrime;
		break;
	case 0:
		zPrime = zOrig;
		break;
	}
	return zPrime;
}
double Processor::calculateSecondDerivitive(SurfaceType modelType, int i_0, int j_0, double theta, std::vector<double> &coef)
{
	int Ni, Nj;
	int iTemp, jTemp;
	int count = 0;

	double x = scans[i_0].getDec(j_0);
	double y = scans[i_0].getRa(j_0);

	double secondDerivative = 0;

	switch (modelType)
	{
	case M10:
		Ni = 3;
		break;
	case M6:
		Ni = 2;
		break;
	case M3:
		Ni = 1;
		break;
	case M0:
		Ni = 0;
		break;
	}

	Nj = Ni;

	for (int i = 0; i <= Ni; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			secondDerivative += coef[count] * ((i*(i - 1))*pow(x, (i - 2))*pow(y, j)*(1 / pow(std::cos(theta), 2)) + 2 * i*j*pow(x, i - 1)*pow(y, j - 1)*(1 / std::cos(theta))*(1 / std::sin(theta)) + j * (j - 1)*pow(x, i)*pow(y, j - 2)*(1 / pow(std::sin(theta), 2)))*(pow(std::cos(theta), 2)*pow(std::sin(theta), 2));
			count++;
		}
	}
	return secondDerivative;
}
double Processor::calculateThetaPerp(int iTemp, int jTemp)
{
	double thetaPerp, x;
	std::vector<double> xVec, yVec, wVec, coef;
	std::vector<bool> checks;

	if (scansInRa)
	{
		xVec = scans[iTemp].getRa();
		yVec = scans[iTemp].getDec();
		x = scans[iTemp].getRa(jTemp);
	}
	else
	{
		xVec = scans[iTemp].getDec();
		yVec = scans[iTemp].getRa();
		x = scans[iTemp].getDec(jTemp);
	}

	checks.resize(xVec.size(), true);
	wVec.resize(xVec.size(), 1.0);

	coef = Tools::quadraticRegressionPivot(checks, wVec, xVec, yVec);

	thetaPerp = -1 / Tools::signedArcTan(coef[1] + 2 * coef[0] * x, 1);
	return thetaPerp;
}
std::vector<std::vector<int>> Processor::calculateScanNumbers(int & count, std::vector<double> & inRange, std::vector<bool> &checks)
{
	int holdI, holdJ;
	bool contained;
	std::vector<int> intFiller;
	std::vector<std::vector<int> > scanNumbers;
	//DETERMINE WHICH MODELING SCHEME CAN BE USED
	for (int i = 0; i < scanNumbers.size(); i++)
	{
		scanNumbers[i].clear();
	}
	scanNumbers.clear();
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		if (checks[i] == false)
		{
			continue;
		}
		count++;
		contained = false;
		holdI = inRange[6 * i + 4];
		holdJ = inRange[6 * i + 5];

		for (int j = 0; j < scanNumbers.size(); j++)
		{
			if (scanNumbers[j][0] == holdI)
			{
				scanNumbers[j].push_back(holdJ);
				contained = true;
			}
		}
		if (!contained)
		{
			scanNumbers.push_back(intFiller);
			scanNumbers[scanNumbers.size() - 1].push_back(holdI);
			scanNumbers[scanNumbers.size() - 1].push_back(holdJ);
		}
	}
	return scanNumbers;
}
std::vector<double> Processor::LSSCollectResults(int i_0, int j_0, double thetaFWHM, std::vector<double> & LSS2DScatter, std::vector<double> &coef, std::vector<bool> & checks, std::vector<double> &inRange)
{
	double x0 = scans[i_0].getDec(j_0);
	double y0 = scans[i_0].getRa(j_0);
	double centerLSSFlux = applySurfaceModel(0, 0, scans[i_0].getLSSData(j_0), coef);
	double correlatedScatter = LSS2DScatter[scans[i_0].getScanNumberInSurvey()];

	bool countAbove, countBelow, rejectPoints, fluxAbove;
	double dd;
	double distSq, wDistSq, pointCount, wPointCount;
	double LMWeight, LMValue;//local model value
	//double centerLSSFlux = applySurfaceModel(0, 0, z0, coef);
	double jLSSFlux;

	double pw; //proximity weight -- cos(M_PI*dist/2*rfiScale)^2
	double dw; //dump weight -- dataDump
	double pwSum = 0, dwSum = 0;
	double pwAvg = 0, dwAvg = 0;
	double pwTemp, dwTemp;

	double pwDelta;
	double pwDeltaCounter;

	double rfiCountHold = 0, pointDistance;
	double iDec, iRa, iFlux, jDec, jRa, jFlux;
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

		iDec = scans[iHoldI].getDec(iHoldJ) - x0;
		iRa = scans[iHoldI].getRa(iHoldJ) - y0;
		iFlux = scans[iHoldI].getLSSData(iHoldJ);

		if (checks[i])
		{
			distSq = 0.0;
			wDistSq = 0.0;
			pointCount = 0.0;
			wPointCount = 0.0;

			pw = inRange[6 * i + 1];

			LMValue = applySurfaceModel(iDec, iRa, iFlux, coef);

			for (int j = 0; j < inRange.size() / 6; j++)
			{
				jHoldI = inRange[6 * j + 4];
				jHoldJ = inRange[6 * j + 5];

				jDec = scans[jHoldI].getDec(jHoldJ) - x0;
				jRa = scans[jHoldI].getRa(jHoldJ) - y0;
				jFlux = scans[jHoldI].getLSSData(jHoldJ);

				pwTemp = inRange[6 * j + 1];
				pwDelta += pow(pwTemp - pwAvg, 2);
				pwDeltaCounter++;

				jLSSFlux = applySurfaceModel(jDec, jRa, jFlux, coef);

				if (wCorrMap)//We determine the first correction factor and the correlation scale for photometry here.
				{
					rfiCountHold = 0.0;
					countAbove = (centerLSSFlux > sqrt(2.0)*correlatedScatter) ? true : false;
					fluxAbove = (jLSSFlux > sqrt(2.0)*correlatedScatter) ? true : false;

					dd = scans[jHoldI].getDataDumps(jHoldJ);
					pointDistance = Tools::getGCDistance(iDec, iRa, jDec, jRa, partSetProcSSS.centerDecDeg)*toDeg;//CHECK THAT fwhmRfi*toDeg is right units

					if (checks[j] && (pointDistance < thetaFWHM*psfFWHM) && ((countAbove && fluxAbove) || (!countAbove && !fluxAbove))) {//Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are (Case 1: above OR Case 2: below) the 2d-scatter. 
						rfiCountHold += 1;
						if ((pointDistance != 0))
						{
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

SurfaceType Processor::determineSurfaceType(std::vector<std::vector<int>> & scanNumbers, MapTypes mapType)
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
std::vector<double> Processor::determineModelCoef(SurfaceType modelType, int pointCount, std::vector<double>& inRange, std::vector<bool> &checks, double alpha, double dec, double ra)
{
	std::vector<double> coef;
	switch (modelType)
	{
	case M10:
		coef = surfaceFit10(inRange, checks, alpha, dec, ra);
		break;
	case M6:
		coef = surfaceFit6(inRange, checks, alpha, dec, ra);
		break;
	case M3:
		coef = surfaceFit3(inRange, checks, alpha, dec, ra);
		break;
	case M0:
		coef.push_back(1);
		break;
	}
	return coef;
}
std::vector<double> Processor::surfaceFit3(std::vector<double>& inRange, std::vector<bool> &checks, double alpha, double i_0, double j_0)
{
	int columnCount = 3, pivotIndex, intersect;
	int scanIndex, dataIndex, iPointCandidates;

	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double w = 0,
		xw = 0, yw = 0,
		xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, wLocalW = 0, wLocalW1 = 0, wLocalWInv = 0, finalAnswer;
	double ww = 0, sumWSquared = 0;
	double weightingFunction;

	
	double withinRfiDec, withinBeamDec, withinRfiRa, withinBeamRa;
	double rfiDistance;
	int iTemp, jTemp;

	std::vector<double> A;
	std::vector<double> b;
	std::vector<double> coef;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);

	for (int i = 0; i < inRange.size() / 6; i++)
	{

		if (checks[i] == false)
		{
			continue;
		}

		intersect = 0;
		iTemp = inRange[4 + i * 6];
		jTemp = inRange[5 + i * 6];

		xHold = (scans[iTemp].getDec(jTemp) - i_0);
		yHold = (scans[iTemp].getRa(jTemp) - j_0);
		zHold = inRange[i * 6];

		wHold = pow(inRange[1 + i * 6], alpha);
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

	coef = Tools::matrixSolver(columnCount, A, b);

	return coef;
}
std::vector<double> Processor::surfaceFit6(std::vector<double>& inRange, std::vector<bool> &checks, double alpha, double i_0, double j_0)
{
	int columnCount = 6, pivotIndex, intersect;
	int scanIndex, dataIndex, iPointCandidates;
	double rfiDistance, distanceToCenter;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;

	double w = 0, xw = 0, yw = 0, xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, finalAnswer;
	double ww = 0, sumWSquared = 0, wLocalW = 0, wLocalW1 = 0, wLocalWInv = 0;
	double weightingFunction;

	//weightingFunction = processingWeightingFunction;
	//weightingFunction = 0.33;

	//double alpha = (log(0.5)) / (log(cos(M_PI*thetaWPrime / (4.0 *thetaFWHM))));
	double withinRfiDec, withinBeamDec, withinRfiRa, withinBeamRa;
	int iTemp, jTemp;

	std::vector<double> A;
	std::vector<double> b;
	std::vector<double> coef;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	for (int i = 0; i < inRange.size() / 6; i++)
	{
		if (checks[i] == false)
		{
			continue;
		}

		intersect = 0;
		iTemp = inRange[4 + i * 6];
		jTemp = inRange[5 + i * 6];

		xHold = (scans[iTemp].getDec(jTemp) - i_0);
		yHold = (scans[iTemp].getRa(jTemp) - j_0);
		zHold = inRange[i * 6];

		wHold = pow(inRange[1 + i * 6], alpha);
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

	coef = Tools::matrixSolver(columnCount, A, b);

	return coef;

}
std::vector<double> Processor::surfaceFit10(std::vector<double>& inRange, std::vector<bool> &checks, double alpha, double i_0, double j_0)
{
	int columnCount = 10, pivotIndex, intersect;
	int scanIndex, dataIndex, iPointCandidates;
	double distanceToCenter, rfiDistance;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower;
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
	std::vector<double> coef;
	std::vector<int> possibles;
	std::vector<std::vector<double> > finalHold;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	double weightingFunction;

	double withinRfiDec, withinBeamDec, withinRfiRa, withinBeamRa;
	int iTemp, jTemp;

	for (int i = 0; i < inRange.size() / 6; i++) //Was i < pointCount
	{
		if (checks[i] == false)
		{
			continue;
		}

		intersect = 0;
		intersect = 0;
		iTemp = inRange[4 + i * 6];
		jTemp = inRange[5 + i * 6];

		xHold = (scans[iTemp].getDec(jTemp) - i_0);
		yHold = (scans[iTemp].getRa(jTemp) - j_0);
		zHold = inRange[i * 6];

		wHold = pow(inRange[1 + i * 6], alpha);
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

	coef = Tools::matrixSolver(columnCount, A, b);

	if (coef[0] != coef[0] || coef[4] != coef[4])
	{
		Debugger::print("Info", "Stop");
	}

	return coef;
}

//lss leveling
void Processor::levelLSSData(std::vector<Survey> & surveys)
{
	PartitionSet partSetI, partSetJ;
	std::vector<Scan> scansI, scansJ;
	int pointCandidates;
	bool inBounds = false;
	double ra, dec, largestGap;
	double zPrime, scatterLSS;
	double mu, sigma, delta;
	std::vector<int> possibles;
	std::vector<double> LSS2DScatter;
	std::vector<double> diff, weights;
	std::vector<double> dubFiller;
	std::vector<std::vector<double>> mu_ij, sigma_ij;

	dubFiller.resize(surveys.size());
	mu_ij.resize(surveys.size(), dubFiller);
	sigma_ij.resize(surveys.size(), dubFiller);

	RCR rcr = RCR(LS_MODE_DL);

	for (int i = 0; i < surveys.size(); i++)
	{
		scansI = surveys[i].getScans();
		LSS2DScatter = surveys[i].getLSS2DScatter();
		for (int j = 0; j < surveys.size(); j++)
		{
			if (i == j)
			{
				continue;
			}
			scansJ = surveys[j].getScans();
			partSetJ = surveys[j].getPartSetProcSSS();


			for (int k = 0; k < scansI.size(); k++)
			{
				for (int l = 0; l < scansI[k].getSize(); l++)
				{
					ra = scansI[k].getRa(l);
					dec = scansI[k].getDec(l);
					largestGap = scansI[k].getLSSThetaGap(l);
					inBounds = checkEdgeCriteria(surveys[i].getProcessingCoordinate(), partSetJ, ra, dec);
					if (inBounds)
					{
						zPrime = surfaceModel(scansJ, partSetJ, ra, dec, largestGap, wScaleBW) - scansI[k].getLSSData(l);
						scatterLSS = LSS2DScatter[k];

						if (zPrime == zPrime)
						{
							diff.push_back(zPrime);
							weights.push_back(1 / std::pow(scatterLSS, 2));
						}
					}
				}
			}

			rcr.performBulkRejection(weights, diff);
			mu_ij[i].push_back(rcr.result.mu);
			sigma_ij[i].push_back(rcr.result.sigma);
		}
	}

	double num = 0, denom = 0;
	double maxDiff = -999999;
	double devNew, devOld, devDiff;
	std::vector<double> deviations;
	deviations.resize(surveys.size(), 0.0);
	//iterate through deviations until the most a single one changes is 0.01
	while (maxDiff > 0.01)
	{
		devDiff = 0;
		for (int i = 1; i < surveys.size(); i++)//always let one of the deviations be zero
		{
			devOld = deviations[i];
			num = 0;
			denom = 0;
			for (int j = i + 1; j < surveys.size(); j++)
			{
				num += (deviations[j] + mu_ij[i][j]) / std::pow(sigma_ij[i][j], 2);
				denom += 1 / std::pow(sigma_ij[i][j], 2);
			}
			deviations[i] = num / denom;
			devNew = deviations[i];
			devDiff = std::abs(devNew - devOld);
			if (devDiff > maxDiff)
			{
				maxDiff = devDiff;
				Debugger::print("Info", maxDiff);
			}
		}
	}
	

	//WE THEN NEED TO SUBTRACT THESE VALUES OFF OF THE GRID AND EXCISE BAD POINTS
	subtractDeviations(surveys, deviations);


	//scansI.removeLSSPoints(checks);
}
void Processor::subtractDeviations(std::vector<Survey> & surveys, std::vector<double> & deviations)
{
	std::vector<Scan> scansI;
	std::vector<double> fluxVec;
	double devI, devMin, min = 999999;
	for (int i = 1; i < deviations.size(); i++)
	{
		if (deviations[i] < min)
		{
			min = deviations[i];
		}
	}
	for (int i = 1; i < surveys.size(); i++)
	{
		scansI = surveys[i].getScans();
		devI = deviations[i];
		for (int j = 0; j < scansI.size(); j++)
		{
			fluxVec = scansI[j].getLSSData();
			for (int k = 0; k < fluxVec.size(); k++)
			{
				fluxVec[k] = fluxVec[k] - (devI - min);
			}
		}
		surveys[i].setScans(scansI);
	}
}
double Processor::surfaceModel(std::vector<Scan> scansJ, PartitionSet partSetJ, double ra, double dec, double largestGap, double wScale)
{
	int pointCandidates, holdI, holdJ, count = 0;
	double distance, scale, alpha;
	double xPrime, yPrime, zPrime;

	SurfaceType modelType;
	std::vector<int> possibles;
	std::vector<double> inRange, coef;
	std::vector<bool> checks;

	std::vector<std::vector<int> > scanNumbers;
	pointCandidates = findPossLSS(dec, ra, psfFWHM, possibles);
	inRange.resize(pointCandidates * 6); //THIS MAY BE CHANGED TO RESERVE INSTEAD OF RESIZE
	checks.resize(pointCandidates, true);

	//CALCULATE EXPECTED SCATTER OVER GIVEN RADIUS
	for (int k = 0; k < pointCandidates; k++)
	{
		holdI = possibles[k * 2];
		holdJ = possibles[k * 2 + 1];

		distance = Tools::getGCDistance(dec, ra, scansJ[holdI].getDec(holdJ), scansJ[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg)*toDeg;

		if (distance < psfFWHM)
		{
			//STORES ALL THE INFORMATION OF THE ACCEPTED DATA POINT
			inRange[count * 6] = scans[holdI].getLSSData(holdJ);
			inRange[count * 6 + 1] = cos((M_PI*distance) / (psfFWHM*2.0));
			inRange[count * 6 + 2] = 999999;// THESE ARE SET SO THAT THE INDEXING OF THE SURFACE FIT FUNCTION STILL WORKS;
			inRange[count * 6 + 3] = 999999;
			inRange[count * 6 + 4] = holdI;
			inRange[count * 6 + 5] = holdJ;
			count++;
		}
	}

	for (int i = 0; i < scanNumbers.size(); i++)
	{
		scanNumbers[i].clear();
	}
	scanNumbers.clear();
	scanNumbers = calculateScanNumbers(count, inRange, checks);
	modelType = determineSurfaceType(scanNumbers, mapType);

	scale = Tools::max(wScale, Tools::min((4.0 / 3.0)*largestGap / (psfFWHM*toRad), 1.0));

	alpha = -1.0*std::log(2) / (std::log(std::cos(M_PI*scale / (4.0*psfFWHM))));

	coef = determineModelCoef(modelType, count, inRange, checks, alpha, dec, ra);

	xPrime = 0.0;// (scans[iTemp].getDec(jTemp) - scans[i_0].getDec(j_0));
	yPrime = 0.0;// (scans[iTemp].getRa(jTemp) - scans[i_0].getRa(j_0));
	zPrime = applySurfaceModel(xPrime, yPrime, NAN, coef);

	return zPrime;
}
bool Processor::checkEdgeCriteria(Coordinates coordinate, PartitionSet partSetProcSSS, double ra, double dec)
{
	bool withinBounds = false;
	MapTypes mapType;

	std::vector<double> coordinateCheck, edgeOneParameters, edgeTwoParameters, edgeThreeParameters, edgeFourParameters;
	mapType = partSetProcSSS.mapType;
	bool tracking = partSetProcSSS.tracking;

	coordinateCheck = determineMappedEdges(ra, dec, tracking, partSetProcSSS, coordinate, mapType);

	ra = coordinateCheck[0];
	dec = coordinateCheck[1];

	mapType = partSetProcSSS.mapType;

	edgeOneParameters = partSetProcSSS.edgeOneParameters;
	edgeTwoParameters = partSetProcSSS.edgeTwoParameters;
	edgeThreeParameters = partSetProcSSS.edgeThreeParameters;
	edgeFourParameters = partSetProcSSS.edgeFourParameters;

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
		edgeRadius = partSetProcSSS.edgeRadius;
		medianDec = partSetProcSSS.medianDec;
		medianRa = partSetProcSSS.medianRa;
		trimSize = partSetProcSSS.trimSize;

		distance = Tools::getModGCDistance(dec, ra, medianDec, medianRa)*toDeg;

		if ((distance / psfFWHM) < ((edgeRadius - trimSize * psfFWHM) / psfFWHM))
		{
			withinBounds = true;
		}
	}
	return withinBounds;
}
std::vector<double> Processor::determineMappedEdges(double LONG, double LATI, bool tracking, PartitionSet compPartSetProc, Coordinates coordinate, MapTypes mapType)
{
	double centerProcLati, centerProcLong, centerMapLati, centerMapLong;
	double unTransLati, unTransLong, unTransLatiVert, unTransLongVert, doAngDeg, undoAngDeg;

	centerProcLati = compPartSetProc.centerLatProcDeg;
	centerProcLong = compPartSetProc.centerLongProcDeg;
	centerMapLati = compPartSetProc.centerDecDeg;
	centerMapLong = compPartSetProc.centerRaDeg;

	if (!(mCoordinate == "galactic" && pCoordinate == GALACTIC) && !(mCoordinate == "equatorial" && pCoordinate == EQUATORIAL)) // if mapping coords != processing coords.
	{
		if (pCoordinate == GALACTIC)
		{	// This is in limbo because Skynet is wrong.
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




//background
/*
bool Processor::findDuplicate(double value, std::vector<double> vectorHold)
{
bool duplicate = false;
for (int i = 0; i < vectorHold.size(); i++)
{
if (value == vectorHold[i])
{
duplicate = true;
}
}
return duplicate;
}
std::vector<double> Processor::autoFixedWRegression(int start, int end, std::vector<bool> checks, std::vector<double> weights, std::vector<double> horizAxis, std::vector<double> vertAxis, double horizAxisCoordinate, double vertAxisCoordinate)
{
//ATTEMPTS QUADRATIC REGRESSION
//IF A COEFFICENT IS UNSTABLE RESORST TO LINEAR
//IF LINEAR FALIS NO MODEL
//LAST INDEX OF toRet CONTAINS MODEL ORDER
std::vector<double> quadCoef, linCoef;
std::vector<double> toRet;
bool linearModel = false;
bool noModel = false;
quadCoef.resize(3, 0.0);
linCoef.resize(2, 0.0);

quadCoef = Tools::fixedQuadraticRegression(start, end, checks, weights, horizAxis, vertAxis, horizAxisCoordinate, vertAxisCoordinate);
linCoef = Tools::fixedRegression(start, end, checks, weights, horizAxis, vertAxis, horizAxisCoordinate, vertAxisCoordinate);

for (int i = 0; i < quadCoef.size(); i++)
{
if (quadCoef[i] != quadCoef[i])
{
linearModel = true;
}
}

for (int i = 0; i < linCoef.size(); i++)
{
if (linCoef[i] != linCoef[i])
{
noModel = true;
}
}

if (linearModel)
{
for (int i = 0; i < quadCoef.size(); i++)
{
quadCoef[i] = NAN;
}
}

if (noModel)
{
for (int i = 0; i < linCoef.size(); i++)
{
linCoef[i] = NAN;
}
}

toRet.insert(toRet.end(), quadCoef.begin(), quadCoef.end());
toRet.insert(toRet.end(), linCoef.begin(), linCoef.end());

return toRet;
}
std::vector<double> Processor::autoWRegression(int start, int end, std::vector<bool> checks, std::vector<double> weights, std::vector<double> horizAxis, std::vector<double> vertAxis)
{
//ATTEMPTS QUADRATIC REGRESSION
//IF A COEFFICENT IS UNSTABLE RESORST TO LINEAR
//IF LINEAR FALIS NO MODEL
//LAST INDEX OF toRet CONTAINS MODEL ORDER
std::vector<double> quadCoef, linCoef;
std::vector<double> toRet;
bool linearModel = false;
bool noModel = false;
quadCoef.resize(3, 0.0);
linCoef.resize(2, 0.0);

quadCoef = Tools::quadraticRegressionPivot(start, end, checks, weights, horizAxis, vertAxis);
linCoef = Tools::regressionPivot(start, end, checks, weights, horizAxis, vertAxis);

for (int i = 0; i < quadCoef.size(); i++)
{
if (quadCoef[i] != quadCoef[i])
{
linearModel = true;
}
}

for (int i = 0; i < linCoef.size(); i++)
{
if (linCoef[i] != linCoef[i])
{
noModel = true;
}
}

if (linearModel)
{
std::cout << "Quadratic Regression Unsuccessful\n";
for (int i = 0; i < quadCoef.size(); i++)
{
quadCoef[i] = NAN;
}
}

if (noModel)
{
std::cout << "Linear Regression Unsuccessful\n";
for (int i = 0; i < linCoef.size(); i++)
{
linCoef[i] = NAN;
}
}

toRet.insert(toRet.end(), quadCoef.begin(), quadCoef.end());
toRet.insert(toRet.end(), linCoef.begin(), linCoef.end());

return toRet;
}
std::vector<double> Processor::applyModel(double x, double x0, std::vector<double> coefs)
{
std::vector<double> modelVal;
modelVal.resize(2, 0.0);

//FIRST ATTEMPT QUADRATIC
modelVal[0] = coefs[0] + coefs[1] * (x - x0) + coefs[2] * pow(x - x0, 2);
//modelVal[0] = coefs[0] + coefs[1] * (x) + coefs[2] * pow(x, 2);
modelVal[1] = 2.0;

//THEN ATTEMPT LINEAR
if (modelVal != modelVal)
{
//std::cout << "Local Model: Failure to use Quadratic Coefficents\n";
modelVal[0] = coefs[3] + coefs[4] * (x - x0);
//modelVal[0] = coefs[3] + coefs[4] * (x);
modelVal[1] = 1.0;
}
else if (modelVal != modelVal)
{
//std::cout << "Local Model: Failure to use Linear Coefficents\n";
modelVal[0] = NAN;
modelVal[1] = 0.0;
}

return modelVal;
}
void Processor::calculateBaselines(bool forward, int start, int end, int size, double scatter, std::vector<double> &angDist, std::vector<double> &flux, std::vector<double> &dumps, std::vector<std::vector<double> > &bgData, std::vector<std::vector<double> > &bgWeights, std::vector<double> &ra, std::vector<double> &dec)
{
std::vector<bool> checks;
std::vector<bool> checksTemp;
std::vector<double> angDistTemp, dumpsTemp;
std::vector<double> raTemp, decTemp;
std::vector<double> zTemp;
std::vector<double> regressionLine;
std::vector<double> holderAVec, holderBVec;
std::vector<double> distinctXValues;
std::vector<double> distinctCorrespondingYValues;
std::ofstream baselineCheck;
double largest, smallest, splus = 999999, holderA, holderB, regressionHorizValue, regressionVertValue, splusTemp = 0, splusSigmaDiff = 999999;
double distinctValueHold;
double distinctYHold;
bool isLinear = false;
bool testBool;
int counter = 0, counterplus, largeIndex, low = 999999, high = -999999, smallIndex, offset, lastChecked;
int checkCount;
std::ofstream testFile;

//MAKE ALL POINTS IN THE ENTIRE SCAN NON-REJECTED FIRST
checks.resize(size, true);

//DEFINE THE ANCHOR POINTS WHEN USING FIXED REGRESSION

//THIS MEANS WE MUST RESTRICT THE QUADRATIC OR LINEAR MODEL TO FIT THROUGH THE regression VALUES
if (forward)
{
checks[start] = false;

//ANCHOR POINTS
regressionHorizValue = 0;
regressionVertValue = flux[start];
}
else
{
checks[end] = false;

//ANCHOR POINTS
regressionHorizValue = angDist[end] - angDist[start];
regressionVertValue = flux[end];
}

distinctXValues.push_back(regressionHorizValue);

angDistTemp.resize(end - start + 1);
zTemp.resize(end - start + 1);
checksTemp.resize(end - start + 1);
dumpsTemp.resize(end - start + 1);
raTemp.resize(end - start + 1);
decTemp.resize(end - start + 1);

//STORE INFORMATION FOR DATA POINTS WITHIN DESIRED BASELINE
for (int i = 0; i < angDistTemp.size(); i++)
{
angDistTemp[i] = angDist[start + i] - angDist[start];
zTemp[i] = flux[start + i];
dumpsTemp[i] = dumps[start + i];
raTemp[i] = ra[start + i];
decTemp[i] = dec[start + i];
}

//REJECT POINTS UNTIL LOCAL MODEL IS WITHIN THE NOISE
while (splus > scatter)
{
checkCount = 0;
distinctXValues.clear();
distinctXValues.push_back(regressionHorizValue);
//COLLECT ORIGINAL FLAGS FOR SUB-SECTION OF POINTS
for (int i = 0; i < checksTemp.size(); i++)
{
checksTemp[i] = checks[start + i];
distinctValueHold = angDist[start + i];
distinctYHold = flux[start + i];
testBool = findDuplicate(distinctValueHold, distinctXValues);
if ((checksTemp[i] == true) && (testBool == false) && (checkCount < 6))
{
distinctXValues.push_back(distinctValueHold);
distinctCorrespondingYValues.push_back(distinctYHold);
checkCount = distinctXValues.size();
}
}

checkCount = distinctXValues.size();
regressionLine = autoFixedWRegression(0, end - start, checksTemp, dumpsTemp, angDistTemp, zTemp, regressionHorizValue, regressionVertValue);

//BASED ON THE MODEL PARAMETERS, REMOVE THE HIGHEST DEVIATING POINT AND REPEAT
splus = 0.0;
largest = 0.0;
counter = 0;
counterplus = 0;
largeIndex = -1;

for (int k = start; k < end + 1; k++)
{
if (checks[k])
{
counter++;
holderAVec = applyModel(angDist[k], angDist[start], regressionLine);
if (flux[k] > holderAVec[0])
{
splus = splus + (flux[k] - holderAVec[0])*(flux[k] - holderAVec[0]);
counterplus++;
}
if (flux[k] - holderAVec[0] > largest)
{
largest = flux[k] - holderAVec[0];
largeIndex = k;
}
}
}

//HALT REJECTION IF THERE ARE ONLY THREE DISTINCT POINTS REMAINING
if (checkCount > (int)holderAVec[1] + 1)
{
splus = sqrt((splus / (counterplus)));
splusTemp = splus;
if ((splus > scatter) && (largeIndex>-1))
{
checks[largeIndex] = false;
}
}
else
{
splusTemp = sqrt((splus) / (counterplus));
splus = 0;
}
}

//NOW "REMOVE" ANCHOR; START ADDING POINTS INTO A MORE FLEXIBLE MODEL
splus = splusTemp;
while (splus < scatter)
{
//DETERMINE THE LOWEST AND HIGHEST NON-REJECTED POINTS
for (int k = start; k < end + 1; k++)
{
if (checks[k])
{
if (k < low)
{
low = k;
}
if (k > high)
{
high = k;
}
}
}

//RESTORE THE FLAGS TO SHOW THE REJECTED AND NON-REJECTED POINTS
checkCount = 0;
distinctXValues.clear();
//distinctXValues.push_back(regressionHorizValue);
//COLLECT ORIGINAL FLAGS FOR SUB-SECTION OF POINTS
for (int i = 0; i < checksTemp.size(); i++)
{
checksTemp[i] = checks[start + i];
distinctValueHold = angDist[start + i];
distinctYHold = flux[start + i];
testBool = findDuplicate(distinctValueHold, distinctXValues);
if ((checksTemp[i] == true) && (testBool == false) && (checkCount < 6))
{
distinctXValues.push_back(distinctValueHold);
distinctCorrespondingYValues.push_back(distinctYHold);
checkCount = distinctXValues.size();
}
}

checkCount = distinctXValues.size();

//REBUILD MODELS
regressionLine = autoWRegression(0, end - start, checksTemp, dumpsTemp, angDistTemp, zTemp);

isLinear = false;
splus = 0.0;
smallest = 999999;
counter = 0;
counterplus = 0;
smallIndex = -1;

if (forward)
{
offset = 1;
}
else
{
offset = 0;
}

for (int k = std::max(low - 1, start + offset); k < std::min(high + 2, end + offset); k++)
{
//CALCULATE VALUES FOR THE MODEL
holderAVec = applyModel(angDist[k], angDist[start], regressionLine);

//DETERMINE DEVIATION FROM MODEL AND SMALLEST *REJECTED* POINT
if (checks[k])
{
counter++;
if (flux[k] > holderAVec[0])
{
splus = splus + (flux[k] - holderAVec[0])*(flux[k] - holderAVec[0]);
counterplus++;
}
}
else
{
if (flux[k] - holderAVec[0] < smallest)
{
smallest = flux[k] - holderAVec[0];
smallIndex = k;
}
}
}

//IF MODEL HAS AT LEAST THREE POINTS, ADD VALUES BACK IN?
if (checkCount >= (int)holderAVec[1] + 1)
{
splus = sqrt((splus / (counterplus)));
if ((splus < scatter) && (smallIndex>-1))
{
checks[smallIndex] = true;
splusSigmaDiff = std::abs(scatter - splus);
lastChecked = smallIndex;
}
else
{
if (splusSigmaDiff < std::abs(scatter - splus))
{
checks[lastChecked] = false;
}
splus = 999999;
}
}
else
{
splus = 999999;
}
}

for (int i = 0; i < checksTemp.size(); i++)
{
checksTemp[i] = checks[start + i];
}

regressionLine = autoWRegression(0, end - start, checksTemp, dumpsTemp, angDistTemp, zTemp);

for (int k = start; k < end + 1; k++)
{
if (checks[k])
{
if (k < low)
{
low = k;
}
if (k > high)
{
high = k;
}
}
}

double mean = 0, stDev = 0, kurtosis = 0, dumpSum = 0;
for (int i = low; i < high + 1; i++)
{
if (checks[i])
{
mean += dumps[i] * angDist[i];
dumpSum += dumps[i];
}
}
mean = mean / dumpSum;

for (int i = low; i < high + 1; i++)
{
if (checks[i])
{
stDev += dumps[i] * pow(angDist[i] - mean, 2);
kurtosis += dumps[i] * pow(angDist[i] - mean, 4);
}
}

stDev = sqrt(stDev / dumpSum);
kurtosis = pow(kurtosis / dumpSum, .25);
bgData.reserve(high - low + 1);
bgWeights.reserve(high - low + 1);

for (int k = low; k < high + 1; k++)
{
holderBVec = applyModel(angDist[k], angDist[start], regressionLine);

if (!((holderBVec[0] > flux[k]) && (checks[k] == 0)))
{
bgData[k].push_back(holderBVec[0]);

if (holderBVec[0] == NAN || holderBVec[0] != holderBVec[0] || stDev == 0 || kurtosis == 0)
{
std::cout << "Bad Local Model\n";
bgWeights[k].push_back(NAN);
}
else if (counter == 1)
{
bgWeights[k].push_back(1.0);
}
else
{
bgWeights[k].push_back(dumpSum / (1 + pow((angDist[k] - mean) / stDev, 2) + pow((angDist[k] - mean) / kurtosis, 4)));
}
}
}
}
void Processor::calculateBG(int scanIndex, double baseline)
{
std::cout << scanIndex << "\n";
int size = scans[scanIndex].getSize();
double scatter = scans[scanIndex].getScatter();
std::vector<int> intFiller;
std::vector<double> dubFiller, data, weights, bg, angDist, flux, dumps;
std::vector<double> ra, dec;
std::vector<std::vector<double> > bgData, bgWeights;

angDist = scans[scanIndex].getAngDist();
flux = scans[scanIndex].getFlux();
dumps = scans[scanIndex].getDataDumps();
ra = scans[scanIndex].getRa();
dec = scans[scanIndex].getDec();

bg.resize(size, 999999);
bgData.resize(size, dubFiller);
bgWeights.resize(size, dubFiller);

int end = 0, start = 0, counter = 0, low, high, trueCount;
double mean, maxDev = 0;

RCR rcr = RCR(LS_MODE_DL);

//FIND THE MINIMUM AND MAXIMUM INDEX VALUES FOR A BASELINE LENGTH STARTING AT INDEX start
for (; start < size - 2; start++)
{
while (end < size - 1 && angDist[end] - angDist[start] < baseline)
{
end++;
}
end--;
if (end == size - 2 && angDist[end + 1] - angDist[start] < baseline)
{
end++;
}
//CALCULATE THE LOCAL MODEL VALUES
if (end != start)
{
calculateBaselines(true, start, end, size, scatter, angDist, flux, dumps, bgData, bgWeights, ra, dec);
}

}
end = size - 1;
start = size - 1;
for (; end > 1; end--)
{
while (angDist[end] - angDist[start] < baseline && start > 0)
{
start--;
}
start++;
if (start == 1 && angDist[end] - angDist[start - 1] < baseline)
{
start--;
}
if (start != end)
{
calculateBaselines(false, start, end, size, scatter, angDist, flux, dumps, bgData, bgWeights, ra, dec);
}
}

for (int i = 0; i < bgData.size(); i++)
{
trueCount = bgData[i].size();
if (trueCount == 0)
{
mean = 999999;
}
else
{
for (int j = 0; j < bgData[i].size(); j++)
{
if (bgData[i][j] == NAN || bgData[i][j] != bgData[i][j] || bgWeights[i][j] == NAN || bgWeights[i][j] != bgWeights[i][j])
{
trueCount--;
}
else
{
data.push_back(bgData[i][j]);
weights.push_back(bgWeights[i][j]);
}
}
if (trueCount == 0)
{
mean = 999999;
}
else if (trueCount == 1)
{
mean = bgData[i][0];
}
else
{
rcr.performBulkRejection(weights, data);
mean = rcr.result.mu;
}
data.clear();
weights.clear();
}
if (trueCount > 0)
{
bg[i] = mean;
}
else
{
bg[i] = 999999;
}
}

for (int k = 0; k < size; k++)
{
if (bg[k] == 999999)
{
low = Tools::max(0, k - 1);
high = Tools::min(k + 1, size - 1);
while (bg[low] == 999999 && low > 1)
{
low--;
}
while (bg[high] == 999999 && high < size - 1)
{
high++;
}
if (bg[low] != 999999 && bg[high] != 999999)
{
bg[k] = bg[low] + (bg[high] - bg[low]) / (high - low)*(k - low);
}
else if (bg[low] == 999999)
{
bg[low] = bg[high];
}
else if (bg[high] == 999999)
{
bg[high] = bg[low];
}
}
}
scans[scanIndex].removeBG(bg);
}
std::vector<double> Processor::calculateBGMulti(int scanIndex, double baseline)
{
std::cout << scanIndex << "\n";
int size = scans[scanIndex].getSize();
double scatter = scans[scanIndex].getScatter();
std::vector<int> intFiller;
std::vector<double> dubFiller, data, weights, bg, angDist, flux, dumps;
std::vector<double> ra, dec;
std::vector<std::vector<double> > bgData, bgWeights;

angDist = scans[scanIndex].getAngDist();
flux = scans[scanIndex].getFlux();
dumps = scans[scanIndex].getDataDumps();
ra = scans[scanIndex].getRa();
dec = scans[scanIndex].getDec();

bg.resize(size, 999999);
bgData.resize(size, dubFiller);
bgWeights.resize(size, dubFiller);

int end = 0, start = 0, counter = 0, low, high, trueCount;
double mean, maxDev = 0;

RCR rcr = RCR(LS_MODE_DL);

//FIND THE MINIMUM AND MAXIMUM INDEX VALUES FOR A BASELINE LENGTH STARTING AT INDEX start
for (; start < size - 2; start++)
{
while (end < size - 1 && angDist[end] - angDist[start] < baseline)
{
end++;
}
end--;
if (end == size - 2 && angDist[end + 1] - angDist[start] < baseline)
{
end++;
}
//CALCULATE THE LOCAL MODEL VALUES
if (end != start)
{
calculateBaselines(true, start, end, size, scatter, angDist, flux, dumps, bgData, bgWeights, ra, dec);
}

}
end = size - 1;
start = size - 1;
for (; end > 1; end--)
{
while (angDist[end] - angDist[start] < baseline && start > 0)
{
start--;
}
start++;
if (start == 1 && angDist[end] - angDist[start - 1] < baseline)
{
start--;
}
if (start != end)
{
calculateBaselines(false, start, end, size, scatter, angDist, flux, dumps, bgData, bgWeights, ra, dec);
}
}

for (int i = 0; i < bgData.size(); i++)
{
trueCount = bgData[i].size();
if (trueCount == 0)
{
mean = 999999;
}
else
{
for (int j = 0; j < bgData[i].size(); j++)
{
if (bgData[i][j] == NAN || bgData[i][j] != bgData[i][j] || bgWeights[i][j] == NAN || bgWeights[i][j] != bgWeights[i][j])
{
trueCount--;
}
else
{
data.push_back(bgData[i][j]);
weights.push_back(bgWeights[i][j]);
}
}
if (trueCount == 0)
{
mean = 999999;
}
else if (trueCount == 1)
{
mean = bgData[i][0];
}
else
{
rcr.performBulkRejection(weights, data);
mean = rcr.result.mu;
}
data.clear();
weights.clear();
}
if (trueCount > 0)
{
bg[i] = mean;
}
else
{
bg[i] = 999999;
}
}

for (int k = 0; k < size; k++)
{
if (bg[k] == 999999)
{
low = Tools::max(0, k - 1);
high = Tools::min(k + 1, size - 1);
while (bg[low] == 999999 && low > 1)
{
low--;
}
while (bg[high] == 999999 && high < size - 1)
{
high++;
}
if (bg[low] != 999999 && bg[high] != 999999)
{
bg[k] = bg[low] + (bg[high] - bg[low]) / (high - low)*(k - low);
}
else if (bg[low] == 999999)
{
bg[low] = bg[high];
}
else if (bg[high] == 999999)
{
bg[high] = bg[low];
}
}
}
return bg; //scans[scanIndex].removeBG(bg);
}
*/
/*
void Processor::performBGSubtractionMulti(Survey &survey, double bgScale)
{
	Output output;
	output.printBgSubtractionInfo(bgScale);
	characterizeSurvey(survey);

	int threadCounter = 0;
	int liveThreads = 0;
	int completedThreads = 0;
	int maxThreads = 8;

	std::vector<std::future<std::vector<double>>> futureVec;
	futureVec.resize(scans.size());
	std::vector<double> bgTemp;

	timer bgClock;
	bgClock.clockTime(true);

	if (bgScale != 0.0)
	{
		for (int i = 0; i < scans.size(); i++)
		{
			//futureVec[i] = std::async(std::launch::async, &Processor::calculateBGMulti, this, i, bgScale*psfFWHM);
		}
		for (int i = 0; i < scans.size(); i++)
		{
			//bgTemp = futureVec[i].get();
			//scans[i].removeBG(bgTemp);
		}
	}
	bgClock.clockTime(false);

	//if (bgScale != 0.0)
	//{
	//	for (int i = 0; i < scans.size(); i++)
	//	{
	//		futureVec[i] = std::async(std::launch::async, &Processor::calculateBGMulti, this, i, bgScale*psfFWHM);
	//		liveThreads++;
	//		threadCounter++;
	//		if (liveThreads >= maxThreads)
	//		{
	//			bgTemp = futureVec[i].get();
	//			scans[i].removeBG(bgTemp);
	//		}
	//		completedThreads += liveThreads;
	//		liveThreads = 0;
	//	}

	//	for (int i = 0; i < scans.size(); i++)
	//	{
	//		if (i >= completedThreads)
	//		{
	//			bgTemp = futureVec[i].get();
	//			scans[i].removeBG(bgTemp);
	//		}
	//	}
	//}

	survey.setScans(scans);
}
*/
/*
void Processor::performBGSubtraction(Survey &survey, double bgScale)
{
	Output output;
	output.printBgSubtractionInfo(bgScale);
	characterizeSurvey(survey);

	if (bgScale != 0.0)
	{
		for (int i = 0; i < scans.size(); i++)
		{
			//calculateBG(i, bgScale*psfFWHM);
		}
	}
	survey.setScans(scans);
}
*/

//OLD FUNCTIONS
//std::vector<double> rfiFit(int pointsInRange, double avgSigma, std::vector<double> &sigmaValues, std::vector<int> &checks, std::vector<double> &inRange)
//{
//	double sigmaTerm = 0, fluxShapeTerm = 0, shapeSquaredTerm = 0, shapeTerm = 0, fluxTerm = 0;
//	double f, g, sigma, sigmaSquare, flux, shape, dump;
//	std::vector<double> solution;
//	solution.resize(2);
//	for (int i = 0; i < pointsInRange; i++)
//	{
//		if (checks[i])
//		{
//			sigma = sigmaValues[i];
//			sigmaSquare = sigma * sigma;
//			flux = inRange[3 * i];
//			shape = inRange[3 * i + 1];
//			dump = inRange[3 * i + 2];
//			sigmaTerm += dump / sigmaSquare;
//			fluxShapeTerm += (dump * flux * shape) / sigmaSquare;
//			shapeSquaredTerm += (dump * shape * shape) / sigmaSquare;
//			shapeTerm += dump * shape / sigmaSquare;
//			fluxTerm += dump * flux / sigmaSquare;
//		}
//	}
//	f = (2 * sigmaTerm * fluxShapeTerm - fluxTerm * shapeTerm) / (2 * sigmaTerm * shapeSquaredTerm - shapeTerm * shapeTerm);
//	g = (fluxTerm - f * shapeTerm) / (2 * sigmaTerm);
//	solution[0] = f;
//	solution[1] = g;
//	return solution;
//}
//void Processor::rfiRemovePoints(int pointsInRange, double avgSigma, std::vector<int> & holdIVec2, std::vector<int> & holdJVec2, std::vector<double> &sigmaValues, std::vector<int> &checks, std::vector<double> &inRange)
//{
//	double largestDiff, largest, splusSigmaDiff, splus = 0;
//	double f, g;
//	double shape, weight, modelHold;
//	int counter, counterplus, largeIndexI, lastChecked;
//	std::vector<double> rfiFitSolution;
//
//
//	while (splus > avgSigma)
//	{
//		largestDiff = 0;
//		splus = 0.0;
//		largest = 0.0;
//		counter = 0;
//		counterplus = 0;
//		largeIndexI = -1;
//		rfiFitSolution = rfiFit(pointsInRange, avgSigma, sigmaValues, checks, inRange);
//		f = rfiFitSolution[0];
//		g = rfiFitSolution[1];
//		if (f > 0)
//		{
//			for (int i = 0; i < pointsInRange; i++)
//			{
//				if (checks[i])
//				{
//					shape = inRange[3 * i];
//					weight = inRange[3 * i + 1];
//					modelHold = f * weight + g;
//					counter++;
//					if (shape > modelHold)
//					{
//						splus = splus + pow(shape - modelHold, 2);
//						counterplus++;
//					}
//					if ((shape - modelHold) > largest)
//					{
//						largest = shape - modelHold;
//						largeIndexI = i;
//					}
//				}
//			}
//		}
//		else
//		{
//			for (int i = 0; i < pointsInRange; i++)
//			{
//				if (checks[i])
//				{
//					shape = inRange[3 * i];
//					weight = inRange[3 * i + 1];
//					modelHold = f * weight + g;
//					counter++;
//					splus = splus + pow(shape - modelHold, 2);
//					counterplus++;
//					if (std::abs(shape - modelHold) > largest)
//					{
//						largest = shape - modelHold;
//						largeIndexI = i;
//					}
//				}
//			}
//		}
//		if (counter > 1)
//		{
//			splus = sqrt((splus / (counterplus)));
//			if ((splus > avgSigma) && (largeIndexI > -1))
//			{
//				checks[largeIndexI] = false;
//				scans[holdIVec2[largeIndexI]].setLocalModelRejections(holdJVec2[largeIndexI], scans[holdIVec2[largeIndexI]].getLocalModelRejections(holdJVec2[largeIndexI]) + 1);
//				lastChecked = largeIndexI;
//				splusSigmaDiff = std::abs(avgSigma - splus);
//			}
//			else
//			{
//				if (splusSigmaDiff < std::abs(avgSigma - splus))
//				{
//					checks[lastChecked] = false;
//					scans[holdIVec2[largeIndexI]].setLocalModelRejections(holdJVec2[largeIndexI], scans[holdIVec2[largeIndexI]].getLocalModelRejections(holdJVec2[largeIndexI]) + 1);
//				}
//				splus = 0;
//			}
//		}
//		else
//		{
//			splus = 0;
//		}
//	}
//}
//std::vector<double> Processor::rfiRemovalMulti(int i_0, int j_0)
//{
//	//std::cout << scatter2d.size() << "\n";
//	std::vector<int> possibles;
//	bool nonRejectedPointInTestScan = false, writeModel = false;
//	int pointsInRange = 0, normalLocalWeightCounter = 0, counter, counterplus, largeIndexI, holdI, holdJ, pointCandidates = 0, pointsNotRejected = 0, pointsInLasso = 0, lastChecked, rfiCountHold = 0;
//	double avgSigma = 0, largestDiff = 0, splus = 999999, fitHold = 0, largest, f, g, shape, weight, distance, angleHold, splusSigmaDiff = 999999, toRad = M_PI / 180.0;
//	double raHold = scans[i_0].getRa(j_0), weightHold, localModelDumpSum, localModelShapeSum;
//	double shapeParam, averageInRangeDec = 0, modelHold, normalLocalWeight = 0, distanceTemp, shapeParamTemp, wResult;
//	std::vector<int> checks, holdIVec, holdJVec, holdIVec2, holdJVec2;
//	std::vector<double> dubFiller, distances, inRange, ra, dec, highs, lows, sigmaValues, rfiFitSolution, results;
//
//	//correlatedWeightMap variables
//	int aboveRfiCount, belowRfiCount;
//	bool countAbove, countBelow, rejectPoints;
//	double distSq, wDistSq, pointCount, wPointCount, stdDev, rejLimit;
//	double correlatedScatter, rfiDistance, pointDistance;
//
//	pointCandidates = findPossSSS(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), fwhmRfi*toDeg, possibles);
//	distances.resize(pointCandidates, 0.0);
//	inRange.reserve(pointCandidates * 2);
//	holdIVec2.reserve(pointCandidates);
//	holdJVec2.reserve(pointCandidates);
//
//	sigmaValues.reserve(pointCandidates);
//
//	//COLLECT POSSIBLES AND STORE INFORMATION ABOUT PROXIMITY, SCATTER, AND FLUX
//	for (int i = 0; i < pointCandidates; i++)
//	{
//		holdI = possibles[2 * i];
//		holdJ = possibles[2 * i + 1];
//
//		distance = Tools::getGCDistance(scans[i_0].getDec(j_0), scans[i_0].getRa(j_0), scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg);
//		distances[i] = distance;
//		if (distance < fwhmRfi)
//		{
//			scans[holdI].setLocalModelInstances(holdJ, scans[holdI].getLocalModelInstances(holdJ) + 1);
//			holdIVec2.push_back(holdI);
//			holdJVec2.push_back(holdJ);
//			inRange.push_back(scans[holdI].getFlux(holdJ));
//			inRange.push_back(pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2));
//			inRange.push_back(scans[holdI].getDataDumps(holdJ));
//			sigmaValues.push_back(scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()]);
//			pointsInRange++;
//			avgSigma += scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()];
//		}
//	}
//
//	checks.resize(pointsInRange, true);
//	avgSigma = (avgSigma / (double)pointsInRange);
//
//
//	//FIT A LOCAL MODEL AND REJECT MOST DEVIANT POINT UNTIL WITHIN TOLERATED SCATTER
//
//	
//	/*while (splus > avgSigma)
//	{
//		largestDiff = 0;
//		splus = 0.0;
//		largest = 0.0;
//		counter = 0;
//		counterplus = 0;
//		largeIndexI = -1;
//		rfiFitSolution = rfiFit(pointsInRange, avgSigma, sigmaValues, checks, inRange);
//		f = rfiFitSolution[0];
//		g = rfiFitSolution[1];
//		if (f > 0)
//		{
//			for (int i = 0; i < pointsInRange; i++)
//			{
//				if (checks[i])
//				{
//					shape = inRange[3 * i];
//					weight = inRange[3 * i + 1];
//					modelHold = f * weight + g;
//					counter++;
//					if (shape > modelHold)
//					{
//						splus = splus + pow(shape - modelHold, 2);
//						counterplus++;
//					}
//					if ((shape - modelHold) > largest)
//					{
//						largest = shape - modelHold;
//						largeIndexI = i;
//					}
//				}
//			}
//		}
//		else
//		{
//			for (int i = 0; i < pointsInRange; i++)
//			{
//				if (checks[i])
//				{
//					shape = inRange[3 * i];
//					weight = inRange[3 * i + 1];
//					modelHold = f * weight + g;
//					counter++;
//					splus = splus + pow(shape - modelHold, 2);
//					counterplus++;
//					if (std::abs(shape - modelHold) > largest)
//					{
//						largest = shape - modelHold;
//						largeIndexI = i;
//					}
//				}
//			}
//		}
//		if (counter > 1)
//		{
//			splus = sqrt((splus / (counterplus)));
//			if ((splus > avgSigma) && (largeIndexI > -1))
//			{
//				checks[largeIndexI] = false;
//				scans[holdIVec2[largeIndexI]].setLocalModelRejections(holdJVec2[largeIndexI], scans[holdIVec2[largeIndexI]].getLocalModelRejections(holdJVec2[largeIndexI]) + 1);
//				lastChecked = largeIndexI;
//				splusSigmaDiff = std::abs(avgSigma - splus);
//			}
//			else
//			{
//				if (splusSigmaDiff < std::abs(avgSigma - splus))
//				{
//					checks[lastChecked] = false;
//					scans[holdIVec2[largeIndexI]].setLocalModelRejections(holdJVec2[largeIndexI], scans[holdIVec2[largeIndexI]].getLocalModelRejections(holdJVec2[largeIndexI]) + 1);
//				}
//				splus = 0;
//			}
//		}
//		else
//		{
//			splus = 0;
//		}
//	}*/
//
//	rfiRemovePoints(pointsInRange, avgSigma, holdIVec2, holdJVec2, sigmaValues, checks, inRange);
//
//	//FIT FINAL LOCAL MODEL
//	rfiFitSolution = rfiFit(pointsInRange, avgSigma, sigmaValues, checks, inRange);
//	f = rfiFitSolution[0];
//	g = rfiFitSolution[1];
//
//	//I DONT GET THIS TO BE HONEST. (thanks for your honesty).
//	int checksTrue = 0;
//	localModelDumpSum = 0;
//	localModelShapeSum = 0;
//	for (int i = 0; i < checks.size(); i++)
//	{
//		if (checks[i])
//		{
//			localModelDumpSum += inRange[3 * i + 2];
//			localModelShapeSum += inRange[3 * i + 1];
//			checksTrue++;
//		}
//	}
//
//	localModelShapeSum = localModelShapeSum / (double)checksTrue;
//	//localModelShapeSum = localModelShapeSum / (double)checks.size();// True;
//
//	//ITERATE THROUGH ALL ORIGINAL COLLECTED POINTS AND STORE THEIR RA, DEC, SCANINDEX, DATA INDEX
//	int trueCount = 0;
//	holdIVec.reserve(distances.size());
//	holdJVec.reserve(distances.size());
//	for (int i = 0; i < pointCandidates; i++)
//	{
//		if (distances[i] < fwhmRfi)
//		{
//			holdI = possibles[2 * i];
//			holdJ = possibles[2 * i + 1];
//			holdIVec.push_back(holdI);
//			holdJVec.push_back(holdJ);
//			if (trueCount != i)
//			{
//				distances[trueCount] = distances[i];//ONLY STORE DISTANCES OF POINTS WITHIN RFI RADIUS, NOT OF ALL COLLECTED POSSIBLES
//			}
//			trueCount++;
//		}
//	}
//	distances.erase(distances.begin() + trueCount, distances.end());
//
//	int iHoldI, iHoldJ, jHoldI, jHoldJ;
//	double dd;
//	results.reserve(holdIVec.size() * 7);
//	for (int i = 0; i < holdIVec.size(); i++)
//	{
//		//rfiCountHold = 0;
//		normalLocalWeightCounter = 0;
//		normalLocalWeight = 0;
//		holdI = holdIVec[i];
//		holdJ = holdJVec[i];
//
//		iHoldI = holdI;
//		iHoldJ = holdJ;
//
//		if (checks[i])
//		{
//			distSq = 0.0;
//			wDistSq = 0.0;
//			pointCount = 0.0;
//			wPointCount = 0.0;
//			distance = distances[i];
//			shapeParam = pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2);
//			fitHold = f * shapeParam + g;
//
//			if (correlatedWeightMap)
//			{
//				//We determine the first correction factor and the correlation scale for photometry here.
//				rfiCountHold = 0.0;
//				countAbove = false, countBelow = false, rejectPoints = true;
//				correlatedScatter = scatter2d[scans[i_0].getSurveyNumber()][scans[i_0].getScanNumberInSurvey()];
//
//				if (scans[i_0].getFlux(j_0) > sqrt(2.0)*correlatedScatter) {
//					countAbove = true;
//				}
//				else {
//					countBelow = true;
//				}
//
//				for (int j = 0; j < holdIVec.size(); j++)
//				{
//					jHoldI = holdIVec[j];
//					jHoldJ = holdJVec[j];
//					dd = scans[jHoldI].getDataDumps(jHoldJ);
//
//					rfiDistance = Tools::getGCDistance(scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), scans[jHoldI].getDec(jHoldJ), scans[jHoldI].getRa(jHoldJ), partSetProcSSS.centerDecDeg)*toDeg;
//					if (rfiDistance < fwhmRfi*toDeg)
//					{
//						if (countAbove && scans[jHoldI].getFlux(jHoldJ) > sqrt(2.0)*correlatedScatter) {
//							rfiCountHold += 1;	//Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are above the 2d-scatter.
//						}
//						if (countBelow && scans[jHoldI].getFlux(jHoldJ) < sqrt(2.0)*correlatedScatter) {
//							rfiCountHold += 1; //Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are below the 2d-scatter. 
//						}
//					}
//					if (checks[j]) //Correlation scale here.
//					{
//						pointDistance = Tools::getGCDistance(scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), scans[jHoldI].getDec(jHoldJ), scans[jHoldI].getRa(jHoldJ), partSetProcSSS.centerDecDeg)*toDeg;
//						if ((pointDistance < fwhmRfi*toDeg) && pointDistance != 0)
//						{
//							if (countAbove && (scans[jHoldI].getFlux(jHoldJ) > sqrt(2.0)*correlatedScatter))
//							{
//								pointCount += 1.0;
//								distSq += pow(pointDistance, 2);
//
//								wPointCount += 1.0*dd;
//								wDistSq += pow(pointDistance, 2)*dd;
//							}
//							if (countBelow && (scans[jHoldI].getFlux(jHoldJ) < sqrt(2.0)*correlatedScatter))
//							{
//								pointCount += 1.0;
//								distSq += pow(pointDistance, 2);
//
//								wPointCount += 1.0*dd;
//								wDistSq += pow(pointDistance, 2)*dd;
//							}
//						}
//					}
//
//					distanceTemp = distances[j];
//					shapeParamTemp = pow(cos((M_PI*distanceTemp) / (fwhmRfi*2.0)), 2);
//					normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2); //sum of proximity weight of this point minus all other points used to make the local model squared.
//					normalLocalWeightCounter++;
//				}
//			}
//			else
//			{
//				for (int j = 0; j < holdIVec.size(); j++)
//				{
//					distanceTemp = distances[j];
//					shapeParamTemp = pow(cos((M_PI*distanceTemp) / (fwhmRfi*2.0)), 2);
//					normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2); //sum of proximity weight of this point minus all other points used to make the local model squared.
//					normalLocalWeightCounter++;
//				}
//				rfiCountHold = 1;
//			}
//
//			if (rfiCountHold != rfiCountHold || rfiCountHold == 0) {
//				rfiCountHold = 1;
//			}
//
//			if (pointsInRange == 1) {
//				wResult = inRange[2];
//			}
//			else if (pointsInRange == 2) {
//				wResult = inRange[2] + inRange[3 * 1 + 2];
//			}
//			else {
//				wResult = (localModelDumpSum / (1.0 + 2.0 * pow((shapeParam - localModelShapeSum), 2) / (normalLocalWeight / normalLocalWeightCounter)));
//				
//			}
//			
//			results.push_back(holdI);
//			results.push_back(holdJ);
//			results.push_back(wResult);
//			results.push_back(fitHold);
//			results.push_back(rfiCountHold);
//			results.push_back(wDistSq);
//			results.push_back(wPointCount);
//
//		}
//	}
//	return results;
//}
//void Processor::extraLocalModels(double decPoint, double raPoint, std::vector<std::vector<std::vector<double> > > &rfiValues, std::vector<std::vector<std::vector<double> > > &rfiWeights, std::vector<std::vector<std::vector<double> > >&rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > & rfiNs, double rfiScale, std::vector<MapTypes> mapHolder)
//{
//	int rfiCountHold = 0;
//	double comDec, comRa;
//	bool extraModels = false;
//	double standardGap = 999999;
//	std::vector<double> originalGaps;
//	originalGaps.resize(standardGaps.size());
//
//	//correlatedWeightMap variables
//	bool rejectPoints;
//	double stdDev, distSq, wDistSq, wPointCount, pointCount, rejLimit;
//	double correlatedScatter, rfiDistance, pointDistance;
//
//	for (int i = 0; i < standardGaps.size(); i++)
//	{
//		originalGaps[i] = standardGaps[i];
//		if (mapHolder[i] == DAISY)
//		{
//			standardGaps[i] = (standardGaps[i] / 2)*scans.size();
//		}
//		if (standardGaps[i] < standardGap)
//		{
//			standardGap = standardGaps[i];
//		}
//	}
//
//	standardGaps = originalGaps;
//
//
//	int extraSteps = std::floor((psfFWHM / standardGap)*pow((rfiScale), 2));
//
//	if (extraSteps % 2 != 0)
//	{
//		extraSteps += 1;
//	}
//
//	if ((rfiScale < 1.0) && (rfiScale > pow((standardGap / psfFWHM), 0.5)))
//	{
//		for (int i = -extraSteps; i <= extraSteps; i++)
//		{
//			for (int j = -extraSteps; j <= extraSteps; j++)
//			{
//				std::vector<int> possibles;
//				bool nonRejectedPointInTestScan = false, writeModel = false;
//				int pointsInRange = 0, counter, normalLocalWeightCounter = 0, counterplus, largeIndexI, holdI, holdJ, pointCandidates = 0, pointsNotRejected = 0, pointsInLasso = 0, lastChecked;
//				double avgSigma = 0, largestDiff = 0, splus = 999999, fitHold = 0, largest, f, g, shape, weight, distance, angleHold, splusSigmaDiff = 999999, toRad = M_PI / 180.0;
//				double weightHold, localModelDumpSum, localModelShapeSum;
//				double shapeParam, averageInRangeDec = 0, modelHold, normalLocalWeight = 0, distanceTemp, shapeParamTemp;
//				std::vector<int> checks, holdIVec, holdJVec;
//				std::vector<double> dubFiller, distances, inRange, ra, dec, highs, lows, sigmaValues, rfiFitSolution;
//
//				double deltaX = i * (standardGap)*(pow((1 / rfiScale), 2) - 1); //In degrees
//				double deltaY = j * (standardGap)*(pow((1 / rfiScale), 2) - 1);
//
//				double sampleWidth = psfFWHM * (1 - pow(rfiScale, 2));
//
//				distance = Tools::getGCDistance((decPoint + deltaX), (raPoint + deltaY), decPoint, raPoint, partSetProcSSS.centerDecDeg)*toDeg;
//
//				if (distance > sampleWidth * 2)
//				{
//					continue;
//				}
//
//				pointCandidates = findPossSSS((decPoint + deltaX), (raPoint + deltaY), fwhmRfi*toDeg, possibles);
//
//				distances.resize(pointCandidates, 0.0);
//				inRange.reserve(pointCandidates * 2);
//				sigmaValues.reserve(pointCandidates);
//
//				for (int i = 0; i < pointCandidates; i++)
//				{
//					holdI = possibles[2 * i];
//					holdJ = possibles[2 * i + 1];
//					distance = Tools::getGCDistance((decPoint + deltaX), (raPoint + deltaY), scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg);
//
//					distances[i] = distance;
//					if (distance < fwhmRfi)
//					{
//						inRange.push_back(scans[holdI].getFlux(holdJ));
//						inRange.push_back(pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2));
//						inRange.push_back(scans[holdI].getDataDumps(holdJ));
//						sigmaValues.push_back(scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()]);
//						pointsInRange++;
//						avgSigma += scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()];
//
//					}
//				}
//
//				checks.resize(pointsInRange, true);
//				avgSigma = (avgSigma / (double)pointsInRange);
//				while (splus > avgSigma)
//				{
//					largestDiff = 0;
//					splus = 0.0;
//					largest = 0.0;
//					counter = 0;
//					counterplus = 0;
//					largeIndexI = -1;
//					rfiFitSolution = rfiFit(pointsInRange, avgSigma, sigmaValues, checks, inRange);
//					f = rfiFitSolution[0];
//					g = rfiFitSolution[1];
//					if (f > 0)
//					{
//						for (int i = 0; i < pointsInRange; i++)
//						{
//							if (checks[i])
//							{
//								shape = inRange[3 * i];
//								weight = inRange[3 * i + 1];
//								modelHold = f * weight + g;
//								counter++;
//								if (shape > modelHold)
//								{
//									splus = splus + pow(shape - modelHold, 2);
//									counterplus++;
//								}
//								if ((shape - modelHold) > largest)
//								{
//									largest = shape - modelHold;
//									largeIndexI = i;
//								}
//							}
//						}
//					}
//					else
//					{
//						for (int i = 0; i < pointsInRange; i++)
//						{
//							if (checks[i])
//							{
//								shape = inRange[3 * i];
//								weight = inRange[3 * i + 1];
//								modelHold = f * weight + g;
//								counter++;
//
//								splus = splus + pow(shape - modelHold, 2);
//								counterplus++;
//
//								if (std::abs(shape - modelHold) > largest)
//								{
//									largest = shape - modelHold;
//									largeIndexI = i;
//								}
//							}
//						}
//					}
//					if (counter > 1)
//					{
//						splus = sqrt((splus / (counterplus)));
//						if ((splus > avgSigma) && (largeIndexI > -1))
//						{
//							checks[largeIndexI] = false;
//							lastChecked = largeIndexI;
//							splusSigmaDiff = std::abs(avgSigma - splus);
//						}
//						else
//						{
//							if (splusSigmaDiff < std::abs(avgSigma - splus))
//							{
//								checks[lastChecked] = false;
//							}
//							splus = 0;
//						}
//					}
//					else
//					{
//						splus = 0;
//					}
//				}
//				rfiFitSolution = rfiFit(pointsInRange, avgSigma, sigmaValues, checks, inRange);
//				localModelDumpSum = 0;
//				localModelShapeSum = 0;
//				int checksTrue = 0;
//				for (int i = 0; i < checks.size(); i++)
//				{
//					if (checks[i])
//					{
//						localModelDumpSum += inRange[3 * i + 2];
//						//localModelShapeSum += pow(inRange[3 * i + 1], 2);
//						localModelShapeSum += inRange[3 * i + 1];
//						checksTrue++;
//					}
//				}
//				//localModelShapeSum = localModelShapeSum / (double)checks.size();
//				localModelShapeSum = localModelShapeSum / (double)checksTrue;
//				f = rfiFitSolution[0];
//				g = rfiFitSolution[1];
//				int trueCount = 0;
//				ra.reserve(distances.size());
//				dec.reserve(distances.size());
//				holdIVec.reserve(distances.size());
//				holdJVec.reserve(distances.size());
//				lows.reserve(distances.size());
//				highs.reserve(distances.size());
//				for (int i = 0; i < pointCandidates; i++)
//				{
//					if (distances[i] < fwhmRfi)
//					{
//						holdI = possibles[2 * i];
//						holdJ = possibles[2 * i + 1];
//						dec.push_back(scans[holdI].getDec(holdJ));
//						ra.push_back(scans[holdI].getRa(holdJ));
//						holdIVec.push_back(holdI);
//						holdJVec.push_back(holdJ);
//						lows.push_back(999999);
//						highs.push_back(-999999);
//						if (trueCount != i)
//						{
//							distances[trueCount] = distances[i];
//						}
//						trueCount++;
//					}
//				}
//				distances.erase(distances.begin() + trueCount, distances.end());
//
//				for (int i = 0; i < holdIVec.size(); i++)
//				{
//					rfiCountHold = 0;
//					normalLocalWeight = 0;
//					normalLocalWeightCounter = 0;
//					holdI = holdIVec[i];
//					holdJ = holdJVec[i];
//					if (checks[i])
//					{
//						distSq = 0.0;
//						wDistSq = 0.0;
//						pointCount = 0.0;
//						wPointCount = 0.0;
//						distance = distances[i];
//						shapeParam = pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2);
//						fitHold = f * shapeParam + g;
//
//						if (correlatedWeightMap)
//						{
//							//We determine the first correction factor and the correlation scale for photometry here.
//							rejectPoints = true;
//							correlatedScatter = scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()];
//
//							//We assume that each grid point is above sqrt(2)*scatter.
//
//							for (int j = 0; j < holdIVec.size(); j++)
//							{
//								rfiDistance = Tools::getGCDistance(scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), scans[holdIVec[j]].getDec(holdJVec[j]), scans[holdIVec[j]].getRa(holdJVec[j]), partSetProcSSS.centerDecDeg)*toDeg;
//								if (rfiDistance < fwhmRfi*toDeg)
//								{
//									if (scans[holdIVec[j]].getFlux(holdJVec[j]) > sqrt(2.0)*correlatedScatter) {
//										//Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are above the 2d-scatter.
//										rfiCountHold += 1;
//									}
//								}
//								if (checks[j]) //Correlation scale here.
//								{
//									pointDistance = Tools::getGCDistance(scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), scans[holdIVec[j]].getDec(holdJVec[j]), scans[holdIVec[j]].getRa(holdJVec[j]), partSetProcSSS.centerDecDeg)*toDeg;
//									if (pointDistance < fwhmRfi*toDeg && pointDistance != 0)
//									{
//										if (scans[holdIVec[j]].getFlux(holdJVec[j]) > sqrt(2.0)*correlatedScatter)
//										{
//											pointCount += 1.0;
//											distSq += pow(pointDistance, 2);
//
//											wPointCount += 1.0*scans[holdIVec[j]].getDataDumps(holdJVec[j]);
//											wDistSq += pow(pointDistance, 2)*scans[holdIVec[j]].getDataDumps(holdJVec[j]);
//										}
//									}
//
//								}
//
//								distanceTemp = distances[j];
//								shapeParamTemp = pow(cos((M_PI*distanceTemp) / (fwhmRfi*2.0)), 2);
//								normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2); //sum of proximity weight of this point minus all other points used to make the local model squared.
//								normalLocalWeightCounter++;
//							}
//						}
//						else
//						{
//							for (int j = 0; j < holdIVec.size(); j++)
//							{
//								distanceTemp = distances[j];
//								shapeParamTemp = pow(cos((M_PI*distanceTemp) / (fwhmRfi*2.0)), 2);
//								normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2); //sum of proximity weight of this point minus all other points used to make the local model squared.
//								normalLocalWeightCounter++;
//							}
//							rfiCountHold = 1;
//						}
//
//						if (rfiCountHold != rfiCountHold || rfiCountHold == 0) {
//							rfiCountHold = 1;
//						}
//						if (pointsInRange == 1) {
//							rfiWeights[holdI][holdJ].push_back(inRange[2]);
//						}
//						else if (pointsInRange == 2) {
//							rfiWeights[holdI][holdJ].push_back(inRange[2] + inRange[3 * 1 + 2]);
//						}
//						else {
//							rfiWeights[holdI][holdJ].push_back((localModelDumpSum / (1.0 + 2.0 * pow((shapeParam - localModelShapeSum), 2) / (normalLocalWeight / normalLocalWeightCounter))));
//						}
//
//						rfiValues[holdI][holdJ].push_back(fitHold);
//						//rfiWeights[holdI][holdJ].push_back((localModelDumpSum / (1.0 + 2.0 * pow((shapeParam - localModelShapeSum), 2) / (normalLocalWeight / normalLocalWeightCounter))));
//						//rfiWeights[holdI][holdJ].push_back((localModelDumpSum / (1.0 + 2.0 * (shapeParam*shapeParam) / (localModelShapeSum))));
//						rfiCounts[holdI][holdJ].push_back(rfiCountHold);
//						rfiDists[holdI][holdJ].push_back(wDistSq);
//						rfiNs[holdI][holdJ].push_back(wPointCount);
//					}
//				}
//			}
//		}
//	}
//}
//void Processor::extraLocalModels(double decPoint, double raPoint, std::vector<std::vector<std::vector<double> > > &rfiValues, std::vector<std::vector<std::vector<double> > > &rfiWeights, std::vector<std::vector<std::vector<double> > >&rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > & rfiNs, double rfiScale, std::vector<MapTypes> mapHolder)
//{
//	int rfiCountHold = 0;
//	double comDec, comRa;
//	bool extraModels = false;
//	double standardGap = 999999;
//	std::vector<double> originalGaps;
//	originalGaps.resize(standardGaps.size());
//
//	//correlatedWeightMap variables
//	bool rejectPoints;
//	double stdDev, distSq, wDistSq, wPointCount, pointCount, rejLimit;
//	double rfiDistance, pointDistance;
//
//	for (int i = 0; i < standardGaps.size(); i++)
//	{
//		originalGaps[i] = standardGaps[i];
//		if (mapHolder[i] == DAISY)
//		{
//			standardGaps[i] = (standardGaps[i] / 2)*scans.size();
//		}
//		if (standardGaps[i] < standardGap)
//		{
//			standardGap = standardGaps[i];
//		}
//	}
//
//	standardGaps = originalGaps;
//
//
//	int extraSteps = std::floor((psfFWHM / standardGap)*pow((rfiScale), 2));
//
//	if (extraSteps % 2 != 0)
//	{
//		extraSteps += 1;
//	}
//
//	if ((rfiScale < 1.0) && (rfiScale > pow((standardGap / psfFWHM), 0.5)))
//	{
//		for (int i = -extraSteps; i <= extraSteps; i++)
//		{
//			for (int j = -extraSteps; j <= extraSteps; j++)
//			{
//				int  holdI, holdJ, pointCandidates = 0;
//				double avgSigma = 0, distance, toRad = M_PI / 180.0;
//				double correlatedScatter, rfiDistance, pointDistance, fluxVal;
//				std::vector<int> checks, possibles;
//				std::vector<double> inRange, rfiFitSolution, results;
//
//				double deltaX = i * (standardGap)*(pow((1 / rfiScale), 2) - 1); //In degrees
//				double deltaY = j * (standardGap)*(pow((1 / rfiScale), 2) - 1);
//
//				double sampleWidth = psfFWHM * (1 - pow(rfiScale, 2));
//
//				distance = Tools::getGCDistance((decPoint + deltaX), (raPoint + deltaY), decPoint, raPoint, partSetProcSSS.centerDecDeg)*toDeg;
//
//				if (distance > sampleWidth * 2)
//				{
//					continue;
//				}
//
//				pointCandidates = findPossSSS((decPoint + deltaX), (raPoint + deltaY), fwhmRfi*toDeg, possibles);
//				inRange.reserve(pointCandidates * 6);
//				checks.resize(pointCandidates, true);
//
//				for (int i = 0; i < pointCandidates; i++)
//				{
//					holdI = possibles[2 * i];
//					holdJ = possibles[2 * i + 1];
//					distance = Tools::getGCDistance((decPoint + deltaX), (raPoint + deltaY), scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg);
//
//					inRange.push_back(scans[holdI].getFlux(holdJ));
//					inRange.push_back(pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2));
//					inRange.push_back(scans[holdI].getDataDumps(holdJ));
//					inRange.push_back(scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()]);
//					inRange.push_back(holdI);
//					inRange.push_back(holdJ);
//
//					avgSigma += scatter2d[scans[holdI].getSurveyNumber()][scans[holdI].getScanNumberInSurvey()];
//				}
//
//				avgSigma = (avgSigma / (double)pointCandidates);
//				rfiRemovePoints(avgSigma, checks, inRange);
//				rfiFitSolution = rfiFit(checks, inRange);
//				results = rfiCollectResults(999999, -999999, rfiFitSolution, checks, inRange);//We assume these data are above the correlated scatter value
//			}
//		}
//	}
//}

//CLEANNED RFI FUNCTIONS
/*
void Processor::autoCentroid(std::vector<std::vector<std::vector<double> > > rfiValues)
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
std::vector<double> Processor::determineCenters(int pointCount, std::vector<double>& inRange, double largestGap, double i_0, double j_0)
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

coef = Tools::matrixSolver(columnCount, A, b);

double x = ((coef[2] * coef[5]) - (2 * coef[1] * coef[4])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));
double y = ((coef[1] * coef[5]) - (2 * coef[2] * coef[3])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));

coordinates.push_back(x);
coordinates.push_back(y);

return coordinates;
}
*/
/*
std::vector<double> rfiFit(std::vector<int> &checks, std::vector<double> &inRange)
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
void Processor::rfiRemovePoints(double avgSigma, std::vector<int> &checks, std::vector<double> &inRange)
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
void Processor::rfiBuildLocalModels(std::vector<std::vector<std::vector<double> > >& rfiWeights, std::vector<std::vector<std::vector<double> > >& rfiValues, std::vector<std::vector<std::vector<double> > >& rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > &rfiNs)
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

//rfiRemovalMulti(decHold, raHold, fluxVal, correlatedScatter);
futureFiller[j] = std::async(std::launch::async, &Processor::rfiRemovalMulti, this, decHold, raHold, fluxVal, correlatedScatter);
}

rfiSetFinalValues(futureFiller, rfiWeights, rfiValues, rfiCounts, rfiDists, rfiNs);
}
}
void Processor::rfiBuildExtraModels(double decPoint, double raPoint, double standardGap, std::vector<std::vector<std::vector<double> > >& rfiWeights, std::vector<std::vector<std::vector<double> > >& rfiValues, std::vector<std::vector<std::vector<double> > >& rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > &rfiNs)
{
int extraSteps;
double deltaX, deltaY, simpleWidth, distance;
double decHold, raHold, fluxVal, correlatedScatter;
double sampleWidth = psfFWHM * (1 - pow(rfiScale, 2));
std::vector<std::future<std::vector<double>>> futureFiller;

extraSteps = std::floor((psfFWHM / standardGap)*pow((rfiScale), 2));
futureFiller.reserve(extraSteps * 2 + 1);
if (extraSteps % 2 != 0)
{
extraSteps += 1;
}

for (int i = -extraSteps; i <= extraSteps; i++)
{
for (int j = -extraSteps; j <= extraSteps; j++)
{
deltaX = i * (standardGap)*(pow((1 / rfiScale), 2) - 1); //In degrees
deltaY = j * (standardGap)*(pow((1 / rfiScale), 2) - 1);
decHold = decPoint + deltaX;
raHold = raPoint + deltaY;
fluxVal = 999999;//WE ASSUME THESE POINTS ARE HIGHER THAN THE SCATTER
correlatedScatter = -999999;
distance = Tools::getGCDistance(decHold, raHold, decPoint, raPoint, partSetProcSSS.centerDecDeg)*toDeg;

if (distance > sampleWidth * 2)
{
continue;
}
futureFiller.push_back(std::async(std::launch::async, &Processor::rfiRemovalMulti, this, decHold, raHold, fluxVal, correlatedScatter));
}
rfiSetFinalValues(futureFiller, rfiWeights, rfiValues, rfiCounts, rfiDists, rfiNs);
}

}
void Processor::rfiSetFinalValues(std::vector<std::future<std::vector<double>>> & futureVec, std::vector<std::vector<std::vector<double> > >& rfiWeights, std::vector<std::vector<std::vector<double> > >& rfiValues, std::vector<std::vector<std::vector<double> > >& rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > &rfiNs)
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
std::vector<double> Processor::rfiRemovalMulti(double dec, double ra, double fluxVal, double correlatedScatter)
{
int  holdI, holdJ, pointCandidates = 0;
double avgSigma = 0, distance, toRad = M_PI / 180.0;
std::vector<int> checks, possibles;
std::vector<double> inRange, rfiFitSolution, results;
pointCandidates = findPossSSS(dec, ra, fwhmRfi*toDeg, possibles);
inRange.reserve(pointCandidates * 6);
checks.resize(pointCandidates, true);

//COLLECT POSSIBLES AND STORE INFORMATION ABOUT PROXIMITY, SCATTER, AND FLUX
for (int i = 0; i < pointCandidates; i++)
{
holdI = possibles[2 * i];
holdJ = possibles[2 * i + 1];

distance = Tools::getGCDistance(dec, ra, scans[holdI].getDec(holdJ), scans[holdI].getRa(holdJ), partSetProcSSS.centerDecDeg);

inRange.push_back(scans[holdI].getFlux(holdJ));
inRange.push_back(pow(cos((M_PI*distance) / (fwhmRfi*2.0)), 2));
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
std::vector<double> Processor::rfiCollectResults(double fluxVal, double correlatedScatter, std::vector<double> &rfiFitSolution, std::vector<int> & checks, std::vector<double> &inRange)
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

double pw; //proximity weight -- cos(M_PI*dist/2*rfiScale)^2
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
//localModelShapeSum += inRange[6 * i + 1];
//localModelDumpSum += inRange[6 * i + 2];

checksTrue++;
}
}

pwAvg = pwSum / (double)checksTrue;

//localModelShapeSum = localModelShapeSum / (double)checksTrue;


results.reserve((inRange.size() / 6) * 7);
for (int i = 0; i < inRange.size() / 6; i++)
{
//rfiCountHold = 0;
//normalLocalWeightCounter = 0;
//normalLocalWeight = 0;

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

//shapeParamTemp = inRange[6 * j + 1];
//normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2); //sum of proximity weight of this point minus all other points used to make the local model squared.
//normalLocalWeight += pow(shapeParamTemp - localModelShapeSum, 2);

//normalLocalWeight += pow(pwTemp - pwAvg, 2);
//normalLocalWeightCounter++;

if (correlatedWeightMap)//We determine the first correction factor and the correlation scale for photometry here.
{
rfiCountHold = 0.0;
countAbove = (centerLMFlux > sqrt(2.0)*correlatedScatter) ? true : false;
fluxAbove = (jLMFlux > sqrt(2.0)*correlatedScatter) ? true : false;

dd = scans[jHoldI].getDataDumps(jHoldJ);
pointDistance = Tools::getGCDistance(iDec, iRa, jDec, jRa, partSetProcSSS.centerDecDeg)*toDeg;//CHECK THAT fwhmRfi*toDeg is right units

if (checks[j] && (pointDistance < fwhmRfi*toDeg) && ((countAbove && fluxAbove) || (!countAbove && !fluxAbove))) {//Includes all points that are within a rfi-beam of LMV && within a rfi-beam of GMV && are (Case 1: above OR Case 2: below) the 2d-scatter.
rfiCountHold += 1;
if ((pointDistance != 0))
{
//pointCount += 1.0;
//distSq += pow(pointDistance, 2);
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
//wResult = (localModelDumpSum / (1.0 + 2.0 * pow((shapeParam - localModelShapeSum), 2) / (normalLocalWeight / normalLocalWeightCounter)));
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
std::vector<std::vector<double> > Processor::rfiBuildGlobalModels(std::vector<std::vector<std::vector<double> > >& rfiWeights, std::vector<std::vector<std::vector<double> > >& rfiValues, std::vector<std::vector<std::vector<double> > >& rfiCounts, std::vector<std::vector<std::vector<double> > > & rfiDists, std::vector<std::vector<std::vector<double> > > &rfiNs)
{
double wSum, dSum, nSum;
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
futureFiller[j] = std::async(std::launch::async, spawnNewRCRThread, wVec, xVec);

}
else if (rfiValues[i][j].size() == 1)
{
rfiSubtracted[i][j] = rfiValues[i][j][0];
scans[i].setGMWeight(j, rfiWeights[i][j][0] / rfiCounts[i][j][0]);
scans[i].setSSSCorrelation(j, 2.0*factor_corr*sqrt(rfiDists[i][j][0] / std::max(rfiNs[i][j][0], 1.0)));
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
scans[i].setGMWeight(j, wSum);
scans[i].setSSSCorrelation(j, 2.0*factor_corr*sqrt(dSum / std::max(nSum, 1.0)));
}
}
}
return rfiSubtracted;
}
*/


//theta gap
/*
std::vector<int> Processor::quadrantSort(double dec, double ra, std::vector<int> &inRange, Quadrant quadrant)
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
void Processor::safeguardEdgeThetaGaps(Composite &composite, double weightScale)
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
void Processor::determineCircleParams(double dec, double ra, Quadrant quadrant, std::vector<int> &inRange, std::vector<double> &circleParams)
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
			std::cout << "NON-SENSIBLE RADIUS CALCULATED" << "\n";
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
double Processor::maxGapQuadrant(std::vector<double> &circleParams)
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
double Processor::circleThetaGap(int i_0, int j_0, bool LSSOn)
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
*/

Processor::~Processor()
{

}