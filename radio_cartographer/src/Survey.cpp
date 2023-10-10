#include "Survey.h"
#include "Tools.h"
#include "GBParser.h"
#include "FourtyParser.h"
#include "OutputFile.h"
#include "RCR.h"
#include "Debugger.h"
#include <time.h>
#include <math.h>
#include <sstream>
#include <iomanip> 

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

LinearModel::LinearModel(std::vector<double> x, std::vector<double> y)
{
	this->x = x;
	this->y = y;
	this->w.resize(x.size(), 1.0);
}
LinearModel::LinearModel(std::vector<double> w, std::vector<double> x, std::vector<double> y)
{
	this->x = x;
	this->y = y;
	this->w = w;
}
std::vector<double> LinearModel::getErrors(std::vector<double> line)
{
	m = line[0], b = line[1];
	//std::cout << m << "\t" << b << "\t" << xBar << "\n";
	double lineY;
	std::vector<double> toRet;
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			lineY = m*(x[i] - xBar) + b;
			toRet.push_back(y[i] - lineY);
		}
	}
	return toRet;
}
void LinearModel::buildModelSpace()
{
	parameterSpace.clear();
	weightSpace.clear();
	double mHold, weightProduct;
	getAverage(); //ADD WEIGHTED AVERAGE
	std::vector<double> m, b, mW, bW, mbW;
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			for (int j = i + 1; j < y.size(); j++)
			{
				if (flags[j] && x[j] != x[i])
				{
					mHold = (y[j] - y[i]) / (x[j] - x[i]);
					m.push_back(mHold);
					b.push_back(y[j] - mHold * (x[j] - xBar));
					weightProduct = w[i] * w[j];
					mW.push_back((.5 * pow(x[j] - x[i], 2.0)) * weightProduct);
					bW.push_back((pow(x[j] - x[i], 2.0) / (pow(x[i] - xBar, 2.0) + pow(x[j] - xBar, 2.0))) * weightProduct);
				}
			}
		}
	}
	parameterSpace.push_back(m);
	parameterSpace.push_back(b);
	weightSpace.push_back(mW);
	weightSpace.push_back(bW);
}
void LinearModel::getAverage()
{
	double dSum = 0, wSum = 0;
	double trueCount = 0;
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			dSum += w[i] * x[i];
			wSum += w[i];
			//trueCount += 1.;
		}
	}
	xBar = dSum / wSum;
	//std::cout << xBar << " XBAR\n";
}
void LinearModel::printData()
{
	std::ofstream outputData;
	outputData.open("15XYDataAfterWegithed.txt");
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			outputData << x[i] << "\t" << y[i] << "\n";
		}
	}
	outputData.close();
}
std::vector<double> LinearModel::regression()
{
	double wxSum = 0, wySum = 0, wSum = 0, wxBar, wyBar, top = 0, bottom = 0;
	std::vector<double> ret;
	ret.resize(5);
	getAverage(); //ADD WEIGHTED AVERAGE
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			wySum += w[i] * y[i];
			wSum += w[i];
		}
	}
	wyBar = wySum / wSum;
	for (int i = 0; i < x.size(); i++)
	{
		if (flags[i])
		{
			top += w[i] * (x[i] - xBar) * (y[i] - wyBar);
			bottom += w[i] * (x[i] - xBar) * (x[i] - xBar);
		}
	}
	m = top / bottom;
	b = wyBar;
	std::vector<double> toRet;
	toRet.push_back(m);
	toRet.push_back(b);
	//std::cout << "M: " << m << " B: " << b << "xBar: " << xBar << "\n";
	return toRet;
}
LinearModel::~LinearModel()
{

}


//constructors
Survey::Survey()
{

}
Survey::Survey(SurveyParameters &sParams, Input input)
{
	setParams(sParams);
	setTelescopeParams(sParams, input);
	determineInputFile(input.filename);
	initializeData(input);
}

void Survey::setTelescopeParams(SurveyParameters &params, Input input) {

	if (input.telescope == "GreenBank-20" || input.telescope == "NRAO20" || input.telescope ==  "NRAO_GBT") {

		if (input.telescope != "NRAO_GBT") this->telescope = TWENTY_METER;
		else this->telescope = GBT;
		
		double telescopeDiameter = this->telescope == GBT ? 105.0 : 20.0;

		params.frequency = input.observedFrequencies[0];

		this->telescopeFrequency = params.frequency;
		this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * telescopeDiameter * M_PI);

		setMappingCoordinate(input.coordinate);

		if (input.mapPattern == "ralongmap")
		{
			this->mapType = RASTER;
			this->scansInRa = true;
		}
		else if (input.mapPattern == "declatmap")
		{
			this->mapType = RASTER;
			this->scansInRa = false;
		}
		else
		{
			this->mapType = DAISY;
			this->scansInRa = false;
		}
	}
	else if (input.telescope == "MightyForty") {
		this->telescope = FOURTY_FOOT;
	}
	else {
		throw "Unsupported telescope: " + input.telescope + "\n";
	}
}

Survey::Survey(SurveyParameters &sParams, SpectralParameters cParams, std::string filename)
{
	setParams(sParams);
	determineInputFile(filename);

	if (!ASCII)
	{
		std::vector<std::vector<double>> sdfitsData;
		PreProcessor sdfits = PreProcessor();
		sdfitsData = sdfits.sdfitsReader(cParams);
		std::cout << "Out of PreProcessor\n";
		setSdfitsParams(sParams, sdfits);
		initializeData(sdfitsData);
		std::cout << "Out of Initialize\n";
	}
	/*else
	{
		if (telescope == TWENTY_METER)
		{
			GBParser gb = GBParser(filename);
			gbData = gb.parseFile();
			setTwentyParams(gb);
			initializeData(gbData);
		}
		else
		{
			FourtyParser gb = FourtyParser(filename);
			gbData = gb.parseFile();
			setFortyParams(gb);
			initializeData(gbData);
		}
	}*/
}
Survey::Survey(SurveyParameters params, std::string fileName, std::string fileNameB)
{
	/*setParams(params);
	this->channel = COMPOSITE;

	if (telescope == TWENTY_METER)
	{
		GBParser gb = GBParser(fileName);
		gbData = gb.parseFile();
		setTwentyParams(gb);
		initializeData(gbData);
	}
	else
	{
		FourtyParser gb = FourtyParser(fileName, fileNameB);
		gbData = gb.parseFile();
		setFortyParams(gb);
		initializeData(gbData);
	}*/
}

//automated preprocessing
void Survey::setParams(SurveyParameters params)
{
	this->calMethod = params.calMethod;
	this->channel = params.channel;
	this->tracking = params.tracking;
	this->t_int = params.forcedTS;
	this->trimSize = params.trimSize;
	this->pCoordinate = params.pCoordinate;
	this->debugging = params.debug;
}
void Survey::setTwentyParams(GBParser &gb)
{
	// Set parameters specific to the 20m

	this->MJD = gb.getParamMJD();
	setMappingCoordinate(gb.getParamCoordinate());

	std::string scanTypeStr = gb.getParamScanType();
	std::string mapTypeStr = gb.getParamMapType();
	double obsTemp = gb.getParamFrequency();

	this->telescopeFrequency = obsTemp*pow(10, 9);
	this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * 20.0 * M_PI);

	if (mapTypeStr == "DAISY")
	{
		this->mapType = DAISY;
	}
	else
	{
		this->mapType = RASTER;
	}
	if (scanTypeStr == "RA")
	{
		this->scansInRa = true;
	}
	else
	{
		this->scansInRa = false;
	}
	if (mCoordinate == "")
	{
		mCoordinate = "equatorial";
		Debugger::print("Warn", "No coordinate system found in header.. forced to equatorial");
	}

}
void Survey::setSdfitsParams(SurveyParameters &params, PreProcessor sdfits)
{
	std::string mapPattern = sdfits.getMapPattern();
	std::string telescopeName = sdfits.getTelescope();

	if (telescopeName == "GreenBank-20" || telescopeName == "NRAO20")
	{
		// Used for RFIMax
		params.tele = TWENTY_METER;
		params.frequency = sdfits.getFrequency() / 1000.0;

		this->MJD = sdfits.getMJDate();
		this->telescope = TWENTY_METER;
		this->telescopeFrequency = sdfits.getFrequency() * pow(10, 6);
		this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * 20.0 * M_PI);

		setMappingCoordinate(sdfits.getCoordinateSystem());

		if (mapPattern == "ralongmap")
		{
			this->mapType = RASTER;
			this->scansInRa = true;
		}
		else if (mapPattern == "declatmap")
		{
			this->mapType = RASTER;
			this->scansInRa = false;
		}
		else
		{
			this->mapType = DAISY;
			this->scansInRa = false;
		}
	}
	else if (telescopeName == "NRAO_GBT")
	{
		params.tele = GBT;
		params.frequency = sdfits.getFrequency() / pow(10, 9);

		this->telescope = GBT;
		this->MJD = sdfits.getMJDate();
		this->telescopeFrequency = sdfits.getFrequency();
		this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * 105.0 * M_PI);

		setMappingCoordinate(sdfits.getCoordinateSystem());

		if (mapPattern == "RALongMap")
		{
			this->mapType = RASTER;
			this->scansInRa = true;
		}
		else if (mapPattern == "DecLatMap")
		{
			this->mapType = RASTER;
			this->scansInRa = true;
		}
		else
		{
			this->mapType = DAISY;
			this->scansInRa = false;
		}
	}
	else // forty foot
	{
		params.tele = FOURTY_FOOT;
		params.frequency = 1.405;

		setMappingCoordinate("equatorial");

		this->scansInRa = false;
		this->mapType = NODDING;
		this->telescope = FOURTY_FOOT;
		this->telescopeFrequency = 1.405 * pow(10, 9);
		this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * 12.192 * M_PI);
	}
}
void Survey::setFortyParams(FourtyParser &gb)
{
	this->scansInRa = false;
	this->mapType = NODDING;
	this->telescopeFrequency = 1.405*pow(10, 9);
	this->psfFWHM = 1.22*299792458.0*180.0 / (telescopeFrequency * 12.192 * M_PI);
}
void Survey::determineInputFile(std::string fileName)
{
	if (fileName.find(".txt") != fileName.npos || fileName.find(".md2") != fileName.npos)
	{
		ASCII = true;
	}
	else if (fileName.find(".fits") != fileName.npos)
	{
		ASCII = false;
	}
	else
	{
		// throw error;
	}
}
void Survey::switchChannels(Channel newChannel)
{
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].switchChannels(newChannel);
	}
}

//file processors
void Survey::initializeData(Input input) {
	std::vector<std::vector<double>> data;
	data.resize(11);

	data[0] = input.times;
	data[1] = input.ras;
	data[2] = input.decs;
	data[3] = input.azimuths;
	data[4] = input.elevations;
	data[5] = input.lContinuum;
	data[6] = input.rContinuum;
	data[7] = input.dataDumps;
	data[8] = input.calibrations;
	data[9] = input.sweepIndices;
	data[10] = input.valids;

	if (telescope == TWENTY_METER || telescope == GBT)
	{
		dataProc(data);
	}
	else
	{
		dataProc40(data);
	}
}

void Survey::initializeData(std::vector<std::vector<double> > data)
{
	if (telescope == TWENTY_METER || telescope == GBT)
	{
		dataProc(data);
	}
	else
	{
		dataProc40(data);
	}
	Debugger::print("Info", "data loaded");
}
void Survey::dataProc(std::vector<std::vector<double> > &data)
{
	if (data.size() == 11)
	{
		formatData11(data);
	}
	else
	{
		//throw error for unrecognized data format;
	}

	//WHEN YOU TAKE OBSERVATIONS WITH SKYNET IN GALACTIC COORDINATES, THE
	//TEXT FILE WILL STILL BE IN RA AND DEC. IN ORDER TO GET THE SQUARE
	//GRID, YOU NEED TO CONVERT THE COORDINATES TO GALACTIC INITIALLY.
	if (mCoordinate == "galactic")
	{
		convertToGalacticInitial();
		zeroCrossCheck();
	}

	int offset = 1;
	std::vector<double> dubFiller;
	
	Output output;
	output.printCalibrationHeader();
	
	if (telescope != GBT)
	{
		gainCalibration(LEFT, channel);
		gainCalibration(RIGHT, channel);
	}

	for (int i = 0; i < fluxL.size(); i++)
	{
		fluxComp.push_back(dubFiller);
		for (int j = 0; j < fluxL[i].size(); j++)
		{
			fluxComp[i].push_back(0.5*(fluxL[i][j] + fluxR[i][j]));
		}
	}	
	if (mapType == DAISY)
	{
		daisyPrelim();
	}
	janskyCalibration(1.0, LEFT);
	janskyCalibration(1.0, RIGHT);
	janskyCalibration(1.0, COMPOSITE);

	for (int i = offset; i < times.size() - offset; i++)
	{
		if (times[i].size() > 1) {
			this->scans.push_back(Scan(times[i], decs[i], ras[i], elevations[i], dataDumps[i], fluxL[i], fluxR[i], fluxComp[i]));
		}
		else {
			std::cout << "[ WARN ] No valid data found for scan " + i << std::endl;
		}
	}

	Debugger::print("Info", "scans loaded");

	if (mapType == DAISY)
	{
		findCenters();
		Debugger::print("Info", "centers loaded");
	}

	for (int i = 0; i < scans.size(); i++)
	{
		if (mapType == DAISY)
		{
			scans[i].undoCosTransform(0, partSetProcSSS.medianDec, partSetProcSSS.medianRa);// in daisyPrelim, to divide scans we transform onto the center of the grid and calculated distances and times to break apart scans.  We want to undo this for generality
			scans[i].updateAngDistTemp(partSetProcSSS.centerDecDeg);
		}
		scans[i].setScanNumberInSurvey(i);
	}



}
void Survey::dataProc40(std::vector<std::vector<double> > &data)
{
	bool lineCheck = true;
	int i;
	std::vector<double> dubFiller;

	//WHEN YOU TAKE OBSERVATIONS WITH SKYNET IN GALACTIC COORDINATES, THE
	//TEXT FILE WILL STILL BE IN RA AND DEC. IN ORDER TO GET THE SQUARE
	//GRID, YOU NEED TO CONVERT THE COORDINATES TO GALACTIC INITIALLY.
	if (mCoordinate == "galactic")
	{
		convertToGalacticInitial();
	}

	times.push_back(dubFiller);
	ras.push_back(dubFiller);
	decs.push_back(dubFiller);
	azimuths.push_back(dubFiller);
	elevations.push_back(dubFiller);
	fluxL.push_back(dubFiller);
	fluxR.push_back(dubFiller);
	dataDumps.push_back(dubFiller);
	calibrationFlags.push_back(dubFiller);

	scanCount = 0;
	i = 0;
	while (data[9][i] == -1)
	{
		times[scanCount].push_back(data[0][i]);
		ras[scanCount].push_back(15.0*data[1][i]);
		decs[scanCount].push_back(data[2][i]);
		azimuths[scanCount].push_back(data[3][i]);
		elevations[scanCount].push_back(data[4][i]);
		fluxL[scanCount].push_back(data[5][i]);
		fluxR[scanCount].push_back(data[6][i]);
		dataDumps[scanCount].push_back(data[7][i]);
		calibrationFlags[scanCount].push_back(data[8][i]);
		i++;
	}
	times.push_back(dubFiller);
	ras.push_back(dubFiller);
	decs.push_back(dubFiller);
	azimuths.push_back(dubFiller);
	elevations.push_back(dubFiller);
	fluxL.push_back(dubFiller);
	fluxR.push_back(dubFiller);
	dataDumps.push_back(dubFiller);
	calibrationFlags.push_back(dubFiller);
	scanCount = 1;
	while (i < data[0].size())
	{
		if (data[9][i] == scanCount - 1)
		{
			times[scanCount].push_back(data[0][i]);
			ras[scanCount].push_back(15.0*data[1][i]);
			decs[scanCount].push_back(data[2][i]);
			azimuths[scanCount].push_back(data[3][i]);
			elevations[scanCount].push_back(data[4][i]);
			fluxL[scanCount].push_back(data[5][i]);
			fluxR[scanCount].push_back(data[6][i]);
			dataDumps[scanCount].push_back(data[7][i]);
			calibrationFlags[scanCount].push_back(data[8][i]);
		}
		else if (data[9][i] == scanCount - 1 && data[10][i] == 0 && mapType != DAISY)
		{

		}
		else if (data[9][i] == scanCount)
		{
			scanCount++;
			times.push_back(dubFiller);
			ras.push_back(dubFiller);
			decs.push_back(dubFiller);
			azimuths.push_back(dubFiller);
			elevations.push_back(dubFiller);
			fluxL.push_back(dubFiller);
			fluxR.push_back(dubFiller);
			dataDumps.push_back(dubFiller);
			calibrationFlags.push_back(dubFiller);

		}
		else if (data[9][i] == -1)
		{
			if (lineCheck == true)
			{
				scanCount++;
				times.push_back(dubFiller);
				ras.push_back(dubFiller);
				decs.push_back(dubFiller);
				azimuths.push_back(dubFiller);
				elevations.push_back(dubFiller);
				fluxL.push_back(dubFiller);
				fluxR.push_back(dubFiller);
				dataDumps.push_back(dubFiller);
				calibrationFlags.push_back(dubFiller);
				lineCheck = false;
			}
			times[scanCount].push_back(data[0][i]);
			ras[scanCount].push_back(15.0*data[1][i]);
			decs[scanCount].push_back(data[2][i]);
			azimuths[scanCount].push_back(data[3][i]);
			elevations[scanCount].push_back(data[4][i]);
			fluxL[scanCount].push_back(data[5][i]);
			fluxR[scanCount].push_back(data[6][i]);
			dataDumps[scanCount].push_back(data[7][i]);
			calibrationFlags[scanCount].push_back(data[8][i]);

		}
		i++;
	}

	correctDec40(decs, times, 10);

	int offset = 1;
	double minDec, maxDec, decTemp;
	minDec = 999999;
	maxDec = -999999;
	for (int i = offset; i < decs.size() - offset; i++)
	{
		for (int j = 0; j < decs[i].size(); j++)
		{
			decTemp = decs[i][j];
			if (decTemp > maxDec)
			{
				maxDec = decTemp;
			}
			if (decTemp < minDec)
			{
				minDec = decTemp;
			}
		}
	}

	gainCalibration(LEFT,channel);
	gainCalibration(RIGHT,channel);
	for (int i = 0; i < fluxL.size(); i++)
	{
		fluxComp.push_back(dubFiller);
		for (int j = 0; j < fluxL[i].size(); j++)
		{
			fluxComp[i].push_back(.5*(fluxL[i][j] + fluxR[i][j]));

		}
	}
	janskyCalibration(1.0, LEFT);
	janskyCalibration(1.0, RIGHT);
	janskyCalibration(1.0, COMPOSITE);

	for (int i = offset; i < times.size() - offset; i++)
	{
		this->scans.push_back(Scan(times[i], decs[i], ras[i], elevations[i], dataDumps[i], fluxL[i], fluxR[i], fluxComp[i]));
	}

	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setScanNumberInSurvey(i);
	}
	Debugger::print("Info", "scans loaded");
}
void Survey::correctDec40(std::vector<std::vector<double> > &decs, std::vector<std::vector<double> > &times, int N)
{
	RCR rcr = RCR(LS_MODE_DL);
	rcr.setMuType(VALUE);
	std::vector<double> fullLineDecs, fullLineTimes;
	std::vector<double> replacementDecs;
	std::vector<double> weights;
	std::vector<double> coef;
	std::vector<bool> checks;

	LinearModel model = LinearModel(weights, fullLineTimes, fullLineDecs);
	rcr.setParametricModel(model);


	checks.resize(2 * N + 1, true);
	weights.resize(2 * N + 1, 1.0);

	int counter = 0;
	int centerIndex;
	int bound;

	//FOR ALL SCANS
	for (int i = 0; i < decs.size(); i++)
	{
		replacementDecs.resize(decs[i].size());

		//CENTER ON EACH DEC
		for (int j = 0; j < decs[i].size(); j++)
		{

			//COLLECT 5 POINTS BEFORE AND AFTER THE CENTER INDEX
			counter = 0;
			centerIndex = j;
			if (j > N && j < decs[i].size() - N)
			{
				fullLineDecs.resize(2 * N + 1);
				fullLineTimes.resize(2 * N + 1);
				for (int k = -N; k <= N; k++)
				{
					fullLineDecs[counter] = decs[i][j + k];
					fullLineTimes[counter] = times[i][j + k];
					counter++;
				}
			}
			else if (j <= N)
			{
				fullLineDecs.resize(N + 1 + j);
				fullLineTimes.resize(N + 1 + j);
				for (int k = -j; k <= N; k++)
				{
					fullLineDecs[counter] = decs[i][j + k];
					fullLineTimes[counter] = times[i][j + k];
					counter++;
				}
			}
			else
			{
				fullLineDecs.resize(N + 1 + decs[i].size() - 1 - j);
				fullLineTimes.resize(N + 1 + decs[i].size() - 1 - j);
				bound = (decs[i].size() - 1 - j);
				for (int k = -N; k <= bound; k++)
				{
					fullLineDecs[counter] = decs[i][j + k];
					fullLineTimes[counter] = times[i][j + k];
					counter++;
				}
			}
			weights.resize(fullLineDecs.size(), 1.0);
			model = LinearModel(weights, fullLineTimes, fullLineDecs);
			rcr.performBulkRejection(weights, fullLineDecs);


			//coef = Tools::regressionPivot(0, 0, checks, weights, fullLineTimes, fullLineDecs);
			fullLineDecs.clear();
			fullLineTimes.clear();
			replacementDecs[j] = model.m * (times[i][j]-model.xBar) + model.b;
		}
		decs[i] = replacementDecs;
		replacementDecs.clear();
	}
}
void Survey::formatData11(std::vector<std::vector<double> > &data)
{
	int i = 0;
	bool lineCheck = true;
	std::vector<double> dubFiller;

	if (mapType != DAISY)
	{
		//HOLD SPACE FOR INITIAL CALIBRATION

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);

		//HOLD SPACE FOR THE FIRST SCAN

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);

		//CALIBRATION DATA DURING TRANSITION

		while (i < data[0].size())
		{
			//STORE ALL INITIAL CALIBRATION DATA AND ZERO SCAN DATA
			if ((data[9][i] == 0 && data[10][i] == 0))
			{
				times[0].push_back(data[0][i]);
				ras[0].push_back(15.0*data[1][i]);
				decs[0].push_back(data[2][i]);
				azimuths[0].push_back(data[3][i]);
				elevations[0].push_back(data[4][i]);
				fluxL[0].push_back(data[5][i]);
				fluxR[0].push_back(data[6][i]);
				dataDumps[0].push_back(data[7][i]);
				calibrationFlags[0].push_back(data[8][i]);
			}
			else if (data[9][i] == 0 && data[10][i] == 1 && data[8][i] == 1)
			{
				times[0].push_back(data[0][i]);
				ras[0].push_back(15.0*data[1][i]);
				decs[0].push_back(data[2][i]);
				azimuths[0].push_back(data[3][i]);
				elevations[0].push_back(data[4][i]);
				fluxL[0].push_back(data[5][i]);
				fluxR[0].push_back(data[6][i]);
				dataDumps[0].push_back(data[7][i]);
				calibrationFlags[0].push_back(data[8][i]);

				times[1].push_back(data[0][i]);
				ras[1].push_back(15.0*data[1][i]);
				decs[1].push_back(data[2][i]);
				azimuths[1].push_back(data[3][i]);
				elevations[1].push_back(data[4][i]);
				fluxL[1].push_back(data[5][i]);
				fluxR[1].push_back(data[6][i]);
				dataDumps[1].push_back(data[7][i]);
				calibrationFlags[1].push_back(data[8][i]);
			}
			else if (data[9][i] == 0 && data[10][i] == 1)
			{
				times[1].push_back(data[0][i]);
				ras[1].push_back(15.0*data[1][i]);
				decs[1].push_back(data[2][i]);
				azimuths[1].push_back(data[3][i]);
				elevations[1].push_back(data[4][i]);
				fluxL[1].push_back(data[5][i]);
				fluxR[1].push_back(data[6][i]);
				dataDumps[1].push_back(data[7][i]);
				calibrationFlags[1].push_back(data[8][i]);
			}
			else
			{
				break;
			}
			i++;
		}

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);
		scanCount = 2;

		while (i < data[0].size())
		{
			//STORE NON-TRANSITION DATA FOR REMAINING SCANS
			if (data[9][i] != scanCount - 1 && data[9][i] != 0)
			{
				times.push_back(dubFiller);
				ras.push_back(dubFiller);
				decs.push_back(dubFiller);
				azimuths.push_back(dubFiller);
				elevations.push_back(dubFiller);
				fluxL.push_back(dubFiller);
				fluxR.push_back(dubFiller);
				dataDumps.push_back(dubFiller);
				calibrationFlags.push_back(dubFiller);
				scanCount++;
			}
			else if (data[9][i] != scanCount - 1 && data[9][i] == 0)
			{
				break;
			}

			if (data[9][i] == scanCount - 1 && data[10][i] == 1)
			{
				times[scanCount].push_back(data[0][i]);
				ras[scanCount].push_back(15.0*data[1][i]);
				decs[scanCount].push_back(data[2][i]);
				azimuths[scanCount].push_back(data[3][i]);
				elevations[scanCount].push_back(data[4][i]);
				fluxL[scanCount].push_back(data[5][i]);
				fluxR[scanCount].push_back(data[6][i]);
				dataDumps[scanCount].push_back(data[7][i]);
				calibrationFlags[scanCount].push_back(data[8][i]);
			}

			i++;
		}

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);
		scanCount++;

		//POST CALIBRATION 
		while (i < data[0].size())
		{
			times[scanCount].push_back(data[0][i]);
			ras[scanCount].push_back(15.0*data[1][i]);
			decs[scanCount].push_back(data[2][i]);
			azimuths[scanCount].push_back(data[3][i]);
			elevations[scanCount].push_back(data[4][i]);
			fluxL[scanCount].push_back(data[5][i]);
			fluxR[scanCount].push_back(data[6][i]);
			dataDumps[scanCount].push_back(data[7][i]);
			calibrationFlags[scanCount].push_back(data[8][i]);

			i++;
		}
	}
	else
	{
		//HOLD SPACE FOR INITIAL CALIBRATION

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);

		//HOLD SPACE FOR ALL SCANS

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);

		//HOLD SPACE FOR END CALIBRATION

		times.push_back(dubFiller);
		ras.push_back(dubFiller);
		decs.push_back(dubFiller);
		azimuths.push_back(dubFiller);
		elevations.push_back(dubFiller);
		fluxL.push_back(dubFiller);
		fluxR.push_back(dubFiller);
		dataDumps.push_back(dubFiller);
		calibrationFlags.push_back(dubFiller);

		scanCount = 0;

		//CALIBRATION DATA DURING TRANSITION
		while (i < data[0].size())
		{
			if (data[9][i] == 0 && data[10][i] == 0 && scanCount == 0)
			{
				times[0].push_back(data[0][i]);
				ras[0].push_back(15.0*data[1][i]);
				decs[0].push_back(data[2][i]);
				azimuths[0].push_back(data[3][i]);
				elevations[0].push_back(data[4][i]);
				fluxL[0].push_back(data[5][i]);
				fluxR[0].push_back(data[6][i]);
				dataDumps[0].push_back(data[7][i]);
				calibrationFlags[0].push_back(data[8][i]);
			}
			else if (data[9][i] == 0 && data[10][i] == 1)
			{
				scanCount++;
				times[1].push_back(data[0][i]);
				ras[1].push_back(15.0*data[1][i]);
				decs[1].push_back(data[2][i]);
				azimuths[1].push_back(data[3][i]);
				elevations[1].push_back(data[4][i]);
				fluxL[1].push_back(data[5][i]);
				fluxR[1].push_back(data[6][i]);
				dataDumps[1].push_back(data[7][i]);
				calibrationFlags[1].push_back(data[8][i]);
			}
			else
			{
				times[2].push_back(data[0][i]);
				ras[2].push_back(15.0*data[1][i]);
				decs[2].push_back(data[2][i]);
				azimuths[2].push_back(data[3][i]);
				elevations[2].push_back(data[4][i]);
				fluxL[2].push_back(data[5][i]);
				fluxR[2].push_back(data[6][i]);
				dataDumps[2].push_back(data[7][i]);
				calibrationFlags[2].push_back(data[8][i]);
			}

			i++;
		}
	}
}

//calibration
void Survey::gainCalibration(Channel chan, Channel janskyChan)
{
	RCR rcr = RCR(LS_MODE_DL);
	rcr.setMuType(VALUE);
	std::vector<bool> flagHolder;
	std::vector<double> lowFluxArrayStart, highFluxArrayStart, lowFluxArrayEnd, highFluxArrayEnd;
	std::vector<double> lowDumpStart, highDumpStart, lowDumpEnd, highDumpEnd, cleanDumpStart, cleanDumpEnd;
	std::vector<double> lowTimeStart, highTimeStart, lowTimeEnd, highTimeEnd, cleanTimeStart, cleanTimeEnd;
	std::vector<double> fluxArrayStart, dumpStart, timeArrayStart, fluxArrayEnd, dumpEnd, timeArrayEnd;
	std::vector<std::vector<double> > flux;

	switch (chan)
	{
	case LEFT:
		flux = fluxL;
		break;
	case RIGHT:
		flux = fluxR;
		break;
	}

	double lowAverageStart, lowAverageEnd, highAverageStart, highAverageEnd, lowStDevStart, lowStDevEnd, highStDevStart, highStDevEnd, deltaStart, deltaEnd, averageTimeStart, averageTimeEnd;

	fluxArrayStart = flux[0];
	dumpStart = dataDumps[0];
	timeArrayStart = times[0];

	fluxArrayEnd = flux[flux.size() - 1];
	dumpEnd = dataDumps[dataDumps.size() - 1];
	timeArrayEnd = times[times.size() - 1];

	int sizeStart = fluxArrayStart.size();
	int sizeEnd = fluxArrayEnd.size();

	for (int i = 0; i < sizeStart; i++)
	{
		if (calibrationFlags[0][i] != 0 && calibrationFlags[0][i] != -1)
		{
			highFluxArrayStart.push_back(fluxArrayStart[i]);
			highDumpStart.push_back(dumpStart[i]);
			highTimeStart.push_back(timeArrayStart[i]);
		}
		else if (calibrationFlags[0][i] == -1)
		{
			continue;
		}
		else
		{
			lowFluxArrayStart.push_back(fluxArrayStart[i]);
			lowDumpStart.push_back(dumpStart[i]);
			lowTimeStart.push_back(timeArrayStart[i]);
		}
	}

	for (int i = 0; i < sizeEnd; i++)
	{
		if (calibrationFlags[calibrationFlags.size() - 1][i] != 0 && calibrationFlags[calibrationFlags.size() - 1][i] != -1)
		{
			highFluxArrayEnd.push_back(fluxArrayEnd[i]);
			highDumpEnd.push_back(dumpEnd[i]);
			highTimeEnd.push_back(timeArrayEnd[i]);
		}
		else if (calibrationFlags[calibrationFlags.size() - 1][i] == -1)
		{
			continue;
		}
		else
		{
			lowFluxArrayEnd.push_back(fluxArrayEnd[i]);
			lowDumpEnd.push_back(dumpEnd[i]);
			lowTimeEnd.push_back(timeArrayEnd[i]);
		}
	}

	if (lowFluxArrayStart.size() < 3 || highFluxArrayStart.size() < 3)
	{
		calMethod = POST;
		Debugger::print("Warn", "No reliable pre-calibration data found");
		Debugger::print("Warn", "Calibration method switched to post-calibration!");
	}
	else if (lowFluxArrayEnd.size() < 3 || highFluxArrayEnd.size() < 3)
	{
		calMethod = PRE;
		Debugger::print("Warn", "No reliable post-calibration data found");
		Debugger::print("Warn", "Calibration method switched to pre-calibration!");
	}

	if ((lowFluxArrayEnd.size() < 3 || highFluxArrayEnd.size() < 3) && (lowFluxArrayStart.size() < 3 || highFluxArrayStart.size() < 3))
	{
		Debugger::print("Error", "No reliable calibration data found!");
		//throw error for no reliable calibration data;
	}

	if (calMethod == INTERPOLATED)
	{
		double mOnStart, mOffStart, mOnEnd, mOffEnd;
		double bOnStart, bOffStart, bOnEnd, bOffEnd;
		double xBarOnStart, xBarOffStart, xBarOnEnd, xBarOffEnd;

		LinearModel model = LinearModel(highDumpStart, highTimeStart, highFluxArrayStart);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(highDumpStart, highFluxArrayStart);
		mOnStart = model.m;
		bOnStart = model.b;
		xBarOnStart = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpStart.push_back(highDumpStart[i]);
				cleanTimeStart.push_back(highTimeStart[i]);
			}
		}

		model = LinearModel(lowDumpStart, lowTimeStart, lowFluxArrayStart);

		rcr.setParametricModel(model);
		rcr.performBulkRejection(lowDumpStart, lowFluxArrayStart);
		mOffStart = model.m;
		bOffStart = model.b;
		xBarOffStart = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpStart.push_back(lowDumpStart[i]);
				cleanTimeStart.push_back(lowTimeStart[i]);
			}
		}

		averageTimeStart = Tools::getMean(cleanDumpStart, cleanTimeStart);
		deltaStart = std::abs((mOnStart*(averageTimeStart - xBarOnStart) + bOnStart) - (mOffStart*(averageTimeStart - xBarOffStart) + bOffStart));

		model = LinearModel(highDumpEnd, highTimeEnd, highFluxArrayEnd);

		rcr.setParametricModel(model);
		rcr.performBulkRejection(highDumpEnd, highFluxArrayEnd);
		mOnEnd = model.m;
		bOnEnd = model.b;
		xBarOnEnd = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpEnd.push_back(highDumpEnd[i]);
				cleanTimeEnd.push_back(highTimeEnd[i]);
			}
		}

		model = LinearModel(lowDumpEnd, lowTimeEnd, lowFluxArrayEnd);

		rcr.setParametricModel(model);
		rcr.performBulkRejection(lowDumpEnd, lowFluxArrayEnd);
		mOffEnd = model.m;
		bOffEnd = model.b;
		xBarOffEnd = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpEnd.push_back(lowDumpEnd[i]);
				cleanTimeEnd.push_back(lowTimeEnd[i]);
			}
		}

		averageTimeEnd = Tools::getMean(cleanDumpEnd, cleanTimeEnd);
		deltaEnd = std::abs((mOnEnd*(averageTimeEnd - xBarOnEnd) + bOnEnd) - (mOffEnd*(averageTimeEnd - xBarOffEnd) + bOffEnd));
	}
	else if (calMethod == PRE)
	{
		double mOnStart, mOffStart, mOnEnd, mOffEnd;
		double bOnStart, bOffStart, bOnEnd, bOffEnd;
		double xBarOnStart, xBarOffStart, xBarOnEnd, xBarOffEnd;

		LinearModel model = LinearModel(highDumpStart, highTimeStart, highFluxArrayStart);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(highDumpStart, highFluxArrayStart);
		mOnStart = model.m;
		bOnStart = model.b;
		xBarOnStart = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpStart.push_back(highDumpStart[i]);
				cleanTimeStart.push_back(highTimeStart[i]);
			}
		}

		model = LinearModel(lowDumpStart, lowTimeStart, lowFluxArrayStart);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(lowDumpStart, lowFluxArrayStart);
		mOffStart = model.m;
		bOffStart = model.b;
		xBarOffStart = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpStart.push_back(lowDumpStart[i]);
				cleanTimeStart.push_back(lowTimeStart[i]);
			}
		}

		averageTimeStart = Tools::getMean(cleanDumpStart, cleanTimeStart);
		deltaStart = std::abs((mOnStart*(averageTimeStart - xBarOnStart) + bOnStart) - (mOffStart*(averageTimeStart - xBarOffStart) + bOffStart));

	}
	else if (calMethod == POST)
	{

		double mOnStart, mOffStart, mOnEnd, mOffEnd;
		double bOnStart, bOffStart, bOnEnd, bOffEnd;
		double xBarOnStart, xBarOffStart, xBarOnEnd, xBarOffEnd;

		LinearModel model = LinearModel(highDumpEnd, highTimeEnd, highFluxArrayEnd);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(highDumpEnd, highFluxArrayEnd);
		mOnEnd = model.m;
		bOnEnd = model.b;
		xBarOnEnd = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpEnd.push_back(highDumpEnd[i]);
				cleanTimeEnd.push_back(highTimeEnd[i]);
			}
		}

		model = LinearModel(lowDumpEnd, lowTimeEnd, lowFluxArrayEnd);
		rcr.setParametricModel(model);
		rcr.performBulkRejection(lowDumpEnd, lowFluxArrayEnd);
		mOffEnd = model.m;
		bOffEnd = model.b;
		xBarOffEnd = model.xBar;

		flagHolder = rcr.result.flags;
		for (int i = 0; i < flagHolder.size(); i++)
		{
			if (flagHolder[i])
			{
				cleanDumpEnd.push_back(lowDumpEnd[i]);
				cleanTimeEnd.push_back(lowTimeEnd[i]);
			}
		}

		averageTimeEnd = Tools::getMean(cleanDumpEnd, cleanTimeEnd);
		deltaEnd = std::abs((mOnEnd*(averageTimeEnd - xBarOnEnd) + bOnEnd) - (mOffEnd*(averageTimeEnd - xBarOffEnd) + bOffEnd));
	}


	int CalMethodPostBool = 0, CalMethodPreBool = 0;
	//check these interpolations
	if (calMethod == INTERPOLATED)
	{
		CalMethodPostBool = 1;
		CalMethodPreBool = 1;
		for (int i = 1; i < times.size() - 1; i++)
		{
			for (int j = 0; j < times[i].size(); j++)
			{
				flux[i][j] = flux[i][j] / (deltaStart + (deltaEnd - deltaStart) / (averageTimeEnd - averageTimeStart)*(times[i][j] - averageTimeStart));
			}
		}
	}
	else if (calMethod == POST)
	{
		CalMethodPostBool = 1;
		for (int i = 1; i < times.size() - 1; i++)
		{
			for (int j = 0; j < times[i].size(); j++)
			{
				flux[i][j] = flux[i][j] / (deltaEnd);
			}
		}
	}
	else
	{
		CalMethodPreBool = 1;
		for (int i = 1; i < times.size() - 1; i++)
		{
			for (int j = 0; j < times[i].size(); j++)
			{
				flux[i][j] = flux[i][j] / (deltaStart);
			}
		}
	}

	switch (chan)
	{
	case LEFT:
		fluxL = flux;
		break;
	case RIGHT:
		fluxR = flux;
		break;
	}

	int janskyLeft = 0, janskyRight = 0;
	if (janskyChan == INTERPOLATED)
	{
		janskyLeft = 1;
		janskyRight = 1;
	}
	else if (janskyChan == LEFT)
	{
		janskyLeft = 1.00000;
	}
	else
	{
		janskyRight = 1.00000;
	}

	std::ofstream outputFile;
	if (CalMethodPreBool == 0)
	{
		averageTimeStart = 0;
		deltaStart = 0;
	}
	if (CalMethodPostBool == 0)
	{
		averageTimeEnd = 0;
		deltaEnd = 0;
	}

	if (chan == LEFT)
	{
		outputFile.open("Output.txt", std::ios_base::app);

		outputFile << janskyLeft << "\t";
		outputFile << CalMethodPreBool << "\t";
		outputFile << std::setprecision(2) << std::fixed << averageTimeStart << "\t";
		outputFile << deltaStart << "\n";

		outputFile << janskyLeft << "\t";
		outputFile << CalMethodPostBool << "\t";
		outputFile << std::setprecision(2) << std::fixed << averageTimeEnd << "\t";
		outputFile << deltaEnd << "\n";

		outputFile.close();
	}

	if (chan == RIGHT)
	{
		outputFile.open("Output.txt", std::ios_base::app);

		outputFile << janskyRight << "\t";
		outputFile << CalMethodPreBool << "\t";
		outputFile << std::setprecision(2) << std::fixed << averageTimeStart << "\t";
		outputFile << deltaStart << "\n";

		outputFile << janskyRight << "\t";
		outputFile << CalMethodPostBool << "\t";
		outputFile << std::setprecision(2) << std::fixed << averageTimeEnd << "\t";
		outputFile << deltaEnd << "\n";

		outputFile.close();
	}
}
void Survey::janskyCalibration(double calFactor, Channel chan)
{
	switch (chan)
	{
	case LEFT:
		for (int i = 0; i < fluxL.size(); i++)
		{
			for (int j = 0; j < fluxL[i].size(); j++)
			{
				fluxL[i][j] *= calFactor;
			}
		}
	case RIGHT:
		for (int i = 0; i < fluxL.size(); i++)
		{
			for (int j = 0; j < fluxL[i].size(); j++)
			{
				fluxR[i][j] *= calFactor;
			}
		}
	case COMPOSITE:
		for (int i = 0; i < fluxL.size(); i++)
		{
			for (int j = 0; j < fluxL[i].size(); j++)
			{
				fluxComp[i][j] *= calFactor;
			}
		}
	}
}

//data pre processing
void Survey::daisyPrelim()
{
	RCR rcr = RCR(SS_MEDIAN_DL);
	double angDistHold;
	std::vector<double> speeds, gapsAtCenter;
	std::vector<std::vector<double> > centers;

	Scan scanTemp = Scan(times[1], decs[1], ras[1], elevations[1], dataDumps[1], fluxL[1], fluxR[1], fluxComp[1]);

	partSetProcSSS.trimSize = trimSize;

	rcr.performBulkRejection(ras[1]);
	partSetProcSSS.medianRa = rcr.result.mu;
	partSetProcSSS.centerRaDeg = rcr.result.mu;

	rcr.performBulkRejection(decs[1]);
	partSetProcSSS.medianDec = rcr.result.mu;
	partSetProcSSS.centerDecDeg = rcr.result.mu;

	scanTemp.cosDecTransform(t_int, partSetProcSSS.centerDecDeg, partSetProcSSS.centerRaDeg, partSetProcSSS.centerDecDeg);
	scanTemp.updateAngDistTemp(partSetProcSSS.centerDecDeg);
	times[1] = scanTemp.getTime();
	ras[1] = scanTemp.getRa();
	decs[1] = scanTemp.getDec();

	speeds.reserve(times[1].size());
	speeds.push_back(999999);
	for (int k = 1; k < times[1].size(); k++)
	{
		angDistHold = Tools::min(Tools::getGCDistance(decs[1][k - 1], ras[1][k - 1], decs[1][k], ras[1][k], partSetProcSSS.centerDecDeg), 1) * toDeg;//GCDistance will assume a transformed coordinate system
		speeds.push_back(angDistHold / (times[1][k] - times[1][k - 1]));
	}
	this->petalPeriod = Tools::findPeriod(times[1], speeds);

	daisySweepBreaker();
}
void Survey::daisySweepBreaker()
{
	std::vector<double> dubFiller;
	std::vector<double> startTimes;
	std::vector<double> temp;
	std::vector<bool> checks;
	int i = 0, previousI;
	bool stop = false;
	double totalObservationTime = times[1][times[1].size() - 1] - times[1][0];
	numberOfPetals = std::round(totalObservationTime / petalPeriod);
	startTimes.reserve(numberOfPetals);
	startTimes.push_back(times[1][i]);
	double maxDist, distance;
	int maxDistId;
	for (int n = 1; n < numberOfPetals; n++)
	{
		maxDist = -999999;
		while (times[1][i] < times[1][0] + petalPeriod / 2.0 + petalPeriod*(n - 1))
		{
			i++;
		}
		while (i  < times[1].size() && times[1][i] < (n * petalPeriod + petalPeriod*.5 + times[1][0]))
		{
			distance = Tools::getGCDistance(decs[1][i], ras[1][i], 0, 0, partSetProcSSS.centerDecDeg);
			if (distance == distance && distance > maxDist)
			{
				maxDist = distance;
				maxDistId = i;
			}
			i++;
		}
		//DETERMINES INDEX OF FARTHEST POINT ON PEDAL
		startTimes.push_back(times[1][maxDistId]);
	}
	//PUSHBACK END OF ALL TIME VALUES INTO START TIME VECTOR
	startTimes.push_back(times[1][times[1].size() - 1]);
	i = 0;
	for (int counter = 1; counter < startTimes.size(); counter++)
	{
		previousI = i;
		//ITERATE THROUGH TIMES VECTOR UNTIL REACHES FARTHEST POINT OF PEDAL
		while (times[1][i] < startTimes[counter])
		{
			i++;
		}
		//CURRENTLY THE TIMES VECTOR HOLDS THREE VECOTRS, ONE FOR START CALIBRATION, ONE FOR END CALIBRATION, ONE FOR IN BETWEEN
		//FOR EACH PEDAL INSERT A NEW dubFiller VECTOR AT THE END OF THE TIMES VECTOR
		times.insert(times.end() - 1, dubFiller);
		//ONE BEFORE THE END OF THE TIMES VECTOR 
		//INSERT times[1][startPedalIndex] to times[1][endPedalIndex]
		times[times.size() - 2].insert(times[times.size() - 2].begin(), times[1].begin() + previousI, times[1].begin() + i);
		ras.insert(ras.end() - 1, dubFiller);
		ras[ras.size() - 2].insert(ras[ras.size() - 2].begin(), ras[1].begin() + previousI, ras[1].begin() + i);
		decs.insert(decs.end() - 1, dubFiller);
		decs[decs.size() - 2].insert(decs[decs.size() - 2].begin(), decs[1].begin() + previousI, decs[1].begin() + i);
		//John Insert
		dataDumps.insert(dataDumps.end() - 1, dubFiller);
		dataDumps[dataDumps.size() - 2].insert(dataDumps[dataDumps.size() - 2].begin(), dataDumps[1].begin() + previousI, dataDumps[1].begin() + i);
		//
		elevations.insert(elevations.end() - 1, dubFiller);
		elevations[elevations.size() - 2].insert(elevations[elevations.size() - 2].begin(), elevations[1].begin() + previousI, elevations[1].begin() + i);
		fluxL.insert(fluxL.end() - 1, dubFiller);
		fluxL[fluxL.size() - 2].insert(fluxL[fluxL.size() - 2].begin(), fluxL[1].begin() + previousI, fluxL[1].begin() + i);
		fluxR.insert(fluxR.end() - 1, dubFiller);
		fluxR[fluxR.size() - 2].insert(fluxR[fluxR.size() - 2].begin(), fluxR[1].begin() + previousI, fluxR[1].begin() + i);
		//John Insert
		fluxComp.insert(fluxComp.end() - 1, dubFiller);
		fluxComp[fluxComp.size() - 2].insert(fluxComp[fluxComp.size() - 2].begin(), fluxComp[1].begin() + previousI, fluxComp[1].begin() + i);
		//
	}

	times[times.size() - 2].insert(times[times.size() - 2].end(), times[1].begin() + i, times[1].end());
	times.erase(times.begin() + 1);
	ras[ras.size() - 2].insert(ras[ras.size() - 2].end(), ras[1].begin() + i, ras[1].end());
	ras.erase(ras.begin() + 1);
	decs[decs.size() - 2].insert(decs[decs.size() - 2].end(), decs[1].begin() + i, decs[1].end());
	decs.erase(decs.begin() + 1);
	//John Insert
	dataDumps[dataDumps.size() - 2].insert(dataDumps[dataDumps.size() - 2].end(), dataDumps[1].begin() + i, dataDumps[1].end());
	dataDumps.erase(dataDumps.begin() + 1);
	//
	elevations[elevations.size() - 2].insert(elevations[elevations.size() - 2].end(), elevations[1].begin() + i, elevations[1].end());
	elevations.erase(elevations.begin() + 1);

	fluxL[fluxL.size() - 2].insert(fluxL[fluxL.size() - 2].end(), fluxL[1].begin() + i, fluxL[1].end());
	fluxL.erase(fluxL.begin() + 1);
	fluxR[fluxR.size() - 2].insert(fluxR[fluxR.size() - 2].end(), fluxR[1].begin() + i, fluxR[1].end());
	fluxR.erase(fluxR.begin() + 1);
	//John Insert
	fluxComp[fluxComp.size() - 2].insert(fluxComp[fluxComp.size() - 2].end(), fluxComp[1].begin() + i, fluxComp[1].end());
	fluxComp.erase(fluxComp.begin() + 1);
	//



}
void Survey::findCenters()
{
	int minIndex;
	double minDist, minimum;
	std::vector<double> decHold, raHold;

	for (int i = 0; i < scans.size(); i++)
	{
		minimum = 999999;
		minIndex = -1;
		decHold = scans[i].getDec();
		raHold = scans[i].getRa();
		for (int j = 0; j < scans[i].getSize(); j++)
		{
			minDist = Tools::getGCDistance(decHold[j], raHold[j], 0, 0, partSetProcSSS.medianDec);
			if (minDist < minimum)
			{
				minimum = minDist;
				minIndex = j;
			}
		}
		scans[i].setCenter(minIndex);
	}
}

//edge calcuations
void Survey::calculateEdgeParameters()
{
	determineEdgeFlags();

	PartitionSet partSetTemp = partSetProcSSS;//survey.getPartSetProcSSS()
	//bool scansInRaTemp = survey.getScanDirection();
	//bool tracking = survey.getTracking();
	//mapType = survey.getMapType();

	double res = psfFWHM * 0.05;
	int minIndex, maxIndex;

	//EDGE ONE = RIGHT; EDGE TWO = TOP; EDGE THREE = LEFT; EDGE FOUR = BOTTOM;

	double m, b, xBar, yBar;
	std::vector<double> edgeOneRa, edgeTwoRa, edgeThreeRa, edgeFourRa;
	std::vector<double> edgeOneDec, edgeTwoDec, edgeThreeDec, edgeFourDec;
	std::vector<double> turnOneDec, turnTwoDec, turnThreeDec, turnFourDec;
	std::vector<double> turnOneRa, turnTwoRa, turnThreeRa, turnFourRa;

	partSetTemp.edgeOneParameters.clear();
	partSetTemp.edgeTwoParameters.clear();
	partSetTemp.edgeThreeParameters.clear();
	partSetTemp.edgeFourParameters.clear();

	if (mapType != DAISY)
	{
		if (mapType == NODDING) 
		{	// Skynet transforms incorrectly. As a result, we do it incorrectly.
			for (int i = 0; i < scans.size(); i++)
			{
				if (tracking)
				{
					if (pCoordinate == EQUATORIAL)
					{
						scans[i].undoDynamicCosDecTransform(0, partSetProcSSS.medianDec, partSetProcSSS.medianRa, scans[i].getDec());
					}
					else
					{
						scans[i].undoCosTransform(0, partSetProcSSS.medianDec, partSetProcSSS.medianRa);
					}
				}
				else
				{
					if (pCoordinate == EQUATORIAL)
					{
						scans[i].undoDynamicCosDecTransform(0, partSetTemp.centerDecDeg, partSetTemp.centerRaDeg, scans[i].getDec());
					}
					else
					{
						scans[i].undoCosTransform(0, partSetTemp.centerDecDeg, partSetTemp.centerRaDeg);
					}
				}
			}
		}

		if (scansInRa)
		{
			for (int i = 0; i < scans.size(); i++)
			{
				minIndex = scans[i].getMinIndex();
				maxIndex = scans[i].getMaxIndex();

				if (i == 0)
				{
					for (int j = 0; j < scans[i].getRawSize(); j++)
					{
						edgeFourRa.push_back(scans[i].getRa(j));
						edgeFourDec.push_back(scans[i].getDec(j));
					}
				}
				if (i == scans.size() - 1)
				{
					for (int j = 0; j < scans[i].getRawSize(); j++)
					{
						edgeTwoRa.push_back(scans[i].getRa(j));
						edgeTwoDec.push_back(scans[i].getDec(j));
					}
				}

				if (scans[i].getEdgePointFlag(minIndex))
				{
					edgeOneRa.push_back(scans[i].getRa(minIndex));
					edgeOneDec.push_back(scans[i].getDec(minIndex));
				}
				if (scans[i].getEdgePointFlag(maxIndex))
				{
					edgeThreeRa.push_back(scans[i].getRa(maxIndex));
					edgeThreeDec.push_back(scans[i].getDec(maxIndex));
				}

				if (scans[i].getTurningPointFlag(minIndex))//define turning points
				{
					turnOneRa.push_back(scans[i].getRa(minIndex));
					turnOneDec.push_back(scans[i].getDec(minIndex));
				}
				if (scans[i].getTurningPointFlag(maxIndex))
				{
					turnThreeRa.push_back(scans[i].getRa(maxIndex));
					turnThreeDec.push_back(scans[i].getDec(maxIndex));
				}

			}

			partSetTemp.edgeLocations.push_back(turnOneRa); // push back turning points
			partSetTemp.edgeLocations.push_back(turnOneDec);
			partSetTemp.edgeLocations.push_back(turnThreeRa);
			partSetTemp.edgeLocations.push_back(turnThreeDec);
		}
		else
		{
			for (int i = 0; i < scans.size(); i++)
			{
				minIndex = scans[i].getMinIndex();
				maxIndex = scans[i].getMaxIndex();

				if (i == 0)
				{
					for (int j = 0; j < scans[i].getRawSize(); j++)
					{
						edgeOneRa.push_back(scans[i].getRa(j));
						edgeOneDec.push_back(scans[i].getDec(j));
					}
				}
				if (i == scans.size() - 1)
				{
					for (int j = 0; j < scans[i].getRawSize(); j++)
					{
						edgeThreeRa.push_back(scans[i].getRa(j));
						edgeThreeDec.push_back(scans[i].getDec(j));
					}
				}

				if (scans[i].getEdgePointFlag(minIndex))
				{
					edgeFourRa.push_back(scans[i].getRa(minIndex));
					edgeFourDec.push_back(scans[i].getDec(minIndex));
				}

				if (scans[i].getEdgePointFlag(maxIndex))
				{
					edgeTwoRa.push_back(scans[i].getRa(maxIndex));
					edgeTwoDec.push_back(scans[i].getDec(maxIndex));
				}

				if (scans[i].getTurningPointFlag(minIndex))
				{
					turnFourRa.push_back(scans[i].getRa(minIndex));
					turnFourDec.push_back(scans[i].getDec(minIndex));
				}
				if (scans[i].getTurningPointFlag(maxIndex))
				{
					turnTwoRa.push_back(scans[i].getRa(maxIndex));
					turnTwoDec.push_back(scans[i].getDec(maxIndex));
				}
			}

			if (mapType == NODDING) // need to transform the turning points for NODDINGs only.
			{						// it's because of how cos(dec) instead of cos(dec_center) will affect it.
				for (int i = 0; i < scans.size(); i++)
				{
					if (tracking)
					{
						turnTwoDec[i] = turnTwoDec[i] - partSetProcSSS.medianDec;
						turnFourDec[i] = turnFourDec[i] - partSetProcSSS.medianDec;
						turnTwoRa[i] = (turnTwoRa[i] - partSetProcSSS.medianRa)*cos(partSetProcSSS.medianDec*toRad);
						turnFourRa[i] = (turnFourRa[i] - partSetProcSSS.medianRa)*cos(partSetProcSSS.medianDec*toRad);
					}
					else
					{
						turnTwoDec[i] = turnTwoDec[i] - partSetTemp.centerDecDeg;
						turnFourDec[i] = turnFourDec[i] - partSetTemp.centerDecDeg;
						turnTwoRa[i] = (turnTwoRa[i] - partSetTemp.centerRaDeg)*cos(partSetTemp.centerDecDeg*toRad);
						turnFourRa[i] = (turnFourRa[i] - partSetTemp.centerRaDeg)*cos(partSetTemp.centerDecDeg*toRad);
					}
				}
			}

			partSetTemp.edgeLocations.push_back(turnTwoRa);
			partSetTemp.edgeLocations.push_back(turnTwoDec);
			partSetTemp.edgeLocations.push_back(turnFourRa);
			partSetTemp.edgeLocations.push_back(turnFourDec);
		}

		//TOP AND BOTTOM EDGES FIRST

		RCR rcrTwo = RCR(LS_MODE_DL);
		LinearModel modelTwo = LinearModel(edgeTwoRa, edgeTwoDec);
		rcrTwo.setParametricModel(modelTwo);

		rcrTwo.performBulkRejection(edgeTwoDec);

		m = modelTwo.m;
		b = modelTwo.b;
		xBar = modelTwo.xBar;

		partSetTemp.edgeTwoParameters.push_back(m);
		partSetTemp.edgeTwoParameters.push_back(b);
		partSetTemp.edgeTwoParameters.push_back(xBar);

		RCR rcrFour = RCR(LS_MODE_DL);
		LinearModel modelFour = LinearModel(edgeFourRa, edgeFourDec);
		rcrFour.setParametricModel(modelFour);

		rcrFour.performBulkRejection(edgeFourDec);

		m = modelFour.m;
		b = modelFour.b;
		xBar = modelFour.xBar;

		partSetTemp.edgeFourParameters.push_back(m);
		partSetTemp.edgeFourParameters.push_back(b);
		partSetTemp.edgeFourParameters.push_back(xBar);

		//LEFT AND RIGHT EDGES; NOTE THE CHANGE IN AXES

		RCR rcrOne = RCR(LS_MODE_DL);
		LinearModel modelOne = LinearModel(edgeOneDec, edgeOneRa);
		rcrOne.setParametricModel(modelOne);

		rcrOne.performBulkRejection(edgeOneRa);

		m = modelOne.m;
		b = modelOne.b;
		yBar = modelOne.xBar;

		partSetTemp.edgeOneParameters.push_back(m);
		partSetTemp.edgeOneParameters.push_back(b);
		partSetTemp.edgeOneParameters.push_back(yBar);

		RCR rcrThree = RCR(LS_MODE_DL);
		LinearModel modelThree = LinearModel(edgeThreeDec, edgeThreeRa);
		rcrThree.setParametricModel(modelThree);

		rcrThree.performBulkRejection(edgeThreeRa);

		m = modelThree.m;
		b = modelThree.b;
		yBar = modelThree.xBar;

		partSetTemp.edgeThreeParameters.push_back(m);
		partSetTemp.edgeThreeParameters.push_back(b);
		partSetTemp.edgeThreeParameters.push_back(yBar);
	}
	else // daisies here.
	{
		partSetTemp.edgeOneParameters.resize(2);
		partSetTemp.edgeTwoParameters.resize(2);
		partSetTemp.edgeThreeParameters.resize(2);
		partSetTemp.edgeFourParameters.resize(2);

		partSetTemp.edgeOneParameters[0] = 0;
		partSetTemp.edgeOneParameters[1] = partSetTemp.minRa;

		partSetTemp.edgeTwoParameters[0] = 0;
		partSetTemp.edgeTwoParameters[1] = partSetTemp.maxDec;

		partSetTemp.edgeThreeParameters[0] = 0;
		partSetTemp.edgeThreeParameters[1] = partSetTemp.maxRa;

		partSetTemp.edgeFourParameters[0] = 0;
		partSetTemp.edgeFourParameters[1] = partSetTemp.minDec;

		double turnRa, turnDec;
		std::vector<double> turnRaVec, turnDecVec;
		for (int i = 0; i < scans.size(); i++) // define turning points for daisies.
		{
			turnRaVec.push_back(scans[i].getRa(0));
			turnDecVec.push_back(scans[i].getDec(0));
		}

		partSetTemp.edgeLocations.push_back(turnRaVec);
		partSetTemp.edgeLocations.push_back(turnDecVec);
		turnRaVec.clear();
		turnDecVec.clear();

		for (int i = 0; i < scans.size(); i++)
		{
			turnRaVec.push_back(scans[i].getRa(scans[i].getRa().size() - 1));
			turnDecVec.push_back(scans[i].getDec(scans[i].getDec().size() - 1));
		}

		partSetTemp.edgeLocations.push_back(turnRaVec);
		partSetTemp.edgeLocations.push_back(turnDecVec);
	}

	setPartSetProcSSS(partSetTemp);

}
void Survey::determineEdgeFlags()
{
	//std::vector<Scan> scans = survey.getScans();
	//PartitionSet proc = survey.getPartSetProcSSS();
	//MapTypes mapType = survey.getMapType();

	double medianDec = partSetProcSSS.medianDec;
	double medianRa = partSetProcSSS.medianRa;

	//bool scansInRaTemp = survey.getScanDirection();

	//SET ALL EDGE POINT FLAGS TO TRUE
	std::vector<bool> edgePointFlagsHold;
	for (int i = 0; i < scans.size(); i++)
	{
		edgePointFlagsHold.resize(scans[i].getRawSize(), 1);//WE USE RAW SIZE BECAUSE WE ARE USING TIME SHIFTED DATA, INCLUDING POINTS REJECTED
		scans[i].setEdgePointFlag(edgePointFlagsHold);
	}

	if (mapType != DAISY)
	{
		std::ofstream testFile;
		RCR rcr = RCR(LS_MODE_DL);
		int indexMax, indexMin;
		double angHoldMax = -999999;
		double angHoldMin = 999999;
		double angHoldAvg, raHold, decHold;
		double scanMaxesAvg, evenScanMaxesAvg, oddScanMaxesAvg;
		std::vector<double> angHold, evenScanMaxes, oddScanMaxes, evenScansMins, oddScansMins, scanMaxes, scanMins, holdTSDec;
		std::vector<int> maxIndices, minIndices;
		std::vector<bool> edgePointFlagsHold, turningPointFlagsHold, boolFiller;

		scanMins.resize(scans.size(), 0.0);
		scanMaxes.resize(scans.size(), 0.0);
		maxIndices.resize(scans.size(), 0);
		minIndices.resize(scans.size(), 0);

		for (int i = 0; i < scans.size(); i++)
		{
			angHoldMax = -999999;
			angHoldMin = 999999;
			angHold.resize(scans[i].getRawSize(), 0.0);

			if (scansInRa)
			{
				angHold = scans[i].getTSRa();
			}
			else
			{
				if (mapType == NODDING)
				{
					holdTSDec = scans[i].getTSDec();
					for (int k = 0; k < scans[i].getSize(); k++)
					{
						holdTSDec[k] = holdTSDec[k] + partSetProcSSS.centerDecDeg;
					}
					angHold = holdTSDec;
				}
				else
				{
					angHold = scans[i].getTSDec();
				}
			}

			//FOR EACH SCAN DETERMINE THE MAXIMUM AND MINIMUM RA/DEC VALUE AND THEIR INDICES

			for (int j = 0; j < scans[i].getRawSize(); j++)
			{
				if (angHoldMax <= angHold[j])
				{
					angHoldMax = angHold[j];
					indexMax = j;
				}

				if (angHoldMin >= angHold[j])
				{
					angHoldMin = angHold[j];
					indexMin = j;
				}
			}

			scanMins[i] = angHoldMin;
			scanMaxes[i] = angHoldMax;
			maxIndices[i] = indexMax;
			minIndices[i] = indexMin;

			//SEPERATE THE MAXES OF EVEN AND ODD SCANS

			if (i % 2 == 0)
			{
				evenScanMaxes.push_back(angHoldMax);
				evenScansMins.push_back(angHoldMin);
			}
			else
			{
				oddScanMaxes.push_back(angHoldMax);
				oddScansMins.push_back(angHoldMin);
			}
		}

		//CALCULATE THE AVERAGES OF ALL MAXIMUMS, THE EVEN MAXIMUMS, AND ODD MAXIMUMS

		rcr.performBulkRejection(scanMaxes);
		scanMaxesAvg = rcr.result.mu;

		rcr.performBulkRejection(evenScanMaxes);
		evenScanMaxesAvg = rcr.result.mu;

		rcr.performBulkRejection(oddScanMaxes);
		oddScanMaxesAvg = rcr.result.mu;

		//std::cout << scanMaxesAvg << "\t" << evenScanMaxesAvg << "\t" << oddScanMaxesAvg << "\n";

		//IF THE EVEN SCANS HAVE A LARGER 

		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setMinIndex(minIndices[i]);
			scans[i].setMaxIndex(maxIndices[i]);

			boolFiller.resize(scans[i].getRawSize(), 1);
			edgePointFlagsHold = boolFiller;
			turningPointFlagsHold.resize(boolFiller.size(), false);


			turningPointFlagsHold[maxIndices[i]] = 1;
			turningPointFlagsHold[minIndices[i]] = 1;

			if (evenScanMaxesAvg > oddScanMaxesAvg)
			{
				if (i % 2 == 0)
				{
					edgePointFlagsHold[maxIndices[i]] = 0; //EVEN POINTS HAVE MOST EXTREME MAX
				}
				else
				{
					edgePointFlagsHold[minIndices[i]] = 0; //ODD POINTS HAVE MOST EXTREME MIN
				}
			}
			else
			{
				if (i % 2 != 0)
				{
					edgePointFlagsHold[maxIndices[i]] = 0;
				}
				else
				{
					edgePointFlagsHold[minIndices[i]] = 0;
				}
			}

			scans[i].setTurningPointFlag(turningPointFlagsHold);
			scans[i].setEdgePointFlag(edgePointFlagsHold);

		}
	}
	else
	{
		RCR rcrRadii = RCR(LS_MODE_DL);
		int centerIndex;
		double distancePreShift, distancePostShift;
		double maxDeltaDist = -999999;
		double distance, edgeRadiusHold;
		bool timeShifted = true;
		std::vector<int> indices, intFiller;
		std::vector<double> radiiHold;
		std::vector<double> deltaDistance, dubFiller;
		std::vector<bool> edgeFlagHold, boolFiller;

		for (int i = 0; i < scans.size(); i++)
		{
			boolFiller.resize(scans[i].getRawSize(), 1);
			scans[i].setEdgePointFlag(boolFiller);
		}

		for (int i = 0; i < scans.size(); i++)
		{
			dubFiller.resize(scans[i].getRawSize(), 0.0);
			intFiller.resize(scans[i].getRawSize(), 0);
			deltaDistance = dubFiller;
			indices = intFiller;
			maxDeltaDist = -999999;

			for (int j = 0; j < scans[i].getRawSize(); j++)
			{
				distancePreShift = Tools::getGCDistance(scans[i].getRawDec(j), scans[i].getRawRa(j), medianDec, medianRa, partSetProcSSS.centerDecDeg);
				distancePostShift = Tools::getGCDistance(scans[i].getTSDec(j), scans[i].getTSRa(j), medianDec, medianRa, partSetProcSSS.centerDecDeg);
				deltaDistance[j] = distancePostShift - distancePreShift;
				if (deltaDistance[j] > maxDeltaDist)
				{
					indices[i] = j;
				}
			}
		}

		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setEdgePointFlag(indices[i], 0);
		}

		for (int i = 0; i < scans.size(); i++)
		{
			centerIndex = scans[i].getCenter();

			for (int j = 0; j < scans[i].getRawSize(); j++)
			{
				if (scans[i].getEdgePointFlag(j))
				{
					distance = Tools::getGCDistance(scans[i].getTSDec(j), scans[i].getTSRa(j), medianDec, medianRa, partSetProcSSS.centerDecDeg)*(180.0 / M_PI);
					radiiHold.push_back(distance);
				}
			}
		}

		rcrRadii.performBulkRejection(radiiHold);
		setEdgeRadius(rcrRadii.result.mu);
	}

	setScans(scans);

}

//coordinate transforms
void Survey::convertToGalacticInitial()
{
	double l, b;
	std::vector<double> bVec;
	std::vector<double> lVec;

	for (int i = 0; i < decs.size(); i++)
	{
		bVec.resize(0);
		lVec.resize(0);

		for (int j = 0; j < decs[i].size(); j++)
		{
			b = Tools::convertToB(ras[i][j], decs[i][j]);
			l = Tools::convertToL(ras[i][j], decs[i][j]);

			if (l >= 360.0) {
				l -= 360.0;
			}

			bVec.push_back(b);
			lVec.push_back(l);
		}
		decs[i] = bVec;
		ras[i] = lVec;
	}
}
void Survey::zeroCrossCheck()
{
	bool cross = false;
	int crossIndex, mp1, mp2;

	for (int i = 1; i < decs.size(); i++)
	{
		// Sweeps can be unstable near edges
		mp1 = floor(decs[i - 1].size() / 2);
		mp2 = floor(decs[i].size() / 2);

		if (scansInRa)
		{
			if (decs[i][mp2] < decs[i - 1][mp1])
			{
				cross = true;
				crossIndex = i;
				break;
			}
		}
		else
		{
			if (ras[i][mp2] < ras[i - 1][mp1])
			{
				cross = true;
				crossIndex = i;
				break;
			}
		}
	}

	// Adjust lat/long to always increase in value

	if (cross)
	{
		for (int i = 0; i < crossIndex; i++)
		{
			for (int j = 0; j < decs[i].size(); j++)
			{
				if (scansInRa)
				{
					decs[i][j] = decs[i][j] - 360.0;
				}
				else
				{
					ras[i][j] = ras[i][j] - 360.0;
				}
			}
		}
	}
}

//setters
void Survey::setRFIScale(double rfiScaleHold)
{
	this->rfiScale = rfiScaleHold;
}
void Survey::setScans(std::vector<Scan> &scansNew)
{
	this->scans = scansNew;
}
void Survey::setPartSetProcSSS(PartitionSet &partSetHold)
{
	this->partSetProcSSS = partSetHold;
}
void Survey::setPartSetProcLSS(PartitionSet &partSetHold)
{
	this->partSetProcLSS = partSetHold;
}
void Survey::setCentroidLocations(std::vector<double> &centroidLocationsHold)
{
	this->centroidLocations = centroidLocationsHold;
}
void Survey::setStandardThetaGap()
{
	if (mapType != DAISY)
	{
		RCR rcr = RCR(SS_MEDIAN_DL);
		double low, high;
		std::vector<double> deltaRa, deltaDec, deltaAng, dec, ra, lows, highs, decHold, raHold, angDistHold;
		deltaRa.reserve(scans.size() - 1);
		deltaDec.reserve(scans.size() - 1);
		deltaAng.reserve(scans.size() * scans[0].getSize());
		dec.reserve(scans.size() * scans[0].getSize());
		ra.reserve(scans.size() * scans[0].getSize());
		for (int i = 0; i < scans.size() - 1; i++)
		{
			deltaRa.push_back(std::abs(scans[i + 1].getRa(0) - scans[i].getRa(0)));
			deltaDec.push_back(std::abs(scans[i + 1].getDec(0) - scans[i].getDec(0)));
		}
		for (int i = 0; i < scans.size(); i++)
		{
			angDistHold = scans[i].getAngDist();
			for (int j = 1; j < scans[i].getSize(); j++)
			{
				deltaAng.push_back(angDistHold[j] - angDistHold[j - 1]);
			}
		}
		for (int i = 0; i < scans.size(); i++)
		{
			decHold = scans[i].getDec();
			raHold = scans[i].getRa();

			for (int j = 0; j < scans[i].getSize(); j++)
			{
				dec.push_back(decHold[j]);
				ra.push_back(raHold[j]);
			}
		}


		//WE OVERWRITE MEDIAN DEC AFTER COSINE DEC TRANSFORMATION
		RCR rcr2 = RCR(LS_MODE_DL);

		rcr.performBulkRejection(dec);
		partSetProcSSS.medianDec = rcr.result.mu;

		rcr.performBulkRejection(ra);
		partSetProcSSS.medianRa = rcr.result.mu;

		rcr2.performBulkRejection(deltaAng);
		this->medianDiffAlongSweeps = rcr2.result.mu;

		minGapThreshold = 999999;
		std::vector<bool> flagsHold = rcr2.result.flags;
		for (int i = 0; i < deltaAng.size(); i++)
		{
			if (i == 100)
			{
				Debugger::print("Info", rcr2.result.mu);
				Debugger::print("Info", rcr2.result.sigma);
			}
			if (deltaAng[i] < minGapThreshold && flagsHold[i] == true)
			{
				this->minGapThreshold = deltaAng[i];
			}
		}

		Debugger::print("Info", "MIN THETA GAP", minGapThreshold);
		Debugger::print("Info", "MU", rcr2.result.mu);
		Debugger::print("Info", "SIGMA", rcr2.result.sigma);

		rcr2.performBulkRejection(deltaRa);
		double medianRaDiff = rcr2.result.mu;

		rcr2.performBulkRejection(deltaDec);
		double medianDecDiff = rcr2.result.mu;

		if (telescope == TWENTY_METER || telescope == GBT)
		{
			this->StandardGap = Tools::max(Tools::min(medianRaDiff, medianDecDiff), medianDiffAlongSweeps);//DEGREES
			for (int i = 0; i < scans.size(); i++)
			{
				scans[i].setIntraScanGap(StandardGap);
			}
		}
		else
		{
			this->StandardGap = medianDiffAlongSweeps;//DEGREES
		}

		//this->minGapThreshold = StandardGap;

		if (scansInRa == false)
		{
			double raStart, raEnd;
			double medianRaStart, medianRaEnd;
			rcr = RCR(SS_MEDIAN_DL);
			std::vector<double> RaVecStart = scans[0].getRa();
			std::vector<double> RaVecEnd = scans[scans.size() - 1].getRa();

			rcr.performBulkRejection(RaVecStart);
			medianRaStart = rcr.result.mu;

			rcr.performBulkRejection(RaVecEnd);
			medianRaEnd = rcr.result.mu;

			for (int i = 0; i < scans.size(); i++)
			{
				scans[i].setInterScanGap(std::abs(medianRaStart - medianRaEnd) / scans.size());
			}
		}
		else
		{
			double decStart, decEnd;
			double medianDecStart, medianDecEnd;
			rcr = RCR(SS_MEDIAN_DL);
			std::vector<double> decVecStart = scans[0].getDec();
			std::vector<double> decVecEnd = scans[scans.size() - 1].getDec();

			rcr.performBulkRejection(decVecStart);
			medianDecStart = rcr.result.mu;

			rcr.performBulkRejection(decVecEnd);
			medianDecEnd = rcr.result.mu;

			for (int i = 0; i < scans.size(); i++)
			{
				scans[i].setInterScanGap(std::abs(medianDecStart - medianDecEnd) / scans.size());
			}
		}
	}
	else
	{
		RCR rcr = RCR(LS_MODE_DL);
		RCR rcrRadii = RCR(LS_MODE_DL);
		double minDist, minimum, centerIndex;
		double angDist1, angDist2;
		double distance;
		std::vector<double> gapsAtCenter, radiiHold, decHold, raHold;

		for (int i = 0; i < scans.size(); i++)
		{
			centerIndex = scans[i].getCenter(); //gapsAtCenter gets average distance from a point away from the center for each scan. DYLAN.
			gapsAtCenter.push_back(((scans[i].getAngDist(centerIndex + 1) - scans[i].getAngDist(centerIndex - 1)) / 2.0));// / psfFWHM);
		}

		rcr.performBulkRejection(gapsAtCenter);
		this->StandardGap = (2 * rcr.result.mu) / (scans.size());//DEGREES...SEE IN updateAngDist THERE IS /toRad

		minGapThreshold = 999999;
		std::vector<bool> flagsHold = rcr.result.flags;
		for (int i = 0; i < gapsAtCenter.size(); i++)
		{
			if ((gapsAtCenter[i]*4/scans.size()) < minGapThreshold && flagsHold[i] == true)
			{
				this->minGapThreshold = gapsAtCenter[i] *4 / scans.size();
			}
		}

		//this->minGapThreshold = 0.0;
		Debugger::print("Info", "MIN THETA GAP", minGapThreshold);
		Debugger::print("Info", "MU", rcr.result.mu * 4 / scans.size());
		Debugger::print("Info", "SIGMA", rcr.result.sigma * 4 / scans.size());

		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setIntraScanGap(StandardGap);
		}

		for (int i = 0; i < scans.size(); i++)
		{
			scans[i].setInterScanGap((M_PI / std::sqrt(2.0))*(partSetProcSSS.edgeRadius / scans.size()));
		}
		Debugger::print("Info", "Standard Gap in BW", StandardGap / psfFWHM);
	}


}
void Survey::setMedianDec(double value)
{
	partSetProcSSS.medianDec = value;
}
void Survey::setMedianRa(double value)
{
	partSetProcSSS.medianRa = value;
}
void Survey::setSurveyNumber(int number)
{
	this->surveyNumber = number;
	for (int i = 0; i < scans.size(); i++)
	{
		scans[i].setSurveyNumber(number);
	}
}
void Survey::setMappingCoordinate(std::string mc)
{
	this->mCoordinate = mc;
}
void Survey::set2DScatterVec(std::vector<double> scatter2dHold)
{
	this->scatter2d = scatter2dHold;
}
void Survey::setLSS2DScatterVec(std::vector<double> scatter2dHold)
{
	this->LSS2DScatter = scatter2dHold;
}
void Survey::setScansInRa(bool inRa)
{
	this->scansInRa = inRa;
}
void Survey::setEdgeRadius(double edgeRadiusVal)
{
	partSetProcSSS.edgeRadius = edgeRadiusVal;
}

void Survey::setClassificationsSSS(std::vector<std::vector<std::vector<int> > > & classVec)
{
	this->classificationsSSS = classVec;
}

void Survey::setClassificationsLSS(std::vector<std::vector<std::vector<int> > > & classVec)
{
	this->classificationsLSS = classVec;
}

void Survey::setTimeShift(double ts)
{
	this->t_int = ts;
}

//getters
MapTypes Survey::getMapType()
{
	return mapType;
}
Channel Survey::getChannel()
{
	return channel;
}
std::string Survey::getMappingCoordinate()
{
	return this->mCoordinate;
}
Coordinates Survey::getProcessingCoordinate()
{
	return pCoordinate;
}
PartitionSet Survey::getPartSetProcSSS()
{
	return partSetProcSSS;
}
PartitionSet Survey::getPartSetProcLSS()
{
	return partSetProcLSS;
}
std::vector<Scan> Survey::getScans()
{
	return scans; 
}
std::vector<double> Survey::getScatter2d()
{
	return scatter2d;
}
std::vector<double> Survey::getLSS2DScatter()
{
	return LSS2DScatter;
}
std::vector<std::vector<std::vector<int> > > Survey::getClassificationsSSS()
{
	return classificationsSSS;
}
std::vector<std::vector<std::vector<int> > > Survey::getClassificationsLSS()
{
	return classificationsLSS;
}
int Survey::getSurveyNumber()
{
	return surveyNumber;
}
double Survey::getTrimSize()
{
	return trimSize;
}
double Survey::getForcedTS()
{
	return forcedTS;
}
double Survey::getRFIScale()
{
	return rfiScale;
}
double Survey::getStandardGap()
{
	return StandardGap;
}
double Survey::getMedianRa()
{
	return partSetProcSSS.medianRa;
}
double Survey::getMedianDec()
{
	return partSetProcSSS.medianDec;
}
double Survey::getMedianLongMap()
{
	return partSetProcSSS.medianLongMap;
}
double Survey::getMedianLatiMap()
{
	return partSetProcSSS.medianLatiMap;
}
double Survey::getDiffAlongSweeps()
{
	return medianDiffAlongSweeps;
}
double Survey::getPSFFWHM()
{
	return psfFWHM;
}
double Survey::getMinGapThreshold()
{
	return minGapThreshold;
}
double Survey::getMedianDiffAlongSweeps()
{
	return medianDiffAlongSweeps;
}
double Survey::getProcEdgeRadius()
{
	return partSetProcSSS.edgeRadius;
}
double Survey::getTimeShift()
{
	return t_int;
}
double Survey::getMJD()
{
	return this->MJD;
}
bool Survey::getScanDirection()
{
	return scansInRa;
}
bool Survey::getTracking()
{
	return tracking;
}


Survey::~Survey()
{
	//dtor
}