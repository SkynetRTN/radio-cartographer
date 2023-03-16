#include "PreProcessor.h"
#include "BackgroundCUDA.h"
#include "Debugger.h"
#include "Tools.h"
#include "RCR.h"
#include <cmath>
#include <future>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

PreProcessor::PreProcessor()
{
	this->hiRes = false;
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

/* 
	Hello. I recommend rewriting PrePrcoessor. It would
	be nice if there were (at least) 3 classes. These
	are: Input, PreProcessor, and Output.

	It would greatly increase followability of the code.
*/
std::vector<std::vector<double>> PreProcessor::sdfitsReader(SpectralParameters cParams)
{
	CCfits::FITS::setVerboseMode(true);
	Spectra spectra;

	setFilename(cParams);
	determineTelescope();
	accessFITS();
	std::cout << "telescope" << telescope << std::endl;
	if (telescope == "GreenBank-20") 
	{
		if (cParams.subScale != 0.0) 
		{
			buildSpectra(spectra, cParams);
			performCleaningMulti(spectra, cParams.subScale);
			//appendColumnData();
		}

		// temp
		const double freqStep = inputData[9][0];
		spectra.setFlux(spectra20);
		spectra.determineFrequency(freqStep, frequency, centerOffset);
		// temp

		averageSpectra(spectra.getFrequencies(), cParams.inclusionBand, cParams.exclusionBand);
	}

	if (telescope == "GreenBank-20" && hiRes) 
	{
		sortHiResData(cParams.receiver);
	}
	else 
	{
		sortLowResData();
	}
	/*std::cout << "out size 2 " << outputData[0].size() << std::endl;
	for (int i = 0; i < outputData.size(); i++)
	{
		for (int j = 0; j < outputData[i].size(); j++)
		{
			if (j%2)
			{
				outputData[i].erase(outputData[i].begin() + j);
			}
		}
	}
	std::cout << "out size 2 " << outputData[0].size() << std::endl;*/

	return outputData;
}

// Input Data
int PreProcessor::accessExtTable20m(std::vector<valarray<double>> &valArraySpectra, int file)
{
	// Declarations
	std::vector<valarray<double>> vasTemp;
	std::vector<std::vector<double>> inputTemp;
	std::vector<double> MJDateHold;
	std::vector<string> columns;

	std::string hdu = "SINGLE DISH"; // name of binary table
	Debugger::print("Info", "Accessing extension table ", file);

	// Open the binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename[file - 1], CCfits::Read, hdu, false));
	CCfits::ExtHDU& table = pInfile->extension(hdu);

	columns = { "UTSECS", "CRVAL2", "CRVAL3", "AZIMUTH", "ELEVATIO", "NSAMPS", "CALSTATE", "SWPINDEX", "SWPVALID", "CDELT1" };
	inputTemp.resize(columns.size());

	table.column("DATA").readArrays(vasTemp, 1, 999999);

	for (int i = 0; i < columns.size(); i++)
	{
		table.column(columns[i]).read(inputTemp[i], 1, 999999);
	}

	if (file == 1)
	{
		valArraySpectra = vasTemp;
		inputData = inputTemp;
	}
	else if (file > 1)
	{
		for (int i = 0; i < columns.size(); i++)
		{
			inputData[i].insert(inputData[i].end(), inputTemp[i].begin(), inputTemp[i].end());
		}
		valArraySpectra.insert(valArraySpectra.end(), vasTemp.begin(), vasTemp.end());
	}
	if (file == filename.size())
	{
		table.column("MJD").read(MJDateHold, 1, 1);
		setMJD(MJDateHold);
	}

	Debugger::print("Info", inputData.size(), inputData[0].size());

	return 0;
}
int PreProcessor::accessExtTable40f(int file)
{
	// Declarations
	std::vector<std::vector<double>> inputTemp;
	std::vector<double> MJDateHold, ffdTemp;
	std::vector<string> columns;

	std::string hdu = "SINGLE DISH"; // name of binary table
	Debugger::print("Info", "Accessing extension table", file);

	// Open the binary table
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename[file - 1], CCfits::Read, hdu, false));
	CCfits::ExtHDU& table = pInfile->extension(hdu);

	// Define column names
	columns = { "UTSECS", "CRVAL2", "CRVAL3", "AZIMUTH", "ELEVATIO", "NSAMPS", "CALSTATE", "SWPINDEX", "SWPVALID", "CDELT1" };
	inputTemp.resize(columns.size());

	table.column("DATA").read(ffdTemp, 1, 999999);

	for (int i = 0; i < columns.size(); i++)
	{
		table.column(columns[i]).read(inputTemp[i], 1, 999999);
	}

	if (file == 1)
	{
		continuum = ffdTemp;
		inputData = inputTemp;
	}
	else if (file > 1)
	{
		for (int i = 0; i < columns.size(); i++)
		{
			inputData[i].insert(inputData[i].end(), inputTemp[i].begin(), inputTemp[i].end());
		}
		continuum.insert(continuum.end(), ffdTemp.begin(), ffdTemp.end());
	}
	if (file == filename.size())
	{
		table.column("MJD").read(MJDateHold, 1, 1);
		setMJD(MJDateHold);
	}

	return 0;
}
int PreProcessor::accessExtTableGBT(int file)
{
	//// Declarations
	std::vector<double> freqGBT, utc, MJDateHold, gbdTemp;
	std::vector<std::string> obsDateVec, odvTemp, obsMode, columns;
	std::vector<std::vector<double>> inputTemp;

	// Open the binary table
	std::string hdu = "SINGLE DISH";
	Debugger::print("Info", "Accessing extension table", file);

	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename[file - 1], CCfits::Read, hdu, false));
	CCfits::ExtHDU& table = pInfile->extension(hdu);

	// Column names. Some intentionally blank here.
	columns = { "", "CRVAL2", "CRVAL3", "AZIMUTH", "ELEVATIO", "NSAMPS", "", "PROCSEQN", "", ""};
	inputTemp.resize(columns.size());

	for (int i = 1; i < 8; i++)
	{
		if (i != 5 && i != 6)
		{
			table.column(columns[i]).read(inputTemp[i], 1, 99999);
		}
	}

	table.column("DATA").read(gbdTemp, 1, 9999999);
	table.column("DATE-OBS").read(odvTemp, 1, inputTemp[1].size());

	if (file == 1)
	{
		continuum = gbdTemp;
		inputData = inputTemp;
		obsDateVec = odvTemp;
	}
	else if (file > 1)
	{
		for (int i = 0; i < columns.size(); i++)
		{
			inputData[i].insert(inputData[i].end(), inputTemp[i].begin(), inputTemp[i].end());
		}
		continuum.insert(continuum.end(), gbdTemp.begin(), gbdTemp.end());
		obsDateVec.insert(obsDateVec.end(), odvTemp.begin(), odvTemp.end());
	}
	if (file == filename.size())
	{
		// Other
		table.column("OBSMODE").read(obsMode, 1, 1);
		table.column("OBSFREQ").read(freqGBT, 1, 1);
		//table.column("MJD").read(MJDateHold, 1, 1);

		// Define MJD and Frequency
		//MJDate = (int)MJDateHold[0];
		MJDate = 58708;
		frequency = freqGBT[0];

		// Format sweeps, pattern, and utc
		formatGBTInput(obsMode[0], obsDateVec, utc);

		// Define vectors for non-GBT specific values
		std::vector<double> v0, v1;
		v0.resize(inputData[1].size(), 0.0);
		v1.resize(inputData[1].size(), 1.0);

		// Define the remainder of inputData
		inputData[0] = utc;
		inputData[5] = v1;
		inputData[6] = v0;
		inputData[8] = v1;
	}

	return 0;
}
int PreProcessor::primaryHeader()
{
	std::string obsDate;
	std::string name = filename[0];
	Debugger::print("Info", "Accessing the primary header");

	// Open the primary header
	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(name, CCfits::Read, true));
	CCfits::PHDU& image = pInfile->pHDU();

	image.readKey("DATE", obsDate);
	image.readKey("OBSMODE", pattern);
	image.readKey("OBSFREQ", frequency);

	formatObservationYear(obsDate);

	if (telescope == "GreenBank-20" || telescope == "GBT")
	{
		if (year >= 2017) {
			image.readKey("COORDREF", system);
			formatCoordinateSystem();
		}
		else {
			system = "forced_equatorial";
		}

		if (telescope != "GBT")
		{
			getHistoryInfo(image.getHistory());
		}
	}

	validator();

	return 0;
}
void PreProcessor::accessFITS()
{
	const int lastFile = filename.size();
	std::vector<valarray<double>> valArraySpectra;

	for (int file = 1; file <= lastFile; file++)
	{
		if (telescope == "NRAO_GBT")
		{
			accessExtTableGBT(file);
		}
		else if (telescope == "Mighty_Forty")
		{
			accessExtTable40f(file);
		}
		else if (telescope == "GreenBank-20")
		{
			accessExtTable20m(valArraySpectra, file);
			if (file == lastFile)
			{
				primaryHeader();
				formatSpectra20(valArraySpectra);
			}
		}
	}
}
void PreProcessor::getHistoryInfo(const std::string history)
{
	// Declarations
	int numberOfRows = 0, row = 0;
	std::vector<std::string> historyVec;
	std::string line, ssLine, resLine;

	// Although history info is shown as multiple rows in header, it is single string
	for (int i = 0; i < history.size(); i++)
	{
		if (history[i] == '\n')
		{
			numberOfRows++;
		}
	}

	// Create vector instead of string
	historyVec.resize(numberOfRows);
	for (int i = 0; i < history.size(); i++)
	{
		if (history[i] == '\n')
		{
			row += 1;
			continue;
		}
		historyVec[row].push_back(history[i]);
	}

	// Find start/stop and resolution lines
	for (int i = 0; i < historyVec.size(); i++)
	{
		line = historyVec[i];
		if (line.find("START,STOP") != line.npos)
		{
			ssLine = line;
		}
		if (line.find("LOWRES") != line.npos)
		{
			resLine = line;;
		}
	}

	setResolution(resLine);
	formatHistoryInput(ssLine);
}

// Spectral Cleaning
void PreProcessor::buildSpectra(Spectra &spectra, SpectralParameters cParams)
{
	const double freqStep = inputData[9][0];

	spectra.setFlux(spectra20);
	spectra.setWeights();
	spectra.determineFrequency(freqStep, frequency, centerOffset);
	spectra.determineBaselines(cParams);
	spectra.exciseDropOuts(spectra20);

	setScatterParams(spectra);
	determineScatterMulti(spectra);
}
void PreProcessor::performCleaningMulti(Spectra spectra, double baseline)
{
	int maxThreads = 50;
	int counter = 0, completedThreads = 0, liveThreads = 0;
	std::vector<std::future<std::vector<double>>> futureVec;
	std::vector<double> results;

	BackgroundCUDA scCuda;
	futureVec.resize(spectra20.size());

	bgOutput.resize(spectra20.size());
	baseline = baseline * pow(10, 6);
	
	for (int i = 0; i < spectra20.size(); i++)
	{
		scCuda = BackgroundCUDA(spectra, i);
		Debugger::print("Info", "Cleaning Spectrum", i);

		futureVec[i] = std::async(std::launch::async, &BackgroundCUDA::calculateBG, scCuda, baseline);
		counter++;
		liveThreads++;

		if (liveThreads >= maxThreads)
		{
			for (int i = completedThreads; i < counter; i++)
			{
				results = futureVec[i].get();
				setBG(results, i);
			}
			completedThreads += liveThreads;
			liveThreads = 0;
		}
	}
	for (int i = completedThreads; i < spectra20.size(); i++)
	{
		results = futureVec[i].get();
		setBG(results, i);
	}

	if (true) // keep continuum data
	{
		spectra.keepBackground(bgOutput, inputData[9][0]);
		spectra.setSpectra20(spectra20);
	}
	else
	{
		removeBG();
	}
}
void PreProcessor::setBG(std::vector<double> bg, int index)
{
	this->bgOutput[index] = bg;
}
void PreProcessor::removeBG()
{
	std::ofstream cc1;// , cc2, cc3, cc4, cc5, cc6, cc7, cc8;
	cc1.open("12344.txt", std::ios_base::app);
	//cc2.open("12346.txt", std::ios_base::app);
	//cc3.open("68432.txt", std::ios_base::app);
	//cc4.open("74126.txt", std::ios_base::app);
	//cc5.open("74128.txt", std::ios_base::app);
	//cc6.open("74130.txt", std::ios_base::app);
	//cc7.open("112558.txt", std::ios_base::app);
	//cc8.open("112560.txt", std::ios_base::app);

	for (int i = 0; i < spectra20.size(); i++)
	{
		for (int j = 0; j < spectra20[i].size(); j++)
		{
			spectra20[i][j] = spectra20[i][j] - bgOutput[i][j];

			if (i == 1)
			{
				cc1 << bgOutput[i][j] << "\n";
			}
			/*if (i == 12346)
			{
				cc2 << bgOutput[i][j] << "\n";
			}
			if (i == 68432)
			{
				cc3 << bgOutput[i][j] << "\n";
			}
			if (i == 74126)
			{
				cc4 << bgOutput[i][j] << "\n";
			}
			if (i == 74128)
			{
				cc5 << bgOutput[i][j] << "\n";
			}
			if (i == 74130)
			{
				cc6 << bgOutput[i][j] << "\n";
			}
			if (i == 112558)
			{
				cc7 << bgOutput[i][j] << "\n";
			}
			if (i == 112560)
			{
				cc8 << bgOutput[i][j] << "\n";
			}*/
		}
	}

	cc1.close();
	//cc2.close();
	//cc3.close();
	//cc4.close();
	//cc5.close();
	//cc6.close();
	//cc7.close();
	//cc8.close();
}

// Scatter
double PreProcessor::calculateScatter(int i)
{
	RCR rcr = RCR(SS_MEDIAN_DL);
	PointToPointFunc ptpf(spectra20[i], fdists[i]);

	rcr.setNonParametricModel(ptpf);
	rcr.performRejection(weights[i], spectra20[i]);

	return 0.8197*rcr.result.sigma;
}
void PreProcessor::setScatterParams(Spectra spectra)
{
	weights = spectra.getWeights();
	fdists = spectra.getFreqDists();
}
void PreProcessor::determineScatterMulti(Spectra &spectra)
{
	int maxThreads = 125;
	int counter = 0, completedThreads = 0, liveThreads = 0;
	std::vector<std::future<double>> futureVec;
	double result;

	futureVec.resize(spectra20.size());

	for (int i = 0; i < spectra20.size(); i++)
	{
		Debugger::print("Info", i);
		futureVec[i] = std::async(std::launch::async, &PreProcessor::calculateScatter, this, i);
		counter++;
		liveThreads++;

		if (liveThreads >= maxThreads)
		{
			for (int i = completedThreads; i < counter; i++)
			{
				result = futureVec[i].get();
				spectra.setScatter(result, i);
			}
			completedThreads += liveThreads;
			liveThreads = 0;
		}
	}
	for (int i = completedThreads; i < spectra20.size(); i++)
	{
		result = futureVec[i].get();
		spectra.setScatter(result, i);
	}
}

// Formatters
void PreProcessor::formatCoordinateSystem()
{
	if (system.find("RA_DEC") != system.npos)
	{
		system = "equatorial";
	}
	else if (system.find("LNG_LAT") != system.npos)
	{
		system = "galactic";
	}
}
void PreProcessor::formatSpectra20(std::vector<valarray<double>> vaSpec)
{

	//intRanges[0] = 770;
	//intRanges[1] = 800;

	std::vector<double> vSpecHold;
	vSpecHold.resize(intRanges[1] - intRanges[0] + 1);

	/*std::ofstream cc1, cc2, cc3, cc4, cc5, cc6, cc7, cc8;
	cc1.open("12344_raw.txt", std::ios_base::app);
	cc2.open("12346_raw.txt", std::ios_base::app);
	cc3.open("68432_raw.txt", std::ios_base::app);
	cc4.open("74126_raw.txt", std::ios_base::app);
	cc5.open("74128_raw.txt", std::ios_base::app);
	cc6.open("74130_raw.txt", std::ios_base::app);
	cc7.open("112558_raw.txt", std::ios_base::app);
	cc8.open("112560_raw.txt", std::ios_base::app);*/

	for (int i = 0; i < vaSpec.size(); i++)
	{   // Convert from valarray to vector -- revese to match frequencies
		std::copy(std::begin(vaSpec[i]) + intRanges[0] - 1, std::begin(vaSpec[i]) + intRanges[1], vSpecHold.begin());
		std::reverse(vSpecHold.begin(), vSpecHold.end());

		//if (i == 12344)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc1 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 12346)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc2 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 68432)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc3 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 74126)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc4 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 74128)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc5 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 74130)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc6 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 112558)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc7 << vSpecHold[j] << "\n";
		//	}
		//}
		//if (i == 112560)
		//{
		//	for (int j = 0; j < vSpecHold.size(); j++)
		//	{
		//		cc8 << vSpecHold[j] << "\n";
		//	}
		//}

		spectra20.push_back(vSpecHold);
	}

	/*cc1.close();
	cc2.close();
	cc3.close();
	cc4.close();
	cc5.close();
	cc6.close();
	cc7.close();
	cc8.close();*/

	// Center offset from new zero-th spot
	if (intRanges[1] > 512)
	{
		centerOffset = (vaSpec[0].size() / 2) - (vaSpec[0].size() - intRanges[1]);
	}
	else
	{
		centerOffset = (vaSpec[0].size() / 2) - intRanges[0] + 1;
	}

	std::cout << "center:\t" << centerOffset << "\n";
	std::cout << "range0:\t" << intRanges[0] << "\n";
	std::cout << "range1:\t" << intRanges[1] << "\n";
	std::cout << "size:\t" << vaSpec[0].size() - intRanges[1] << "\n";
	std::cout << "size:\t" << vaSpec[0].size() / 2 << "\n";

}
void PreProcessor::formatGBTInput(std::string mode, std::vector<string> dates, std::vector<double> &utc)
{	
	int modeEnd = 0;
	double hdbl, mdbl, sdbl;
	std::string hours, minutes, seconds;
	utc.resize(dates.size());

	for (int i = 0; i < dates.size(); i++)
	{
		// Shift the sweep number
		inputData[7][i] = inputData[7][i] - 1;

		// Convert hh:mm::ss to utc
		hours = dates[i][11];
		hours += dates[i][12];

		minutes = dates[i][14];
		minutes += dates[i][15];

		for (int j = 17; j < 22; j++)
		{
			seconds += dates[i][j];
		}
		hdbl = ::atof(hours.c_str()) * 3600.00;
		mdbl = ::atof(minutes.c_str()) * 60.0;
		sdbl = ::atof(seconds.c_str());

		utc[i] = hdbl + mdbl + sdbl;

		hours.clear(), minutes.clear(), seconds.clear();
	}

	while (modeEnd < mode.size())
	{
		if (mode.at(modeEnd) == ':')
		{
			break;
		}
		pattern += mode[modeEnd];
		modeEnd++;
	}
}
void PreProcessor::formatHistoryInput(std::string hl)
{
	int startIndex, k = 0, whiteSpace = 1;
	intRanges.resize(2);

	while (k < hl.size() && whiteSpace < 3)
	{
		if (hl.at(k) == ' ')
		{
			whiteSpace++;
		}
		k++;
	}

	// Isolate the start and stop channels
	if (hl.size() != 0)
	{
		startIndex = k + 1;
		std::string temp = hl.substr(startIndex, 8);
		size_t pos = temp.find(',');

		std::string start = temp.substr(0, pos);
		std::string stop = temp.substr(pos + 2, temp.size());

		// Convert to integer
		intRanges[0] = std::stoi(start);
		intRanges[1] = std::stoi(stop); 
	}
}
void PreProcessor::formatObservationYear(std::string date)
{
	int count = 0;
	std::string yH;

	// Format the year
	while (date[count] != '-')
	{
		yH.push_back(date[count]);
		count++;
	}
	year = ::atof(yH.c_str());
}

// Misc
int PreProcessor::averageSpectra(std::vector<std::vector<double>> frequencies, std::vector<double> inclusionBand, std::vector<double> exclusionBand)
{
	int stopIndex;
	double spectraSum, inRangeSum;
	continuum.resize(spectra20.size());

	bool outOfNotches;
	double notchStart, notchStop;
	double bandStart = inclusionBand[0] * pow(10, 6);
	double bandStop  = inclusionBand[1] * pow(10, 6);

	if (exclusionBand.size() == 0) {
        exclusionBand.push_back(999999.0);
		exclusionBand.push_back(0.0);
	}

	for (int i = 0; i < spectra20.size(); i++)
	{
		spectraSum = 0.0;
		inRangeSum = 0.0;
		
		stopIndex = 0;
		notchStart = exclusionBand[0] * pow(10, 6);
		notchStop  = exclusionBand[1] * pow(10, 6);	
		outOfNotches = false;

		for (int j = 0; j < spectra20[i].size(); j++)
		{
			// Check for next notch, if not last
			if (frequencies[i][j] > notchStop && !outOfNotches)
			{
				for (int k = 0; k < exclusionBand.size(); k++)
				{
					if (frequencies[i][j] > exclusionBand[k] * pow(10, 6)) {
						stopIndex = k + 1;
					}
				}
				try {
					notchStop = exclusionBand.at(stopIndex + 1) * pow(10, 6);
					notchStart = exclusionBand.at(stopIndex) * pow(10, 6);
				}
				catch (const std::out_of_range& e) {
					outOfNotches = true;
				}
			}

			// Contribute to sum if in selected band and not in a notch
			if ((frequencies[i][j] > bandStart  && frequencies[i][j] < bandStop) &&
				(frequencies[i][j] < notchStart || frequencies[i][j] > notchStop)) 
			{
				spectraSum += spectra20[i][j];
				inRangeSum++;
			}
		}
		continuum[i] = spectraSum / inRangeSum;
	}

	return 0;
}
void PreProcessor::sortHiResData(Receiver receiver)
{
	// Input flux data is stored in a single column. We need to 
	// separate flux into left/right channels for compatabilty in Survey.

	// Declarations
	bool firstIter = true;
	int fourIndex = 0, left, right, ls, rs;
	std::vector<double> dubFiller;
	dubFiller.resize(continuum.size() / 4);
	outputData.resize(11, dubFiller);

	if (receiver == XX)
	{
		rs = 0;
		ls = 2;
	}
	else
	{
		rs = 1;
		ls = 3;
	}

	for (int i = 0; i < 11; i++)
	{
		left = ls;
		right = rs;

		if (i == 1 || i == 5 || i == 6)
		{
			continue;
		}
		for (int j = 0; j < inputData[0].size(); j += 4)
		{
			if (firstIter)
			{
				outputData[1][j / 4] = inputData[1][j] / 15.0;
				outputData[5][j / 4] = continuum[left];
				outputData[6][j / 4] = continuum[right];
				left += 4;
				right += 4;
			}
			if (i < 5)
			{
				outputData[i][j / 4] = inputData[i][j];
			}
			if (i > 6)
			{
				outputData[i][j / 4] = inputData[i - 2][j];
			}
		}
		firstIter = false;
	}
}
void PreProcessor::sortLowResData()
{
	// Input flux data is stored in a single column. We need to 
	// separate flux into left/right channels for compatabilty in Survey.

	// Declarations
	std::vector<double> dubFiller;
	dubFiller.resize((continuum.size() / 2));
	outputData.resize(11, dubFiller);

	for (int i = 0; i < 11; i++)
	{
		if (i == 1 || i == 5 || i == 6)
		{
			continue;
		}
		for (int j = 0; j < continuum.size() / 2; j++)
		{
			if (i == 0)
			{
				outputData[1][j] = inputData[1][2 * j] / 15.0;
				outputData[5][j] = continuum[2 * j + 1];
				outputData[6][j] = continuum[2 * j];
			}
			if (i < 5)
			{
				outputData[i][j] = inputData[i][2 * j];
			}
			if (i > 6)
			{
				outputData[i][j] = inputData[i - 2][2 * j];
			}
		}
	}

	/*std::ofstream cc;
	cc.open("outputData.txt", std::ios_base::app);
	for (int i = 0; i < 11; i++)
	{
		for (int j = 0; j < continuum.size() / 2; j++)
		{
			cc << outputData[i][j] << "\n";
		}
	}
	cc.close();*/

	/*int place = 0;
	int place245 = 0, place39 = 0, place357 = 0, place217 = 0;

	for (int i = 0; i < outputData[0].size(); i++)
	{
		if (outputData[9][i] == 245 - 1)
		{
			if (place245 == 144 || place245 == 146)
			{
				outputData[6][place] = 999;
				std::cout << "place245\t" << 2 * place << "\n";
			}
			place245++;
		}*/
		//if (outputData[9][i] == 39 - 1)
		//{
		//	if (place39 == 77 || place39 == 78)
		//	{
		//		outputData[6][place] = 999;
		//		std::cout << "place39\t" << 2 * place << "\n";
		//	}
		//	place39++;
		//}
		//if (outputData[9][i] == 357 - 1)
		//{
		//	if (place357 == 76 + 1 || place357 == 77 + 1)
		//	{
		//		outputData[6][place] = 999;
		//		std::cout << "place357\t" << 2 * place << "\n";
		//	}
		//	place357++;
		//}
		//if (outputData[9][i] == 217 - 1)
		//{
		//	if (place217 == 79)
		//	{
		//		outputData[6][place] = 999;
		//		std::cout << "place217\t" << 2 * place << "\n";
		//	}
		//	place217++;
		//}
		//place++;
	//}
	
	// Not sure why I ever did this?
	// Convert hours to degrees below
	//if (telescope == "MightyForty")
	//{
	//	std::transform(outputData[1].begin(), outputData[1].end(), outputData[1].begin(),
	//		std::bind1st(std::multiplies<double>(), 15.0));
	//}

}
void PreProcessor::determineTelescope()
{
	// Open the primary header
	std::cout << filename[0] << "\n";

	std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename[0], CCfits::Read, true));
	CCfits::PHDU& image = pInfile->pHDU();

	// Get the telescope name
	image.readKey("TELESCOP", telescope);
}
void PreProcessor::appendColumnData()
{
	// This function needs to seriously be thought over.
	// If the code fails for any reason when appending
	// columns, the input data file is wiped clean.
	// This function should either be removed or
	// another method needs to be found.

	//string hduName("SINGLE DISH"); // DO NOT CHANGE THIS NAME! 
	//std::valarray<double> valSpectrum, valBG;
	//valSpectrum.resize(spectra20[0].size());
	//valBG.resize(spectra20[0].size());
	//std::vector<std::valarray<double>> valSpectra20, vecValBG;

	//// Open sdfits file and open its binary table. This is risky. We are editing original data files!
	//std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename, CCfits::Write, hduName, false));
	//CCfits::ExtHDU& table = pInfile->extension(hduName);

	//// Add empty columns to sdfits file.
	//table.addColumn(CCfits::Tdouble, "BACKGROUND", spectra20[0].size());
	//table.addColumn(CCfits::Tdouble, "CLEANED", spectra20[0].size());

	//// We convert from vector<vector<>> to vector<valarray<>> here.
	//for (int i = 0; i < spectra20.size(); i++)
	//{
	//	std::copy(spectra20[i].begin(), spectra20[i].end(), std::begin(valSpectrum));
	//	valSpectra20.push_back(valSpectrum);

	//	std::copy(bgOutput[i].begin(), bgOutput[i].end(), std::begin(valBG));
	//	vecValBG.push_back(valBG);
	//}

	//// Write data to the columns.
	//table.column("CLEANED").writeArrays(valSpectra20, 1);
	//table.column("BACKGROUND").writeArrays(vecValBG, 1);

}

// Validators
void PreProcessor::validator()
{
	try {
		validateYear();
		validatePattern();
	}
	catch (const char* msg)
	{
		std::cerr << msg << std::endl;
		exit(1);
	}
	try {
		validateFrequency();
		validateSystem();
	}
	catch (const char* msg)
	{
		std::cerr << msg << std::endl;
	}
}
void PreProcessor::validateYear()
{
	if (year <= 2014.0)
	{
		throw "error: Observations taken prior to 2015 are not supported\n";
	}
}
void PreProcessor::validateSystem()
{
	if (system == "forced_equatorial")
	{
		system = "equatorial";
		throw "warning: No coordinate system info was found in the header.. The observation is assumed to be mapped in equatorial coordinates..\n";
	}
}
void PreProcessor::validatePattern()
{
	if (pattern == "track")
	{
		throw "error: Track observations are not supported\n";
	}
	if (pattern == "onoff")
	{
		throw "error: On-off observations are not supported\n";
	}
}
void PreProcessor::validateFrequency()
{
	if (frequency == 1550.0 && year < 2016)
	{
		frequency = 1395.0;
		//throw warning
	}
}

// Getters
int PreProcessor::getMJDate()
{
	return this->MJDate;
}
bool PreProcessor::getScanningDirection()
{
	return this->scanningDirection;
}
double PreProcessor::getYear()
{
	return this->year;
}
double PreProcessor::getFrequency()
{
	return this->frequency;
}
std::string PreProcessor::getTelescope()
{
	return this->telescope;
}
std::string PreProcessor::getMapPattern()
{
	return this->pattern;
}
std::string PreProcessor::getCoordinateSystem()
{
	return this->system;
}

// Setters
void PreProcessor::setMJD(std::vector<double> MJDHold)
{
	// Convert from double to int
	MJDate = (int)MJDHold[0];
}
void PreProcessor::setResolution(std::string l)
{
	if (l.find("LOWRES") != l.npos)
	{
		this->hiRes = false;
	}
	else
	{
		this->hiRes = true;
	}
}
void PreProcessor::setFilename(SpectralParameters sp)
{
	if (sp.files.size() == 0)
	{
		// throw exception
	}
	else
	{
		this->filename = sp.files;
	}
}

PreProcessor::~PreProcessor()
{
	//dtor
}


// Spectra
Spectra::Spectra()
{

}

// Baseline
void Spectra::determineBaselines(SpectralParameters sp)
{
	std::vector<double> blFiller;
	blFiller.resize(flux[0].size(), sp.subScale * pow(10, 6));

	baselines.resize(flux.size());

	for (int i = 0; i < flux.size(); i++)
	{
		baselines[i] = blFiller;
	}

	if (sp.modSubZones.size() != 0)
	{
		for (int i = 0; i < blFiller.size(); i++)
		{
			for (int k = 0; k < sp.modSubZones.size(); k += 2)
			{
				if (frequencies[0][i] >= sp.modSubZones[k] * pow(10, 6))
				{
					if (frequencies[0][i] <= sp.modSubZones[k + 1] * pow(10, 6))
					{
						blFiller[i] = sp.modSubScale * pow(10, 6);
					}
				}
			}
		}
		for (int i = 0; i < baselines.size(); i++)
		{
			baselines[i] = blFiller;
		}
	}

	velocityToFrequency(sp);

	//std::ofstream cc;
	/*cc.open("baselines.txt", std::ios_base::app);
	for (int z = 0; z < baselines[1000].size(); z++)
	{
		cc << baselines[1000][z] << "\n";
	}
	cc.close();*/
}
void Spectra::determineFrequency(double fStep, double obsFreq, double centerOffset)
{
	std::vector<double> fdTemp, fTemp;
	fdTemp.resize(flux[0].size());
	fTemp.resize(flux[0].size());

	freqDists.resize(flux.size());
	frequencies.resize(flux.size());
	fStep = std::abs(fStep);

	// Zeroth Frequency
	const double zf = (obsFreq * pow(10.0, 6.0) - (centerOffset) * fStep) + fStep;

	for (int j = 0; j < fdTemp.size(); j++)
	{
		fdTemp[j] = fStep * (double)j;
		fTemp[j] = zf + fdTemp[j];
		//std::cout << fTemp[j] << "\n";
	}

	for (int i = 0; i < freqDists.size(); i++)
	{
		freqDists[i] = fdTemp;
		frequencies[i] = fTemp;
	}

	std::cout << "Frequency size: " << frequencies[0].size() << "\n";

	// Needed if first or last point is excised later
	freqStart = frequencies[0][0];
	freqEnd = frequencies[0][frequencies[0].size() - 1];
	freqDistEnd = freqDists[0][freqDists[0].size() - 1];
}
void Spectra::velocityToFrequency(SpectralParameters sp)
{
	if (sp.velocity != 0.0)
	{
		std::vector<double> vRange(2);

		const double v = sp.velocity;
		const double c = 299792.458;

		const double f_em = 1420405751.7667;

		vRange[0] = f_em * (1 - v / c);
		vRange[1] = f_em * (1 + v / c);

		for (int i = 0; i < flux.size(); i++)
		{
			for (int j = 0; j < flux[i].size(); j++)
			{
				if (frequencies[i][j] >= vRange[0] && frequencies[i][j] <= vRange[1])
				{
					baselines[i][j] = 0.0;
				}
			}
		}
	}
}

// Background
void Spectra::interpolateDropOuts(std::vector<std::vector<double>> &bgOutput)
{
	int low, high;
	std::vector<double> bg;

	for (auto i : excisedIndices)
	{
		bg = bgOutput[i];

		for (int k = 0; k < bg.size(); k++)
		{
			if (bg[k] == 999999)
			{
				low = Tools::max(0, k - 1);
				high = Tools::min(k + 1, bg.size() - 1);
				while (bg[low] == 999999 && low > 1)
				{
					low--;
				}
				while (bg[high] == 999999 && high < bg.size() - 1)
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
		bgOutput[i] = bg;
	}
}
void Spectra::exciseDropOuts(std::vector<std::vector<double>> &spectra20)
{
	int excisedSize, index;
	bool excised;

	for (int i = 0; i < flux.size(); i++)
	{
		index = 0;
		excised = false;
		excisedSize = flux[i].size();

		while (index < excisedSize)
		{
			if (flux[i][index] < 0)
			{
				flux[i].erase(flux[i].begin() + index);
				weights[i].erase(weights[i].begin() + index);
				freqDists[i].erase(freqDists[i].begin() + index);
				frequencies[i].erase(frequencies[i].begin() + index);
				baselines[i].erase(baselines[i].begin() + index);

				spectra20[i].erase(spectra20[i].begin() + index);

				excised = true;
				excisedSize = excisedSize - 1;
			}
			else
			{
				index++;
			}
		}

		if (excised)
		{
			excisedIndices.push_back(i);
		}
	}
}
void Spectra::keepBackground(std::vector<std::vector<double>> &bgOutput, double freqStep)
{
	int numOfPoints;
	double delta;

	freqStep = std::abs(freqStep);

	std::ofstream cc1;//, cc2, cc3, cc4, cc5, cc6, cc7, cc8;
	cc1.open("12344.txt", std::ios_base::app);
	//cc2.open("12346.txt", std::ios_base::app);
	//cc3.open("68432.txt", std::ios_base::app);
	//cc4.open("74126.txt", std::ios_base::app);
	//cc5.open("74128.txt", std::ios_base::app);
	//cc6.open("74130.txt", std::ios_base::app);
	//cc7.open("112558.txt", std::ios_base::app);
	//cc8.open("112560.txt", std::ios_base::app);

	if (excisedIndices.size() != 0)
	{
		for (auto i : excisedIndices)
		{
			numOfPoints = flux[i].size();

			if (freqDists[i][0] != 0)
			{
				insertPoints(bgOutput, freqStep, i, 0, "start");
				numOfPoints += 1;
			}
			if (frequencies[i][frequencies[i].size() - 1] != freqEnd)
			{
				insertPoints(bgOutput, freqStep, i, frequencies.size() - 1, "end");
				numOfPoints += 1;
			}

			for (int j = 1; j < numOfPoints; j++)
			{
				delta = freqDists[i][j] - freqDists[i][j - 1];

				while (delta > freqStep)
				{
					insertPoints(bgOutput, freqStep, i, j, "");

					j += 1;
					numOfPoints += 1;
					delta = freqDists[i][j] - freqDists[i][j - 1];
				}
			}
		}
		interpolateDropOuts(bgOutput);
	}

	for (int i = 0; i < flux.size(); i++)
	{
		for (int j = 0; j < flux[i].size(); j++)
		{
			if (baselines[i][j] != 0.0)
			{				
				flux[i][j] = flux[i][j] - bgOutput[i][j];
			}
			if (i == 1)
			{
				cc1 << flux[i][j] << "\t";
				cc1 << bgOutput[i][j] << "\n";
			}
			if (i == 12346)
			{
				//cc2 << flux[i][j] << "\n";
				//cc2 << bgOutput[i][j] << "\n";
			}
			if (i == 68432)
			{
				//cc3 << flux[i][j] << "\n";
				//cc3 << bgOutput[i][j] << "\n";
			}
			if (i == 74126)
			{
				//cc4 << flux[i][j] << "\n";
				//cc4 << bgOutput[i][j] << "\n";
			}
			if (i == 74128)
			{
			//	cc5 << flux[i][j] << "\n";
				//cc5 << bgOutput[i][j] << "\n";
			}
			if (i == 74130)
			{
				//cc6 << flux[i][j] << "\n";
				//cc6 << bgOutput[i][j] << "\n";
			}
			if (i == 112558)
			{
				//cc7 << flux[i][j] << "\n";
				//cc7 << bgOutput[i][j] << "\n";
			}
			if (i == 112560)
			{
				//cc8 << flux[i][j] << "\n";
				//cc8 << bgOutput[i][j] << "\n";
			}
		}
	}

	cc1.close();
	//cc2.close();
	//cc3.close();
	//cc4.close();
	//cc5.close();
	//cc6.close();
	//cc7.close();
	//cc8.close();
}
void Spectra::insertPoints(std::vector<std::vector<double>> &bgOutput, double freqStep, int i, int j, std::string position)
{
	if (position == "end")
	{
		freqDists[i].push_back(freqDistEnd);
		frequencies[i].push_back(freqEnd);
		
		weights[i].push_back(1.0);
		flux[i].push_back(999999);
		baselines[i].push_back(999999);
		bgOutput[i].push_back(999999);
	}
	else
	{
		if (position == "start")
		{
			freqDists[i].insert(freqDists[i].begin(), 0);
			frequencies[i].insert(frequencies[i].begin(), freqStart);
		}
		else
		{
			freqDists[i].insert(freqDists[i].begin() + j, freqDists[i][j - 1] + freqStep);
			frequencies[i].insert(frequencies[i].begin() + j, frequencies[i][j - 1] + freqStep);
		}

		weights[i].insert(weights[i].begin() + j, 1.0);
		flux[i].insert(flux[i].begin() + j, 999999);
		baselines[i].insert(baselines[i].begin() + j, 999999);
		bgOutput[i].insert(bgOutput[i].begin() + j, 999999);
	}
}

// Getters
double Spectra::getScatter(int i)
{
	return this->scatters[i];
}
std::vector<double> Spectra::getFlux(int i)
{
	return this->flux[i];
}
std::vector<double> Spectra::getWeight(int i)
{
	return this->weights[i];
}
std::vector<double> Spectra::getFreqDist(int i)
{
	return this->freqDists[i];
}
std::vector<double> Spectra::getBaseline(int i)
{
	return this->baselines[i];
}
std::vector<std::vector<double>> Spectra::getFreqDists()
{
	return this->freqDists;
}
std::vector<std::vector<double>> Spectra::getWeights()
{
	return this->weights;
}
std::vector<std::vector<double>> Spectra::getFrequencies()
{
	return this->frequencies;
}

// Setters
void Spectra::setWeights()
{
	std::vector<double> filler;
	filler.resize(flux[0].size(), 1.0);
	weights.resize(flux.size());

	for (int i = 0; i < flux.size(); i++)
	{
		weights[i] = filler;
	}

	this->scatters.resize(flux.size());
}
void Spectra::setSpectra20(std::vector<std::vector<double>> &spectra20)
{
	spectra20.clear();
	spectra20 = flux;
}
void Spectra::setScatter(double value, int i)
{
	this->scatters[i] = value;
}
void Spectra::setFlux(std::vector<std::vector<double>> spectra20)
{
	flux = spectra20;
}

Spectra::~Spectra()
{

}