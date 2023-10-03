#include "GBParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <float.h>

GBParser::GBParser()
{
	this->hasFile = false;
	this->hasData = false;
}

GBParser::GBParser(std::string dataFile)
{
	this->dataFile = dataFile;
	this->hasFile = true;
	this->hasData = false;
}

GBParser::GBParser(std::vector<std::string> data)
{
	this->fullLineSet = data;
	this->hasData = true;
	this->hasFile = false;
}


void GBParser::setDataFile(std::string dataFile)
{
	this->dataFile = dataFile;
	this->hasFile = true;
}
Input GBParser::parseFile()
{
	try
	{
		if (hasFile)
		{	
			readFile();
		}
		if (hasData || hasFile)
		{
			removeLines();
			splitLines();
			findParams();
			reflectMatrix();
		}
		else
		{
			throw "No datafile or data supplied\n";
		}
	}
	catch (const char* e)
	{
		std::cout << e;
	}

	Input input(dataFile);
	input.convertFromColumnData(dataInColumns);

	// Transfer the meta data
	input.telescope = "GreenBank-20";
	input.scanType = paramScanType;
	input.mapPattern = paramMapType;
	input.coordinate = paramCoordinate;
	input.observedFrequencies.emplace_back(paramFrequency);

	return input;
}


void GBParser::readFile()
{
	std::string tempString;
	std::ifstream fileObj(dataFile.c_str());
	while (!fileObj.eof())
	{
		std::getline(fileObj, tempString);
		fullLineSet.push_back(tempString);	
	}
}
void GBParser::removeLines()
{
	hashCount = 0;
	while (fullLineSet[hashCount][0] == '#')
	{
		hashCount++;
	}	
	noHashLineCount = fullLineSet.size() - hashCount;
	noHashLineSet.resize(noHashLineCount);
	std::copy(fullLineSet.begin() + hashCount, fullLineSet.end(), noHashLineSet.begin());

	//STORES THE LINES STARTING WITH A HASH MARK
	hashLineSet.resize(hashCount);
	std::copy(fullLineSet.begin(), fullLineSet.begin() + hashCount, hashLineSet.begin());
}
void GBParser::splitLines()
{
	bool stop;
	int k, startIndex, stopIndex;
	double hold;
	std::string line, tempString;
	dataPoints.resize(noHashLineCount);
	for (int i = 0; i < noHashLineCount; i++)
	{
		line = noHashLineSet[i];
		if (line.find('\n') != line.npos)
		{
			line.pop_back();
		}
		if (line.find('\r') != line.npos)
		{
			line.pop_back();
		}
		stop = false;
		k = 0;
		while (k < line.size())
		{
			startIndex = k;
			while (k < line.size() && line.at(k) != ' ')
			{
				k++;
			}
			stopIndex = k;
			std::stringstream convert(std::string(line, startIndex, stopIndex - startIndex));
			if (!(convert >> hold))
			{
				hold = NULL;
			}
			dataPoints[i].push_back(hold);
			while (k < line.size() && line.at(k) == ' ')
			{
				k++;
			}
		}
	}
		
	paramsPerLine = dataPoints[0].size();
	k = noHashLineCount - 1;
	while (dataPoints[k].size() != paramsPerLine)
	{
		dataPoints.pop_back();
		k--;
	}
	dataPointCount = k + 1;

}
void GBParser::findParams()
{
	bool stop;
	int k, startIndex, stopIndex;
	double hold;
	std::string line, stLine, cLine, fLine, tLine;
	for (int i = 0; i < hashCount; i++)
	{
		line = hashLineSet[i];
		if (line.find("SCANTYPE") != line.npos)
		{
			stLine = line;
			if (stLine.find("daisy") != stLine.npos)
			{
				paramMapType = "DAISY";
			}
			else
			{
				paramMapType = "RASTER";
				if (stLine.find("ralongmap") != line.npos)
				{
					paramScanType = "RA";
				}
				else
				{
					paramScanType = "DEC";
				}
			}
		}
		if (line.find("OBSFREQ") != line.npos)
		{
			fLine = line;
		}
		if (line.find("MJD") != line.npos)
		{
			tLine = line;
		}
		if (line.find("COORDREF") != line.npos)
		{
			if (line.find("RA_DEC_COORD") != line.npos)
			{
				paramCoordinate = "equatorial";
			}
			else if (line.find("LNG_LAT_COORD") != line.npos)
			{
				paramCoordinate = "galactic";
			}
		}
	}
	k = 0;
	while (k < fLine.size() && fLine.at(k) != '=')
	{
		k++;
	}
	startIndex = k + 1;
	std::string tempObsFreq = fLine.substr(startIndex, fLine.size());
	paramFrequency = ::atof(tempObsFreq.c_str()) / 1000.0;

	k = 0;
	while (k < tLine.size() && tLine.at(k) != '=')
	{
		k++;
	}
	startIndex = k + 1;
	std::string tempObsMJD = tLine.substr(startIndex, fLine.size());
	paramMJD = ::atof(tempObsMJD.c_str());

}

void GBParser::reflectMatrix()
{
	bool crosscheck = false; 
	double minRA = DBL_MAX;
	double maxRA = -1.0 * DBL_MAX;

	std::vector<double> dubFiller;
	dubFiller.resize(dataPointCount);
	dataInColumns.resize(paramsPerLine, dubFiller);
	for (int i = 0; i < dataPointCount; i++)
	{
		for (int j = 0; j < paramsPerLine; j++)
		{
			dataInColumns[j][i] = dataPoints[i][j];
		}
		if (dataInColumns[0][i] > 86400.0)
		{
			dataInColumns[0][i] -= 86400.0;
		}
		if (i > 0 && dataInColumns[0][i] < dataInColumns[0][i - 1])
		{
			crosscheck = true;
		}
		if (crosscheck)
		{
			dataInColumns[0][i] += 86400.0;
		}
		if (dataInColumns[1][i] < minRA)
		{
			minRA = dataInColumns[1][i];
		}
		if (dataInColumns[1][i] > maxRA)
		{
			maxRA = dataInColumns[1][i];
		}

	}
	if (minRA < 0.1 && maxRA > 23.9)
	{
		for (int i = 0; i < dataPointCount; i++)
		{
			if (dataInColumns[1][i] < 12)
			{
				dataInColumns[1][i] += 24;
			}
		}
	}

}

//getters
double GBParser::getParamFrequency()
{
	return this->paramFrequency;
}
std::string GBParser::getParamCoordinate()
{
	return this->paramCoordinate;
}
std::string GBParser::getParamScanType()
{
	return this->paramScanType;
}
std::string GBParser::getParamMapType()
{
	return this->paramMapType;
}
double GBParser::getParamMJD()
{
	return this->paramMJD;
}

GBParser::~GBParser()
{
}
