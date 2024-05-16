#include "io\observation\FortyFootParser.h"
#include <sstream>
#include <fstream>
#include <iostream>

FourtyParser::FourtyParser()
{
	this->hasFile = false;
	this->hasData = false;
}
FourtyParser::FourtyParser(std::string dataFile)
{
	this->dataFile = dataFile;
	this->hasFile = true;
	this->hasData = false;
	this->twoFiles = false;
}
FourtyParser::FourtyParser(std::vector<std::string> data)
{
	this->fullLineSet = data;
	this->hasData = true;
	this->hasFile = false;
}
FourtyParser::FourtyParser(std::string dataFile, std::string dataFileB)
{
	this->dataFile = dataFile;
	this->dataFileB = dataFileB;
	this->hasFile = true;
	this->hasData = false;
	this->twoFiles = true;
}

void FourtyParser::setDataFile(std::string dataFile)
{
	this->dataFile = dataFile;
	this->hasFile = true;
}
Observation FourtyParser::parseFile()
{
	try
	{
		if (hasFile)
		{
			readFile();
		}
		if (hasData || hasFile)
		{
			parseData();
			reflectMatrix40();
		}
		else
		{
			throw "No datafile or data supplied\n";
		}
	}
	catch (const char* e)
	{
		std::cerr << e;
		std::exit(EXIT_FAILURE);
	}

	Observation input(dataFile);
	input.convertFromColumnData(dataInColumns);

	// Meta data
	input.telescope = "MightyForty";

	return input;
}

void FourtyParser::reflectMatrix40()
{
	std::vector<double> dubFiller;
	dataCount = dataPoints.size();
	dubFiller.resize(dataCount);
	dataInColumns.resize(11, dubFiller);
	for (int i = 0; i < dataCount; i++)
	{
		for (int j = 0; j < 11; j++)
		{
			if (j == 0)
			{
				dataInColumns[0][i] = dataPoints[i][j];
				dataInColumns[1][i] = (dataPoints[i][j]) / 3600.0;
			}
			if (j == 1)
			{
				dataInColumns[2][i] = dataPoints[i][1];
			}
			if (j == 2)
			{
				dataInColumns[5][i] = dataPoints[i][j];
				if (twoFiles)
				{
					dataInColumns[6][i] = dataPointsB[i][j];
				}
				else
				{
					dataInColumns[6][i] = dataPoints[i][j];
				}
			}


			if (j == 3)
			{
				dataInColumns[8][i] = dataPoints[i][j];
			}
			if (j == 4)
			{
				dataInColumns[9][i] = dataPoints[i][j];
			}

		}

		dataInColumns[3][i] = 0;
		dataInColumns[4][i] = 0;
		dataInColumns[7][i] = 1;
		dataInColumns[10][i] = 0;

	}
}

void FourtyParser::readFile()
{
	std::string tempString;
	std::ifstream fileObj(dataFile.c_str());
	while (!fileObj.eof())
	{
	std::getline(fileObj, tempString);
	if (tempString.find('\r') != tempString.npos)
	{
		tempString.pop_back();
	}
	fullLineSet.push_back(tempString);
	}

	if (twoFiles)
	{
		std::string tempStringB;
		std::ifstream fileObjB(dataFileB.c_str());
		while (!fileObjB.eof())
		{
			std::getline(fileObjB, tempStringB);
			fullLineSetB.push_back(tempStringB);
		}
	}
}

void FourtyParser::removeLines()
{
	int hashCount = 0;
	while (fullLineSet[hashCount][0] == '#')
	{
		hashCount++;
	}
	noHashLineCount = fullLineSet.size() - hashCount;
	noHashLineSet.resize(noHashLineCount);
	std::copy(fullLineSet.begin() + hashCount, fullLineSet.end(), noHashLineSet.begin());
}
void FourtyParser::splitLines()
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
			while (k < line.size() && line.at(k) != '\t')
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
		
			while (k < line.size() && line.at(k) == '\t')
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

void FourtyParser::reflectMatrix()
{
	std::vector<double> dubFiller;
	dubFiller.resize(dataPointCount);
	dataInColumns.resize(paramsPerLine, dubFiller);
	for (int i = 0; i < dataPointCount; i++)
	{
		for (int j = 0; j < paramsPerLine; j++)
		{
			dataInColumns[j][i] = dataPoints[i][j];
		}
	}
}


double FourtyParser::doubleConvert(std::string input)
{
	double hold;
	std::stringstream convert(input);
	if (!(convert >> hold))
	{
		hold = NULL;
	}
	return hold;
}
void FourtyParser::parseData()
{
	bool stop = false;
	bool crosscheck = false, crosscheckB = false;
	int astCount = 0;
	int dataCount = 0;
	int doubleAstCount = 0;
	int scanCount = 0;
	int i = 0;

	std::vector<double> filler;
	filler.resize(5);

	while (!stop)
	{
		if (fullLineSet[i] == "*")
		{
			astCount++;
			i++;
			if (fullLineSet[i] == "*")
			{
				astCount = 0;
				i++;
				doubleAstCount++;
			}
		}
		if (doubleAstCount == 2)
		{
			stop = true;
			break;
		}
		dataPoints.push_back(filler);
		dataPoints[dataCount][0] = doubleConvert(fullLineSet[i]);
		dataPoints[dataCount][1] = doubleConvert(fullLineSet[i + 1]);
		dataPoints[dataCount][2] = doubleConvert(fullLineSet[i + 2]);

		if (dataCount > 0 && dataPoints[dataCount][0] < dataPoints[dataCount - 1][0])
		{
			crosscheck = true;
		}
		if (crosscheck)
		{
			dataPoints[dataCount][0] += 86400.0;
		}

		if (astCount == 0)
		{
			dataPoints[dataCount][3] = 1;
			dataPoints[dataCount][4] = -1;
		}
		if (astCount == 1)
		{
			dataPoints[dataCount][3] = 0;
			dataPoints[dataCount][4] = -1;
		}
		if (astCount >= 2)
		{
			dataPoints[dataCount][3] = -1;
			dataPoints[dataCount][4] = astCount - 2;
		}

		dataCount = dataCount + 1;
		i += 3;
	}


	if (twoFiles)
	{
		bool stopB = false;
		int astCountB = 0;
		int dataCountB = 0;
		int doubleAstCountB = 0;
		int scanCountB = 0;
		int iB = 0;
		std::vector<double> fillerB;
		fillerB.resize(5);
		while (!stopB)
		{
			if (fullLineSetB[iB] == "*")
			{
				astCountB++;
				iB++;
				if (fullLineSetB[iB] == "*")
				{
					astCountB = 0;
					iB++;
					doubleAstCountB++;
				}
			}
			if (doubleAstCountB == 2)
			{
				stopB = true;
				break;
			}
			dataPointsB.push_back(fillerB);
			dataPointsB[dataCountB][0] = doubleConvert(fullLineSetB[iB]);
			dataPointsB[dataCountB][1] = doubleConvert(fullLineSetB[iB + 1]);
			dataPointsB[dataCountB][2] = doubleConvert(fullLineSetB[iB + 2]);


			if (dataCountB > 0 && dataPointsB[dataCountB][0] < dataPointsB[dataCountB - 1][0])
			{
				crosscheckB = true;
			}
			if (crosscheckB)
			{
				dataPointsB[dataCountB][0] += 86400.0;
			}

			if (astCountB == 0)
			{
				dataPointsB[dataCountB][3] = 1;
				dataPointsB[dataCountB][4] = -1;
			}
			if (astCountB == 1)
			{
				dataPointsB[dataCountB][3] = 0;
				dataPointsB[dataCountB][4] = -1;
			}
			if (astCountB >= 2)
			{
				dataPointsB[dataCountB][3] = -1;
				dataPointsB[dataCountB][4] = astCountB - 2;
			}
			dataCountB = dataCountB + 1;
			iB += 3;
		}
	}
}


FourtyParser::~FourtyParser()
{
}
