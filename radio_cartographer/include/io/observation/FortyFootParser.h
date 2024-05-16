#pragma once
#include "io\Observation.h"
#include <string>
#include <vector>

class FourtyParser
{
public:
	FourtyParser();
	FourtyParser(std::string);
	FourtyParser(std::vector<std::string>);
	FourtyParser(std::string, std::string);
	void setDataFile(std::string);
	Observation parseFile();
	~FourtyParser();

private:
	bool hasFile;
	bool hasData;
	bool twoFiles;
	int noHashLineCount;
	int paramsPerLine;
	int dataPointCount;
	int telescope;
	int dataCount;

	std::string dataFile;
	std::string dataFileB;
	std::vector<std::string> fullLineSet;
	std::vector<std::string> fullLineSetB;
	std::vector<std::string> noHashLineSet;

	std::vector<std::vector<double> > dataPoints;
	std::vector<std::vector<double>> dataPointsB;
	std::vector<std::vector<double> > dataInColumns;


	void readFile();
	void removeLines();
	void splitLines();
	void reflectMatrix();
	void reflectMatrix40();
	void parseData();
	double doubleConvert(std::string);
};


