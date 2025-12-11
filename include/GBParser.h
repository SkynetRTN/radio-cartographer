#pragma once
#include <string>
#include <vector>

class GBParser
{
public:
	GBParser();
	GBParser(std::string);
	GBParser(std::vector<std::string>);
	void setDataFile(std::string);
	std::vector<std::vector<double> > parseFile();

	//getters
	double getParamFrequency();
	double getParamMJD();
	std::string getParamCoordinate();
	std::string getParamScanType();
	std::string getParamMapType();
	
	~GBParser();
private:
	bool hasFile;
	bool hasData;

	int hashCount;//DYLAN
	int noHashLineCount;
	int paramsPerLine;
	int dataPointCount;

	//File Parameters
	std::string paramMapType;  //DYLAN
	std::string paramScanType; //DYLAN
	std::string paramCoordinate; //DYLAN
	double paramFrequency;	   //DYLAN
	double paramMJD;

	std::string dataFile;	
	std::vector<std::string> fullLineSet;
	std::vector<std::string> noHashLineSet;
	std::vector<std::string> hashLineSet;//DYLAN

	std::vector<std::vector<double> > dataPoints;
	std::vector<std::vector<double> > dataInColumns;

	void readFile();
	void removeLines();
	void findParams();//DYLAN
	void splitLines();
	void reflectMatrix();

};

