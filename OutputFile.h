#pragma once
#include "Survey.h"
#include <fstream>

class Output
{
public:
	Output();

	std::ofstream outputFile;

	void printTelescopeInfo(SurveyParameters);
	void printProcessingHeader();
	void printCalibrationHeader();
	void printWScale(double);

	void printPhotometryHolder();
	void printPhotometry(std::vector<double>);
	void printBgSubtractionInfo(double);
	void printTimeShiftInfo(bool, double);
	void printRfiInfo(double);

	std::string wScale;
	std::vector<std::string> parameters;

	~Output();
};
