#include "FunctionalForm.h"
#include <iostream>
FunctionalForm::FunctionalForm()
{
}

void FunctionalForm::setTrueVec(std::vector<bool> &flags, std::vector<double> &w, std::vector<double> &y)
{

	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueWVec, trueYVec;
	this->flags = flags;

	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueWVec.resize(trueCount);
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueWVec[currentIndex] = (w[i]);
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;

			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	trueW = trueWVec;
	indices = indicesVec;
}
void FunctionalForm::setTrueVec(std::vector<bool> &flags, std::vector<double> &y)
{
	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueYVec;
	this->flags = flags;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;
			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	indices = indicesVec;

}
void FunctionalForm::buildModelSpace()
{
	parameterSpace.push_back(trueY);
}
std::vector<double> FunctionalForm::regression()
{
	std::vector<double> ret;
	return ret;
}
std::vector<double> FunctionalForm::getErrors(std::vector<double> x)
{
	std::vector<double> ret;
	return ret;
}
void FunctionalForm::printData()
{

}
void FunctionalForm::setModel(std::vector<double> x)
{

}
FunctionalForm::~FunctionalForm()
{
}
