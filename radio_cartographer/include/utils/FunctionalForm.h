#pragma once
#include <vector>


class FunctionalForm
{
public:
	FunctionalForm();

	virtual void setTrueVec(std::vector<bool>&, std::vector<double>&, std::vector<double>&);
	virtual void setTrueVec(std::vector<bool>&, std::vector<double>&);
	virtual void buildModelSpace();
	virtual std::vector<double> regression();
	virtual std::vector<double> getErrors(std::vector<double>);
	virtual void setModel(std::vector<double>);

	virtual void printData();
	virtual ~FunctionalForm();

	std::vector<bool> flags;
	std::vector<int> indices;
	std::vector<double> trueW, trueY;
	std::vector<std::vector<double> > parameterSpace, weightSpace;
	std::vector<double> modelSpaceW;
};

