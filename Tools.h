#pragma once
#include <vector>
#include <iostream>
#include <fstream>



namespace Tools
{
	
	void sort(std::vector<double>&);
	void sortAll(int, int, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	void sortAll(int, int, std::vector<double>&, std::vector<int>&, std::vector<int>&);
	void sortAll(int, int, std::vector<int>&, std::vector<int>&, std::vector<bool>&);
	void bubbleSort2(std::vector<double> &, std::vector<int> &, std::vector<double> &, std::vector<int> &, std::vector<int> &, std::vector<double> &);
	void bubbleSort3(std::vector<std::vector<double> > &);

	//coordinate conversions
	double convertToDec(double, double); //DYLAN
	double convertToRa(double, double);  //DYLAN
	double convertToB(double, double);   //DYLAN
	double convertToL(double, double);   //DYLAN

	//statistics tools
	double max(double, double);
	double min(double, double);
	double get68th(std::vector<double>&);
	double getMean(std::vector<double>);
	double getMean(std::vector<double>, std::vector<double>);
	double getMedian(std::vector<double>&);
	double getStDev(double, double, std::vector<double>);
	double getWStDev(double, double, std::vector<double>, std::vector<double>);
	double erfc(double);

	//math tools
	double signedArcTan(double, double);
	double arcCos(double);
	std::vector<std::vector<double> > pivotSystem(int, std::vector<double>&, std::vector<double>&);
	std::vector<double> matrixSolver(int, std::vector<double>, std::vector<double>);

	//regression without pivoting
	//std::vector<double> fixedRegression(int, int, std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double, double);
	//std::vector<double> fixedQuadraticRegression(int, int, std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double, double);
	
	//regression with pivoting
	std::vector<double> regressionPivot(std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::vector<double> quadraticRegressionPivot(std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	std::vector<double> fixedRegressionPivot(std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double, double);
	std::vector<double> fixedQuadraticRegressionPivot(std::vector<bool>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, double, double);

	//misc
	double findPeriod(std::vector<double>&, std::vector<double>&);
	std::vector<double> daisyAngleBuilder(int, int, double, std::vector<double>, std::vector<double>);
	std::vector<std::vector<double> > crossCorrelate(std::vector<double>, std::vector<double>);
	void lasso(int, std::vector<int> &, std::vector<double> &, std::vector<double> &, std::vector<double> &, std::vector<double> &);
	double getGCDistance(double, double, double, double,double);
	double getModGCDistance(double, double, double, double);
	double getPythDistance(double, double, double, double);
	int determinePixel(double, double, double);
	int determineNearestIndex(int, std::vector<double> &, std::vector<double> &);
	void interpolate(std::vector<double>&, double);
}

