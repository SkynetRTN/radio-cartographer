#include "Map.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>
#include <cmath>

Map::Map()
{
}

// i COORDINATE CORESPONDS TO DEC

// MAP PREPROCESSING
void Map::reserveGrid(double maxDec, double minDec, double maxRa, double minRa, double res, double centerDec, double centerRa)
{
	this->maxDec = maxDec;
	this->minDec = minDec;
	this->maxRa = maxRa;
	this->minRa = minRa;
	this->centerRa = centerRa;
	this->centerDec = centerDec;
	this->resolution = res;

	int numberOfPixels = (int)(maxDec - minDec) / res;
	std::vector<Pixel> pixelFiller;
	Pixel pixelTemp;

	grid.reserve(numberOfPixels);

	for (double i = 0; i <= maxDec - minDec; i += res)
	{
		grid.push_back(pixelFiller);
		for (double j = minRa; j <= maxRa; j += res)
		{
			pixelTemp.setDec(i);
			pixelTemp.setRa(j);
			grid[round((1.0 / res) * i)].push_back(pixelTemp);
		}
	}
}

// MAP PRINTING
void Map::printSSSWeightMap()
{
	double hold;
	std::ofstream file;

	file.open("weight.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getSSSWeight();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getSSSWeight();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printSSSWeight2Map()
{
	double hold;
	std::ofstream file;

	file.open("weight.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getSSSWeight2();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getSSSWeight2();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printSSSCorrelationMap()
{
	double hold;
	std::ofstream file;

	file.open("correlation.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getSSSCorrelation();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getSSSCorrelation();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printSSSScaleMap()
{
	double hold;
	std::ofstream file;

	file.open("scale.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getSSSScale();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getSSSScale();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printSSSProcFluxMap()
{
	double hold;
	std::ofstream file;

	file.open("main.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getSSSProcFlux();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getSSSProcFlux();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
	std::cout << "[ SUCCESS ] Generating main map" << std::endl;
}

void Map::printLSSWeightMap()
{
	double hold;
	std::ofstream file;
	file.open("weight.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSWeight();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSWeight();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printLSSWeight2Map()
{
	double hold;
	std::ofstream file;
	file.open("weight.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSWeight2();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSWeight2();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printLSSCorrelationMap()
{
	double hold;
	std::ofstream file;
	file.open("correlation.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSCorrelation();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSCorrelation();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printLSSScaleMap()
{
	double hold;
	std::ofstream file;
	file.open("scale.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSScale();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSScale();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printLSSProcFluxMap()
{
	double hold;
	std::ofstream file;
	file.open("main.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSProcFlux();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSProcFlux();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}

void Map::printLSSMap()
{
	double hold;
	std::ofstream file;
	file.open("main.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getLSSFlux();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getLSSFlux();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printProcPathMap()
{
	double hold;
	std::ofstream file;
	file.open("path.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getProcPath();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getProcPath();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();

	std::cout << "[ SUCCESS ] Generating path map" << std::endl;
}
void Map::printScanNumberMap()
{
	double hold;
	std::ofstream file;
	file.open("scanNumber.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getScanNumber();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getScanNumber();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}
void Map::printIndexNumberMap()
{
	double hold;
	std::ofstream file;
	file.open("indexNumber.txt");
	for (int i = 0; i < grid.size(); i++)
	{
		hold = grid[i][(grid[i].size() - 1)].getIndexNumber();
		if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
		{
			file << hold;
		}
		else
		{
			file << "nan";
		}
		for (int j = 1; j < grid[i].size(); j++)
		{
			hold = grid[i][(grid[i].size() - 1) - j].getIndexNumber();
			if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX)
			{
				file << "\t" << hold;
			}
			else
			{
				file << "\tnan";
			}
		}
		file << "\n";
	}
	file.close();
}

void Map::printFits(const std::string fitsName, FitsFile fitsType)
{
	/*
	long naxis = 2;
	long naxes[2] = {  grid[0].size(), grid.size() };
	std::auto_ptr<FITS> pFits(0);
	long& vectorLength = naxes[0];
	long& numberOfRows = naxes[1];
	long nelements(1);
	long  fpixel(1);
	nelements = naxes[0] * naxes[1];// Find the total size of the array. calculating naxes[0]*naxes[1])

	std::valarray<double> row(vectorLength);
	std::valarray<double> array(nelements);

	try
	{
		pFits.reset(new FITS(fitsName, DOUBLE_IMG, naxis, naxes));
	}
	catch (FITS::CantCreate)
	{
		return;
	}

	for (int i = 0; i < grid.size(); i++)
	{
		for (int j = 0; j < grid[0].size(); j++)
		{
			switch (fitsType)
			{
			case (MAIN_SMALL):
				array[i*grid[0].size() + j] = getProcFlux(i, j);
				break;
			case (MAIN_LARGE) :
				array[i*grid[0].size() + j] = getProcFlux(i, j);
				break;
			case (PATH) :
				array[i*grid[0].size() + j] = getProcPath(i, j);
				break;
			case (SCALE) :
				array[i*grid[0].size() + j] = getScale(i, j);
				break;
			case (WEIGHT) :
				array[i*grid[0].size() + j] = getWeight(i, j);
				break;
			case (WEIGHT_TEMP) :
				array[i*grid[0].size() + j] = getWeight2(i, j);
				break;
			}
		}
	}

	pFits->pHDU().write(fpixel, nelements, array);

	std::cout << pFits->pHDU() << std::endl;
	return;
	*/
}

// MAP GETTERS
double Map::getMaxDec()
{
	return maxDec;
}
double Map::getMinDec()
{
	return minDec;
}
double Map::getMaxRa()
{
	return maxRa;
}
double Map::getMinRa()
{
	return minRa;
}
double Map::getCenterRa()
{
	return this->centerRa;
}
double Map::getCenterDec()
{
	return this->centerDec;
}
double Map::getResolution()
{
	return resolution;
}

// GRID GETTERS
int Map::getSize(int i)
{
	if (i == 0)
	{
		return grid.size();
	}
	else if (i == 1)
	{
		return grid[i].size();
	}
	else
	{
		std::cout << "ERROR: GRID INDEX OUT OF BOUNDS!\n";
		std::cout << "PRESS ENTER TO ABORT\n";
		// throw error here
	}
}
bool Map::getCentroidFlag(int i, int j)
{
	return grid[i][j].getCentroidFlag();
}
double Map::getDec(int i, int j)
{
	return grid[i][j].getDec();
}
double Map::getRa(int i, int j)
{
	return grid[i][j].getRa();
}
double Map::getProcPath(int i, int j)
{
	return grid[i][j].getProcPath();
}

double Map::getSSSProcFlux(int i, int j)
{
	return grid[i][j].getSSSProcFlux();
}
double Map::getSSSScale(int i, int j)
{
	return grid[i][j].getSSSScale();
}
double Map::getSSSCorrelation(int i, int j)
{
	return grid[i][j].getSSSCorrelation();
}
double Map::getSSSWeight(int i, int j)
{
	return grid[i][j].getSSSWeight();
}
double Map::getSSSWeight2(int i, int j)
{
	return grid[i][j].getSSSWeight2();
}

double Map::getLSSProcFlux(int i, int j)
{
	return grid[i][j].getLSSProcFlux();
}
double Map::getLSSScale(int i, int j)
{
	return grid[i][j].getLSSScale();
}
double Map::getLSSCorrelation(int i, int j)
{
	return grid[i][j].getLSSCorrelation();
}
double Map::getLSSWeight(int i, int j)
{
	return grid[i][j].getLSSWeight();
}
double Map::getLSSWeight2(int i, int j)
{
	return grid[i][j].getLSSWeight2();
}

int Map::getScanNumber(int i, int j)
{
	return grid[i][j].getScanNumber();
}
int Map::getIndexNumber(int i, int j)
{
	return grid[i][j].getIndexNumber();
}

// GRID SETTERS
void Map::setCentroidFlag(int i, int j, bool value)
{
	grid[i][j].setCentroidFlag(value);
}
void Map::setDec(int i, int j, double value)
{
	grid[i][j].setDec(value);
}
void Map::setRa(int i, int j, double value)
{
	grid[i][j].setRa(value);
}
void Map::setLSSFlux(int i, int j, double value)
{
	grid[i][j].setLSSFlux(value);
}
void Map::setProcPath(int i, int j, double value)
{
	grid[i][j].setProcPath(value);
}

void Map::setSSSProcFlux(int i, int j, double value)
{
	grid[i][j].setSSSProcFlux(value);
}
void Map::setSSSScale(int i, int j, double value)
{
	grid[i][j].setSSSScale(value);
}
void Map::setSSSCorrelation(int i, int j, double value)
{
	grid[i][j].setSSSCorrelation(value);
}
void Map::setSSSWeight(int i, int j, double value)
{
	grid[i][j].setSSSWeight(value);
}
void Map::setSSSWeight2(int i, int j, double value)
{
	grid[i][j].setSSSWeight2(value);
}

void Map::setLSSProcFlux(int i, int j, double value)
{
	grid[i][j].setLSSProcFlux(value);
}
void Map::setLSSScale(int i, int j, double value)
{
	grid[i][j].setLSSScale(value);
}
void Map::setLSSCorrelation(int i, int j, double value)
{
	grid[i][j].setLSSCorrelation(value);
}
void Map::setLSSWeight(int i, int j, double value)
{
	grid[i][j].setLSSWeight(value);
}
void Map::setLSSWeight2(int i, int j, double value)
{
	grid[i][j].setLSSWeight2(value);
}

void Map::setScanNumber(int i, int j, int value)
{
	return grid[i][j].setScanNumber(value);
}
void Map::setIndexNumber(int i, int j, int value)
{
	return grid[i][j].setIndexNumber(value);
}

Map::~Map()
{
}