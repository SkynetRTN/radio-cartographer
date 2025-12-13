#include "Map.h"
#include <CCfits/CCfits>
#include <float.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <math.h>
#include <memory>
#include <valarray>

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

void Map::printFits(const std::string &fitsName, FitsFile fitsType,
                    const std::vector<std::pair<std::string, double>> &headerDoubles,
                    const std::vector<std::pair<std::string, std::string>> &headerStrings)
{
        if (grid.empty() || grid[0].empty())
        {
                std::cerr << "Map grid is empty; skipping FITS write for " << fitsName << std::endl;
                return;
        }

        const long width = static_cast<long>(grid.front().size());
        const long height = static_cast<long>(grid.size());
        const long naxis = 2;
        long naxes[2] = { width, height };

        std::unique_ptr<CCfits::FITS> pFits;
        try
        {
                pFits = std::make_unique<CCfits::FITS>(fitsName, CCfits::DOUBLE_IMG, naxis, naxes);
        }
        catch (const CCfits::FITS::CantCreate &)
        {
                std::cerr << "Unable to create FITS file: " << fitsName << std::endl;
                return;
        }

        std::valarray<double> array(static_cast<std::size_t>(width * height), std::numeric_limits<double>::quiet_NaN());

        auto valueForPixel = [this, fitsType](std::size_t decIdx, std::size_t raIdx) -> double
        {
                switch (fitsType)
                {
                case MAIN_SMALL:
                        return getSSSProcFlux(decIdx, raIdx);
                case MAIN_LARGE:
                        return getLSSProcFlux(decIdx, raIdx);
                case PATH:
                        return getProcPath(decIdx, raIdx);
                case SCALE:
                        return getSSSScale(decIdx, raIdx);
                case WEIGHT:
                        return getSSSWeight(decIdx, raIdx);
                case WEIGHT_TEMP:
                        return getSSSWeight2(decIdx, raIdx);
                case CORRELATION:
                        return getSSSCorrelation(decIdx, raIdx);
                default:
                        return std::numeric_limits<double>::quiet_NaN();
                }
        };

        for (std::size_t decIdx = 0; decIdx < grid.size(); ++decIdx)
        {
                const std::size_t rowWidth = grid[decIdx].size();
                for (std::size_t raIdx = 0; raIdx < rowWidth; ++raIdx)
                {
                        const double value = valueForPixel(decIdx, raIdx);
                        array[decIdx * width + raIdx] = std::isfinite(value) ? value : std::numeric_limits<double>::quiet_NaN();
                }
        }

        const long fpixel = 1;
        pFits->pHDU().write(fpixel, static_cast<long>(array.size()), array);

        auto centerPixel = [](long size) -> long
        {
                if (size % 2 == 0)
                {
                        return size / 2;
                }
                return (size - 1) / 2 + 1;
        };

        std::map<std::string, double> doubleKeys = {
                {"CRPIX1", static_cast<double>(centerPixel(width))},
                {"CRPIX2", static_cast<double>(centerPixel(height))},
                {"CRVAL1", centerRa},
                {"CRVAL2", centerDec},
                {"CDELT1", resolution},
                {"CDELT2", resolution},
        };

        for (const auto &entry : headerDoubles)
        {
                doubleKeys[entry.first] = entry.second;
        }

        std::map<std::string, std::string> stringKeys = {
                {"CTYPE1", "RA---TAN"},
                {"CTYPE2", "DEC--TAN"},
                {"CUNIT1", "deg"},
                {"CUNIT2", "deg"},
        };

        for (const auto &entry : headerStrings)
        {
                stringKeys[entry.first] = entry.second;
        }

        CCfits::PHDU &hdu = pFits->pHDU();
        for (const auto &entry : doubleKeys)
        {
                hdu.addKey(entry.first, entry.second, "");
        }
        for (const auto &entry : stringKeys)
        {
                hdu.addKey(entry.first, entry.second, "");
        }
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
                if (grid.empty())
                {
                        return 0;
                }
                return grid[0].size();
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