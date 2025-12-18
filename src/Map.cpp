#include "Map.h"
#include <cmath>
#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>

Map::Map() {}

// i COORDINATE CORESPONDS TO DEC

// MAP PREPROCESSING
void Map::reserveGrid(double maxDec, double minDec, double maxRa, double minRa,
                      double res, double centerDec, double centerRa) {
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

  for (double i = 0; i <= maxDec - minDec; i += res) {
    grid.push_back(pixelFiller);
    for (double j = minRa; j <= maxRa; j += res) {
      pixelTemp.setDec(i);
      pixelTemp.setRa(j);
      grid[round((1.0 / res) * i)].push_back(pixelTemp);
    }
  }
}

// MAP PRINTING
void Map::printSSSWeightMap() {
  double hold;
  std::ofstream file;

  file.open("weight.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getSSSWeight();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getSSSWeight();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printSSSWeight2Map() {
  double hold;
  std::ofstream file;

  file.open("weight.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getSSSWeight2();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getSSSWeight2();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printSSSCorrelationMap() {
  double hold;
  std::ofstream file;

  file.open("correlation.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getSSSCorrelation();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getSSSCorrelation();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printSSSScaleMap() {
  double hold;
  std::ofstream file;

  file.open("scale.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getSSSScale();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getSSSScale();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printSSSProcFluxMap() {
  double hold;
  std::ofstream file;

  file.open("main.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getSSSProcFlux();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getSSSProcFlux();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}

void Map::printLSSWeightMap() {
  double hold;
  std::ofstream file;
  file.open("weight.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSWeight();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSWeight();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printLSSWeight2Map() {
  double hold;
  std::ofstream file;
  file.open("weight.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSWeight2();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSWeight2();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printLSSCorrelationMap() {
  double hold;
  std::ofstream file;
  file.open("correlation.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSCorrelation();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSCorrelation();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printLSSScaleMap() {
  double hold;
  std::ofstream file;
  file.open("scale.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSScale();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSScale();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printLSSProcFluxMap() {
  double hold;
  std::ofstream file;
  file.open("main.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSProcFlux();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSProcFlux();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}

void Map::printLSSMap() {
  double hold;
  std::ofstream file;
  file.open("main.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getLSSFlux();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getLSSFlux();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printProcPathMap() {
  double hold;
  std::ofstream file;
  file.open("path.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getProcPath();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getProcPath();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printScanNumberMap() {
  double hold;
  std::ofstream file;
  file.open("scanNumber.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getScanNumber();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getScanNumber();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}
void Map::printIndexNumberMap() {
  double hold;
  std::ofstream file;
  file.open("indexNumber.txt");
  for (int i = 0; i < grid.size(); i++) {
    hold = grid[i][(grid[i].size() - 1)].getIndexNumber();
    if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
      file << hold;
    } else {
      file << "nan";
    }
    for (int j = 1; j < grid[i].size(); j++) {
      hold = grid[i][(grid[i].size() - 1) - j].getIndexNumber();
      if (hold == hold && hold <= DBL_MAX && hold >= -DBL_MAX) {
        file << "\t" << hold;
      } else {
        file << "\tnan";
      }
    }
    file << "\n";
  }
  file.close();
}

// Helper to get data for specific layers
void getLayerData(Map *map, FitsFile type, long nelements, long naxes0,
                  std::valarray<double> &array) {
  for (int i = 0; i < map->getSize(0); i++) {
    for (int j = 0; j < map->getSize(1); j++) {
      double val = NAN;
      switch (type) {
      case MAIN_SMALL:
        val = map->getSSSProcFlux(i, j);
        break;
      case PATH:
        val = map->getProcPath(i, j);
        break;
      case SCALE:
        val = map->getSSSScale(i, j);
        break;
      case WEIGHT:
        val = map->getSSSWeight(i, j);
        break;
      case WEIGHT_TEMP:
        val = map->getSSSCorrelation(i, j);
        break; // Mapping Correlation to this for loop
      default:
        break;
      }
      // Reverse j-index to match FITS convention (CDELT1 is negative)
      int reversed_j = (map->getSize(1) - 1) - j;
      array[i * naxes0 + reversed_j] = val;
    }
  }
}

void Map::printFits(const std::string fitsName,
                    std::map<std::string, std::string> headerInfo,
                    Map *rawMap) {
  using namespace CCfits;

  if (grid.empty() || grid[0].empty()) {
    std::cerr << "Map grid is empty, cannot write FITS." << std::endl;
    return;
  }

  long naxis = 2;
  long naxes[2] = {(long)grid[0].size(),
                   (long)grid.size()}; // NAXIS1=RA, NAXIS2=DEC
  long nelements = naxes[0] * naxes[1];
  std::valarray<double> array(nelements);

  try {
    // ! overwrites existing file
    std::auto_ptr<FITS> pFits(
        new FITS("!" + fitsName, DOUBLE_IMG, naxis, naxes));

    // Write Primary HDU (Main Flux)
    getLayerData(this, MAIN_SMALL, nelements, naxes[0], array);
    pFits->pHDU().write(1, nelements, array);

    // Helper lambda or macro not ideal in C++98, so just function call for WCS
    // setting But we need to set keys on pHDU or ExtHDU. Both inherit from HDU.
    // Let's refactor WCS setting to a helper function if possible, or duplicate
    // logic for clarity.

    HDU *currentHDU = &pFits->pHDU();

    // WCS for Primary
    double centerPix1 = (naxes[0] - 1) / 2.0 + 1.0;
    double centerPix2 = (naxes[1] - 1) / 2.0 + 1.0;

    currentHDU->addKey("CTYPE1", "RA---TAN", "Right Ascension");
    currentHDU->addKey("CRVAL1", centerRa, "Reference value (Center RA)");
    currentHDU->addKey("CRPIX1", centerPix1, "Reference pixel (Center X)");
    currentHDU->addKey("CDELT1", -resolution, "Pixel scale (deg)");
    currentHDU->addKey("CUNIT1", "deg", "");
    currentHDU->addKey("CTYPE2", "DEC--TAN", "Declination");
    currentHDU->addKey("CRVAL2", centerDec, "Reference value (Center Dec)");
    currentHDU->addKey("CRPIX2", centerPix2, "Reference pixel (Center Y)");
    currentHDU->addKey("CDELT2", resolution, "Pixel scale (deg)");
    currentHDU->addKey("CUNIT2", "deg", "");
    currentHDU->addKey("CENTERX", centerPix1, "");
    currentHDU->addKey("CENTERY", centerPix2, "");
    currentHDU->addKey("CENTERRA", centerRa, "");
    currentHDU->addKey("CENTERDE", centerDec, "");
    currentHDU->addKey("PIXLSIZE", resolution, "");
    currentHDU->addKey("BEAM", 1.0, "Beam size in pixels (approx)");

    // Metadata for Primary
    // color map key for visualization purposes for Afterglow Access:
    // https://afterglow.skynet.unc.edu/
    currentHDU->addKey("AGCMAP", "Rainbow", "Map color for Afterglow");
    currentHDU->addKey("EXTNAME", "main", "Extension name");

    for (std::map<std::string, std::string>::const_iterator it =
             headerInfo.begin();
         it != headerInfo.end(); ++it) {
      currentHDU->addKey(it->first, it->second, "");
    }

    // Layers to write as Extensions
    std::vector<std::pair<std::string, FitsFile>> layers;
    layers.push_back(std::make_pair("path", PATH));
    layers.push_back(std::make_pair("scale", SCALE));
    layers.push_back(std::make_pair("weight", WEIGHT));
    layers.push_back(std::make_pair("correlation", WEIGHT_TEMP));

    for (size_t l = 0; l < layers.size(); ++l) {
      std::string extName = layers[l].first;
      FitsFile type = layers[l].second;

      getLayerData(this, type, nelements, naxes[0], array);

      std::vector<long> extNaxes;
      extNaxes.push_back(naxes[0]);
      extNaxes.push_back(naxes[1]);

      ExtHDU *imageExt = pFits->addImage(extName, DOUBLE_IMG, extNaxes);
      imageExt->write(1, nelements, array);
      // overwrite phdu for Afterglow color map
      imageExt->addKey("AGCMAP", "Gray", "Map color for Afterglow");

      // Metadata for Extension
      for (std::map<std::string, std::string>::const_iterator it =
               headerInfo.begin();
           it != headerInfo.end(); ++it) {
        imageExt->addKey(it->first, it->second, "");
      }
    }

    if (rawMap != nullptr) {
      long rawWidth = rawMap->getSize(1);
      long rawHeight = rawMap->getSize(0);
      long rawNelements = rawWidth * rawHeight;
      std::valarray<double> rawArray(rawNelements);

      getLayerData(rawMap, MAIN_SMALL, rawNelements, rawWidth, rawArray);

      std::vector<long> rawExtNaxes;
      rawExtNaxes.push_back(rawWidth);
      rawExtNaxes.push_back(rawHeight);

      ExtHDU *imageExt = pFits->addImage("raw", DOUBLE_IMG, rawExtNaxes);
      imageExt->write(1, rawNelements, rawArray);
      imageExt->addKey("AGCMAP", "Rainbow", "Map color for Afterglow");

      for (std::map<std::string, std::string>::const_iterator it =
               headerInfo.begin();
           it != headerInfo.end(); ++it) {
        if (it->first == "RCRAW") {
          imageExt->addKey("RCRAW", 0, "RC:  Is this the raw file");
        } else {
          imageExt->addKey(it->first, it->second, "");
        }
      }
    }
  } catch (FITS::CantCreate &e) {
    std::cerr << "CCfits Error: Cannot create file " << fitsName << " : "
              << e.message() << std::endl;
    return;
  } catch (FitsException &e) {
    std::cerr << "CCfits Exception: " << e.message() << std::endl;
    return;
  }
}

// MAP GETTERS
double Map::getMaxDec() { return maxDec; }
double Map::getMinDec() { return minDec; }
double Map::getMaxRa() { return maxRa; }
double Map::getMinRa() { return minRa; }
double Map::getCenterRa() { return this->centerRa; }
double Map::getCenterDec() { return this->centerDec; }
double Map::getResolution() { return resolution; }

// GRID GETTERS
int Map::getSize(int i) {
  if (i == 0) {
    return grid.size();
  } else if (i == 1) {
    return grid[i].size();
  } else {
    std::cout << "ERROR: GRID INDEX OUT OF BOUNDS!\n";
    std::cout << "PRESS ENTER TO ABORT\n";
    // throw error here
  }
}
bool Map::getCentroidFlag(int i, int j) { return grid[i][j].getCentroidFlag(); }
double Map::getDec(int i, int j) { return grid[i][j].getDec(); }
double Map::getRa(int i, int j) { return grid[i][j].getRa(); }
double Map::getProcPath(int i, int j) { return grid[i][j].getProcPath(); }

double Map::getSSSProcFlux(int i, int j) { return grid[i][j].getSSSProcFlux(); }
double Map::getSSSScale(int i, int j) { return grid[i][j].getSSSScale(); }
double Map::getSSSCorrelation(int i, int j) {
  return grid[i][j].getSSSCorrelation();
}
double Map::getSSSWeight(int i, int j) { return grid[i][j].getSSSWeight(); }
double Map::getSSSWeight2(int i, int j) { return grid[i][j].getSSSWeight2(); }

double Map::getLSSProcFlux(int i, int j) { return grid[i][j].getLSSProcFlux(); }
double Map::getLSSScale(int i, int j) { return grid[i][j].getLSSScale(); }
double Map::getLSSCorrelation(int i, int j) {
  return grid[i][j].getLSSCorrelation();
}
double Map::getLSSWeight(int i, int j) { return grid[i][j].getLSSWeight(); }
double Map::getLSSWeight2(int i, int j) { return grid[i][j].getLSSWeight2(); }

int Map::getScanNumber(int i, int j) { return grid[i][j].getScanNumber(); }
int Map::getIndexNumber(int i, int j) { return grid[i][j].getIndexNumber(); }

// GRID SETTERS
void Map::setCentroidFlag(int i, int j, bool value) {
  grid[i][j].setCentroidFlag(value);
}
void Map::setDec(int i, int j, double value) { grid[i][j].setDec(value); }
void Map::setRa(int i, int j, double value) { grid[i][j].setRa(value); }
void Map::setLSSFlux(int i, int j, double value) {
  grid[i][j].setLSSFlux(value);
}
void Map::setProcPath(int i, int j, double value) {
  grid[i][j].setProcPath(value);
}

void Map::setSSSProcFlux(int i, int j, double value) {
  grid[i][j].setSSSProcFlux(value);
}
void Map::setSSSScale(int i, int j, double value) {
  grid[i][j].setSSSScale(value);
}
void Map::setSSSCorrelation(int i, int j, double value) {
  grid[i][j].setSSSCorrelation(value);
}
void Map::setSSSWeight(int i, int j, double value) {
  grid[i][j].setSSSWeight(value);
}
void Map::setSSSWeight2(int i, int j, double value) {
  grid[i][j].setSSSWeight2(value);
}

void Map::setLSSProcFlux(int i, int j, double value) {
  grid[i][j].setLSSProcFlux(value);
}
void Map::setLSSScale(int i, int j, double value) {
  grid[i][j].setLSSScale(value);
}
void Map::setLSSCorrelation(int i, int j, double value) {
  grid[i][j].setLSSCorrelation(value);
}
void Map::setLSSWeight(int i, int j, double value) {
  grid[i][j].setLSSWeight(value);
}
void Map::setLSSWeight2(int i, int j, double value) {
  grid[i][j].setLSSWeight2(value);
}

void Map::setScanNumber(int i, int j, int value) {
  return grid[i][j].setScanNumber(value);
}
void Map::setIndexNumber(int i, int j, int value) {
  return grid[i][j].setIndexNumber(value);
}

Map::~Map() {}