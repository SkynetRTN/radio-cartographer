#include "Analysis.h"
#include "Tools.h"
#include "RCR.h"
#include "OutputFile.h"
#include "Debugger.h"
#include <iostream>
#include <math.h>
#include <limits>
#include <random>
#include <future>

static double toRad = M_PI / 180.0;
static double toDeg = 180.0 / M_PI;

Analysis::Analysis(Composite composite, MapParameters mParams)
{
	this->rfiScaleBW = mParams.rfiScale;
	this->scans = composite.getScans();
	this->psfFWHM = composite.getPSFFWHM();
	this->compPartSetProc = composite.getCompPartSetProcSSS();
	this->weightScale = mParams.processedWeightScale;
	this->M10Processing = mParams.m10PlusProcessing;
	this->pixelSize = mParams.pixelSize;
}

// Callers
void Analysis::photometry(Map &mapHold, PhotoParams pParams)
{
	// Map
	map = mapHold;

	// Delcarations
	double maxPixelDec, maxPixelRa, centerPixelDec, centerPixelRa;

	// More Declarations
	std::vector<double> coordinates;
	std::vector<double> result, results;
	std::vector<double> coordinatesTemp;
	std::vector<double> coordinatesFull;

	// Define relevant parameters
	int maxRaCoordinate = map.getSize(1);
	determinePixelParameters(maxPixelRa, maxPixelDec, centerPixelDec, centerPixelRa);

	// Output statements
	Debugger::print("Info", "psfFWHM", psfFWHM);
	Debugger::print("Info", "Resolution", map.getResolution());
	Debugger::print("Info", "Max RA Coordinate", maxRaCoordinate);
	Debugger::print("Info", "Map Size", map.getSize(0));

	// Edges and Radii
	Dimensions dim;
	determineDimensions(coordinates, dim, pParams);

	switch (pParams.centroidType)
	{
	case CentroidMethods::CENTER:
		coordinates = autoCentroid(centerPixelDec, centerPixelRa, dim.rSearch);
		results = photometer(0, coordinates, dim);
		break;

	case CentroidMethods::COORDINATES:
		coordinates = pParams.coordinatesPixels;
		
		for (int i = 0; i < coordinates.size() / 2; i++)
		{
			coordinates[2 * i] =  coordinates[2 * i] - 1;
			coordinates[2 * i + 1] = maxRaCoordinate - coordinates[2 * i + 1];
		}

		for (int i = 0; i < coordinates.size() / 2; i++)
		{
			coordinatesTemp = autoCentroid(coordinates[2 * i], coordinates[2 * i + 1], dim.rSearch);
			coordinatesFull.insert(coordinatesFull.end(), coordinatesTemp.begin(), coordinatesTemp.end());
		}

		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			result = photometer(i, coordinatesFull, dim);

			for (int j = 0; j < result.size(); j++)
			{
				results.push_back(result[j]);
			}
		}
		break;

	case CentroidMethods::BRIGHTEST:

		coordinates.push_back(0.0);
		coordinates.push_back(0.0);

		for (int i = 0; i < pParams.numberOfSources; i++)
		{
			coordinatesTemp = autoCentroid(0.0, 0.0, dim.rSearch);
			coordinatesFull.insert(coordinatesFull.end(), coordinatesTemp.begin(), coordinatesTemp.end());
		}

		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			result = photometer(i, coordinatesFull, dim);

			for (int j = 0; j < result.size(); j++)
			{
				results.push_back(result[j]);
			}
		}
		break;
	}

	Output output;
	output.printPhotometry(results);
}
void Analysis::photometryMulti(Map &mapHold, PhotoParams pParams)
{
	// Map
	map = mapHold;

	// Declarations
	double maxPixelDec, maxPixelRa, centerPixelDec, centerPixelRa;

	// More Declarations
	std::vector<double> coordinates;
	std::vector<double> result, results;
	std::vector<double> coordinatesTemp;
	std::vector<double> coordinatesFull;
	std::vector<std::future<std::vector<double>>> futureVec;

	// Define relevant parameters
	int maxRaCoordinate = map.getSize(1);
	determinePixelParameters(maxPixelRa, maxPixelDec, centerPixelDec, centerPixelRa);

	// Output statements
	Debugger::print("Info", "psfFWHM", psfFWHM);
	Debugger::print("Info", "Resolution", map.getResolution());
	Debugger::print("Info", "Max RA Coordinate", maxRaCoordinate);
	Debugger::print("Info", "Map Size", map.getSize(0));

	// Edges and Radii
	Dimensions dim;
	determineDimensions(coordinates, dim, pParams);

	switch (pParams.centroidType)
	{
	case CentroidMethods::CENTER:
		// Perform photometry on source near center
		coordinates = autoCentroid(centerPixelDec, centerPixelRa, dim.rSearch);
		results = photometer(0, coordinates, dim);
		break;

	case CentroidMethods::COORDINATES:
		// Perform photometry on sources with provided pixel locations
		coordinates = pParams.coordinatesPixels;

		for (int i = 0; i < coordinates.size() / 2; i++)
		{
			coordinates[2 * i] = coordinates[2 * i] - 1;
			coordinates[2 * i + 1] = maxRaCoordinate - coordinates[2 * i + 1];
		}

		for (int i = 0; i < coordinates.size() / 2; i++)
		{
			coordinatesTemp = autoCentroid(coordinates[2 * i], coordinates[2 * i + 1], dim.rSearch);
			coordinatesFull.insert(coordinatesFull.end(), coordinatesTemp.begin(), coordinatesTemp.end());
			Debugger::print("Info", "Finished centroiding..");
		}

		futureVec.resize(coordinatesFull.size() / 2);

		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			futureVec[i] = std::async(std::launch::async, &Analysis::photometer, this, i, coordinatesFull, dim);
		}
		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			result = futureVec[i].get();
			results.insert(results.end(), result.begin(), result.end());
		}

		break;

	case CentroidMethods::BRIGHTEST:
		// Perform photometry on brightest source in image
		coordinates.push_back(0.0);
		coordinates.push_back(0.0);

		for (int i = 0; i < pParams.numberOfSources; i++)
		{
			coordinatesTemp = autoCentroid(0.0, 0.0, dim.rSearch);
			coordinatesFull.insert(coordinatesFull.end(), coordinatesTemp.begin(), coordinatesTemp.end());
		}

		futureVec.resize(coordinatesFull.size() / 2);
		
		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			futureVec[i] = std::async(std::launch::async, &Analysis::photometer, this, i, coordinatesFull, dim);
		}
		for (int i = 0; i < coordinatesFull.size() / 2; i++)
		{
			result = futureVec[i].get();
			results.insert(results.end(), result.begin(), result.end());
		}
		break;
	}

	Output output;
	output.printPhotometry(results);
}
std::vector<double> Analysis::photometer(int index, std::vector<double> coordinates, Dimensions dim)
{
	// Declarations
	std::vector<double> metrics, center;
	std::vector<double> annulusParams, apertureParams;
	center.resize(4);

	// Center of centroid (Pixels)
	center[0] = coordinates[2 * index + 1]; // centerRaPix
	center[1] = coordinates[2 * index];     // centerDecPix

	// Center of centroid (Degrees)
	center[2] = map.getRa((int)round(center[1]), (int)round(center[0])); // centRaDeg
	center[3] = map.getDec((int)round(center[1]), (int)round(center[0])) + dim.minDec; // centDecDeg

	// Metrics
	annulusParams = annulusCalculations(dim, center);
	apertureParams = apertureCalculations(dim, center, annulusParams[0]);
	photometryMetrics(dim, annulusParams, apertureParams, metrics);

	// Output statements
	Debugger::print("Info", "Centroid", index);
	Debugger::print("Info", "Dec", center[1]);
	Debugger::print("Info", "Ra", center[0]);
	Debugger::print("Info", "Total Error Bar", metrics[4]);

	// Photometry Results
	std::vector<double> results;

	// These will be different for Skynet integration
	int maxRaCoordinate = map.getSize(1);
	results.push_back((double)maxRaCoordinate - center[0]);
	results.push_back(center[1] + 1);
	results.push_back(annulusParams[0]);
	results.push_back(annulusParams[1]);
	results.push_back(annulusParams[3]);
	results.push_back(apertureParams[0]);
	results.push_back(metrics[0]);
	results.push_back(apertureParams[2]);
	results.push_back(metrics[2]);
	results.push_back(metrics[1]);
	results.push_back(metrics[3]);
	results.push_back(metrics[4]);

	// From Annulus:
	// mu, sigma, wAvg, uncertaintyBG
	// // From Aperture:
	// // px, sourceWSum, mean, peakRadius, peakSum
	// // // From Metrics:
	// // // apertureSigma - metrics[0];
	// // // apertureMu = metrics[1];

	return results;
}

// Determiners
void Analysis::determineDimensions(std::vector<double> coords, Dimensions &dim, PhotoParams pParams)
{
	// Define edge of map
	dim.minRa = map.getMinRa();
	dim.maxRa = map.getMaxRa();
	dim.minDec = map.getMinDec();
	dim.maxDec = map.getMaxDec();

	// Define resolution
	dim.res = map.getResolution();

	// Determine radii
	dim.rIn = pParams.innerRadius / pixelSize;
	dim.rOut = pParams.outerRadius / pixelSize;
	dim.rSearch = Tools::getPythDistance(0, 0, 0, Tools::determinePixel(psfFWHM, 0, map.getResolution()));
}
void Analysis::determinePixelParameters(double &maxPixRa, double &maxPixDec, double &centPixDec, double &centPixRa)
{
	maxPixRa   = Tools::determinePixel(map.getMaxRa(), map.getMinRa(), map.getResolution());
	maxPixDec  = Tools::determinePixel(map.getMaxDec(), map.getMinDec(), map.getResolution());
	centPixDec = Tools::determinePixel((map.getMaxDec() + map.getMinDec()) / 2.0, map.getMinDec(), map.getResolution());
	centPixRa  = Tools::determinePixel((map.getMaxRa() + map.getMinRa()) / 2.0, map.getMinRa(), map.getResolution());// (map.getMaxRa() + map.getMinRa()) / 2.0;
}
double Analysis::determineAnnulusSigma(double mu, double sigma)
{
	double sigma_C;

	if (M10Processing)
	{	// Equation (D.9)
		sigma_C = sigma * log(4.0 * exp(-4.2 * mu / sigma) + exp(1.0));
	}
	else
	{
		sigma_C = sigma; // No correction
	}

	return sigma_C;
}
double Analysis::determineAnnulusMu(std::vector<double> f, std::vector<double> w)
{
	double mean = 0, wSum = 0;

	for (int i = 0; i < f.size(); i++)
	{
		mean += w[i]*f[i];
		wSum += w[i];
	}
	
	return mean / wSum;
}
std::vector<double> Analysis::determineCorrFactor(double peakSum, double peakRadius, double thetaAperPix, double sigma)
{
	// Declarations
	std::vector<double> cVec;
	double zPeak, zPeakGauss, zPeakCos;

	// Re-definitions for clarity
	double thetaRFI = rfiScaleBW;
	double thetaMin = weightScale;
	double thetaPix = pixelSize;

	double thetaAper = thetaAperPix * pixelSize;
	double tApMod = Tools::min(thetaAper * 2, 2.5) / 2.5;

	// Convert peak radius
	peakRadius = peakRadius * pixelSize; // [pixel]*[BW/pixel] = [BW]

	zPeakGauss = peakSum / ((1 / pow(thetaPix, 2))*(M_PI*pow(peakRadius, 2) / 2 + (cos(M_PI*peakRadius) / M_PI) + peakRadius * sin(M_PI*peakRadius) - (1 / M_PI)));
	zPeakCos = peakSum / ((1.13309 - 1.13309*exp(-2.77259*pow(peakRadius, 2))) / pow(thetaPix, 2));
	zPeak = zPeakCos; // Not sure why this is hard-coded \_('.')_/

	// Equation (22)
	double factor = pow(thetaRFI, 2.64*pow(tApMod, -0.54)) + 0.34*pow(tApMod, 0.73) * pow(thetaMin, 1.66*pow(tApMod, -0.18)); 

	// Equation (19)
	double f_lim = 0.22*pow(tApMod, 0.52)*pow(zPeak / (1000.0*sigma), -1.20*pow(tApMod, -0.39))*factor;

	// Equation (18)
	double f_corr = 1.0 + pow((pow(f_lim, -0.82) + 0.052*pow(tApMod, -0.11)*pow(factor, -0.29)), -1.22);

	// Equation (24)
	double uncert_f_lim = 0.082*pow(tApMod, 0.32)*pow(zPeak / (1000.0*sigma), -1.35*pow(tApMod, -0.08))*factor;

    // Equation (23)
	double uncert_corr = pow((pow(uncert_f_lim, -0.82) + 0.065*pow(tApMod, -0.02)*pow(factor, -0.29)), -1.22);

	cVec.push_back(f_corr);
	cVec.push_back(uncert_corr);
	cVec.push_back(2.0); // Correction to total error bar

	return cVec;		
}
std::vector<std::vector<double>> Analysis::determineNRValues(std::vector<double> r, std::vector<double> d, std::vector<double> f, std::vector<double> w, std::vector<bool> flags)
{
	// Declaration
	std::vector<std::vector<double>> nonRejectedValues;
	nonRejectedValues.resize(4);

	// Store non-rejected values
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i] == 1)
		{
			nonRejectedValues[0].push_back(f[i]);
			nonRejectedValues[1].push_back(w[i]);
			nonRejectedValues[2].push_back(d[i]);
			nonRejectedValues[3].push_back(r[i]);
		}
		else // Visualize points that have been rejected
		{
			//map.setProcFlux(d[i], r[i], std::numeric_limits<double>::infinity());
		}
	}

	// nrValues[0] = flux
	// nrValues[1] = weight
	// nrValues[2] = dec
	// nrValues[3] = ra

	return nonRejectedValues;
}

// Annulus and Aperture
std::vector<double> Analysis::annulusCalculations(Dimensions dim, std::vector<double> center)
{
	// Declarations
	int raPix, decPix;
	double w, w2;
	double minRa, maxRa;
	double minDec, maxDec;
	double mean = 0, wSum = 0;
	double radius2, flux, weight;

	// More Declarations
	std::vector<bool> flags;
	std::vector<double> fluxes, weights, decs, ras;
	
	// Define search radius
	double sR = dim.rOut * pixelSize*psfFWHM + 2.0*(dim.res);

	// Define search area
	minRa = Tools::max(center[2] - sR, dim.minRa);
	maxRa = Tools::min(center[2] + sR, dim.maxRa);
	minDec = Tools::max(center[3] - sR, dim.minDec);
	maxDec = Tools::min(center[3] + sR, dim.maxDec);
	
	// Get annulus fluxes, ras, decs, and weights
	for (double i = minDec; i <= maxDec; i += dim.res)
	{
		decPix = Tools::determinePixel(i, dim.minDec, dim.res);

		for (double j = minRa; j <= maxRa; j += dim.res)
		{
			raPix = Tools::determinePixel(j, dim.minRa, dim.res);

			// Distance between point and center
			radius2 = Tools::getPythDistance(center[1], center[0], decPix, raPix);
			flux = map.getSSSProcFlux(decPix, raPix);
			
			if ((radius2 > dim.rIn) && (radius2 < dim.rOut) && !(flux != flux))
			{	// Inside the annulus && flux is a number
				w = map.getSSSWeight(decPix, raPix);
				w2 = map.getSSSWeight2(decPix, raPix);

				weight = (rfiScaleBW != 0.0) ? w2 : w;
				//weight = (w2 != w2) ? w : w2;

				// All annulus values
				weights.push_back(weight);
				fluxes.push_back(flux);
				decs.push_back(decPix);
				ras.push_back(raPix);
			}
		}
	}

	// RCR reject outliers
	RCR rcr = RCR(ES_MODE_DL);
	rcr.performBulkRejection(weights, fluxes);
	flags = rcr.result.flags;

	// Declarations
	double mu, sigma;
	std::vector<double> bgMetrics, params;
	std::vector<std::vector<double>> nrValues;

	// Determine non-rejected values
	nrValues = determineNRValues(ras, decs, fluxes, weights, flags);

	// Equation (D.7)
	mu = determineAnnulusMu(nrValues[0], nrValues[1]);

	// Appendix D calculations
	bgMetrics = sigmaPlus(mu, nrValues);

	// Equation (D.8)
	sigma = determineAnnulusSigma(mu, bgMetrics[0]);

	// Format return vector
	params.resize(4);
	params[0] = mu;			  // Equation (D.7)
	params[1] = sigma;		  // Equation (D.8)
	params[2] = bgMetrics[1]; // Average Weight
	params[3] = bgMetrics[2]; // Equation (D.12)

	return params;
}
std::vector<double> Analysis::apertureCalculations(Dimensions dim, std::vector<double> center, double annulusMu)
{
	// Distance declarations
	double distPix, distCheck, radius2;

	// Dimension declarations
	double raPix, decPix;
	double minRa, maxRa;
	double minDec, maxDec;

	// More declarations
	double px = 0, sourceWSum = 0;
	double flux, weight;
	double peakSum = 0.0;
	double mean = 0.0;

	// Even more declarations
	double w, w2, apCf;
	std::vector<double> decAp, raAp, wAp;

	// Define search radius
	double sR = dim.rIn * pixelSize*psfFWHM + 2.0*(dim.res);
	double peakRadius = Tools::min(dim.rIn, (0.7 / pixelSize));

	// Define search area
	minRa = Tools::max(center[2] - sR, dim.minRa);
	maxRa = Tools::min(center[2] + sR, dim.maxRa);
	minDec = Tools::max(center[3] - sR, dim.minDec);
	maxDec = Tools::min(center[3] + sR, dim.maxDec);

	// Store ra, dec, weights from aperture
	for (double i = minDec; i <= maxDec; i += dim.res)
	{
		decPix = Tools::determinePixel(i, dim.minDec, dim.res);

		for (double j = minRa; j <= maxRa; j += dim.res)
		{
			raPix = Tools::determinePixel(j, dim.minRa, dim.res);
			radius2 = Tools::getPythDistance(center[1], center[0], decPix, raPix);

			flux = map.getSSSProcFlux(decPix, raPix);

			if (radius2 < peakRadius && !(flux != flux))
			{
				peakSum += flux - annulusMu;
			}

			if (radius2 < dim.rIn && !(flux != flux))
			{
				mean += flux - annulusMu;
				px++; // pixel count
				
				// Get appropriate weight value
				w = map.getSSSWeight(decPix, raPix);
				w2 = map.getSSSWeight2(decPix, raPix);

				weight = (rfiScaleBW != 0.0) ? w2 : w;

				// Store aperture values
				decAp.push_back(decPix);
				raAp.push_back(raPix);
				wAp.push_back(weight);
			}
		}
	}

	// Go through points in aperture
	for (int i = 0; i < decAp.size(); i++)
	{
		apCf = 0.0; // APerture Correction Factor
		distCheck = map.getSSSCorrelation(decAp[i], raAp[i]) / (pixelSize * 2.0);

		for (int j = 0; j < decAp.size(); j++)
		{
			distPix = Tools::getPythDistance(decAp[i], raAp[i], decAp[j], raAp[j]);

			if (distPix <= distCheck)
			{
				apCf += 1.0;
			}
		}
		sourceWSum += apCf / wAp[i];
	}

	// Format return vector
	std::vector<double> params;
	params.resize(5);
	params[0] = px;
	params[1] = sourceWSum;
	params[2] = mean;
	params[3] = peakRadius;
	params[4] = peakSum;

	return params;
}

// Calculations
void Analysis::photometryMetrics(Dimensions dim, std::vector<double> anParams, std::vector<double> apParams, std::vector<double> &metrics)
{
	// Equation numbers are as of June 22, 2018

	// Definitions for clarity
	double annulusSigma = anParams[1];
	double annulusWAvg = anParams[2];
	double annulusStdev = anParams[3];

	double px = apParams[0];
	double apertureSourceWSum = apParams[1];
	double apertureMean = apParams[2];
	double peakRadius = apParams[3];
	double peakSum = apParams[4];	

	// Declarations
	double apertureSigma, totalErrorBar;
	double internalErrorBar, externalErrorBar;

	// Equation (D.11)
	apertureSigma = annulusSigma * sqrt(annulusWAvg * apertureSourceWSum);

	// Equation (17) without correction
	internalErrorBar = sqrt(pow(apertureSigma, 2) + pow((px*annulusStdev), 2));

	// Corrections
	std::vector<double> cV;
	double apertureMu_C, internalErrorBar_C;

	// '_C' represents a corrected value
	cV = determineCorrFactor(peakSum, peakRadius, dim.rIn, annulusSigma);

	// Apply corrections
	apertureMu_C = apertureMean * cV[0];
	externalErrorBar = apertureMean * cV[1];		// Equation (23)
	internalErrorBar_C = internalErrorBar * cV[2];  // Equation (17)

	// Total Error Bar
	totalErrorBar = std::sqrt(std::pow(externalErrorBar, 2) + std::pow(internalErrorBar_C, 2));

	// Define return vector
	metrics.resize(5);
	metrics[0] = apertureSigma;      // Equation (D.11)
	metrics[1] = apertureMu_C;
	metrics[2] = internalErrorBar_C; // Equation (17)
	metrics[3] = externalErrorBar;	 // Equation (23)
	metrics[4] = totalErrorBar;
}
std::vector<double> Analysis::sigmaPlus(double mu, std::vector<std::vector<double>> nrValues)
{
	// Declarations and Definitions
	double distPix, correlationScale, nCnt;
	double top = 0, wSum = 0, sigmaSum = 0.0, wCnt = 0;

	// Non-rejected flux, weights, decs (px), ras (px)
	std::vector<double> z = nrValues[0];
	std::vector<double> w = nrValues[1];
	std::vector<double> decPix = nrValues[2];
	std::vector<double> raPix  = nrValues[3];

	// Sum over non-rejected, z > mu values
	for (int i = 0; i < w.size(); i++)
	{	
		if ((z[i] - mu) > 0.0)
		{	
			top += w[i] * (z[i] - mu)*(z[i] - mu);
			wSum += w[i];
			wCnt += 1.0;
			nCnt = 0.0;

			correlationScale = map.getSSSCorrelation(decPix[i], raPix[i]) / (pixelSize * 2.0);

			for (int j = 0; j < w.size(); j++)
			{	// Check within correlation scale	
				distPix = Tools::getPythDistance(decPix[i], raPix[i], decPix[j], raPix[j]);

				if (distPix <= correlationScale)
				{
					nCnt += 1.0; // Correction in Equation (D.12)
				}
			}
			sigmaSum += nCnt * w[i] * w[i] * (z[i] - mu)*(z[i] - mu);
		}
	}

	// Format return vector
	std::vector<double> container;
	container.resize(3);

	container[0] = sqrt(top / wSum); // Equation (D.8) without correction
	container[1] = wSum / wCnt;	     // Average weight
	container[2] = sqrt(sigmaSum) / wSum;  // Equation (D.12) // Uncertainty in BG

	return container;
}

// Centroids
std::vector<double> Analysis::findBrightestPixel()
{
	int maxDec, maxRa;
	double maxFlux = -999999, centroidRa = -999999, centroidDec = -999999;
	std::vector<double> results;
	results.resize(2);
	for (int i = 0; i < map.getSize(0); i++)
	{
		for (int j = 0; j < map.getSize(1); j++)
		{
			if (map.getSSSProcFlux(i, j) > maxFlux && map.getCentroidFlag(i,j) == false)
			{
				maxFlux = map.getSSSProcFlux(i, j);
				maxDec = i;
				maxRa = j;
			}
		}
	}
	results[0] = maxDec;
	results[1] = maxRa;
	return results;
}
std::vector<double> Analysis::autoCentroid(int yPixel, int xPixel, double searchRadius)
{
	//yPixel is dec pixel
	//xPixel is ra pixel
	//maps are called using (decPixel, raPixel)
	double maxFlux = -999999, centroidRa = -999999, centroidDec = -999999;
	std::vector<double> coordinates, inRange;
	double distancePixel, distanceDeg;
	double initialDec, initialRa;

	centroidDec = yPixel;
	centroidRa = xPixel;

	if (xPixel == 0.0 && yPixel == 0.0)
	{
		coordinates = findBrightestPixel();
		centroidDec = coordinates[0];
		centroidRa = coordinates[1];
	}

	for (int k = 0; k < 5; k++)
	{
		for (int i = 0; i < map.getSize(0); i++)
		{
			for (int j = 0; j < map.getSize(1); j++)
			{
				distancePixel = Tools::getPythDistance(centroidDec, centroidRa, i, j);
				distanceDeg = Tools::getGCDistance(map.getDec(centroidDec, centroidRa), map.getRa(centroidDec, centroidRa), map.getDec(i, j), map.getRa(i, j), map.getCenterDec())*toDeg;

				if (distancePixel < searchRadius)
				{
					if (map.getSSSProcFlux(i, j) != map.getSSSProcFlux(i, j))
					{
						continue;
					}
					else
					{
						map.setCentroidFlag(i, j, true);
						inRange.push_back(map.getSSSProcFlux(i, j));
						inRange.push_back(i);
						inRange.push_back(j);
						inRange.push_back(cos((M_PI*distanceDeg) / (psfFWHM*2.0)));
					}
				}
			}
		}

		coordinates = determineCenters(inRange.size() / 4, inRange, map.getSSSScale(centroidDec, centroidRa), (double)centroidDec, (double)centroidRa);

		centroidDec = (int)std::round(coordinates[0]); // These get rounded to ints later.
		centroidRa = (int)std::round(coordinates[1]);// ERROR: ROUNDING CHECK
	}
	
	return coordinates;
}
std::vector<double> Analysis::determineCenters(int pointCount, std::vector<double>& inRange, double scale, double i_0, double j_0)
{
	int columnCount = 6, pivotIndex;
	double pivot, swap, multiplier, toRad = M_PI / 180.0;
	double weightingFunctionTemp, weightPower;
	double w = 0, xw = 0, yw = 0, xxw = 0, yyw = 0, xyw = 0,
		xxxw = 0, xyyw = 0, xxyw = 0, yyyw = 0,
		xxxxw = 0, xxyyw = 0, xxxyw = 0, yyyyw = 0, xyyyw = 0,
		zw = 0, xzw = 0, yzw = 0, xxzw = 0, yyzw = 0, xyzw = 0,
		xHold, yHold, zHold, wHold, finalAnswer;
	double ww = 0, sumWSquared = 0;
	weightingFunctionTemp = scale;
	double alpha = (log(0.5)) / (log(cos(M_PI*weightingFunctionTemp / 4.0)));

	std::vector<double> A;
	std::vector<double> b;
	std::vector<double> coef;
	std::vector<double> coordinates;
	A.resize(columnCount * columnCount, 0.0);
	b.resize(columnCount, 0.0);
	for (int i = 0; i < pointCount; i++)
	{

		//CHECK THAT THESE MATCH UP
		zHold = inRange[4 * i];
		xHold = (inRange[4 * i + 1] - i_0);
		yHold = (inRange[4 * i + 2] - j_0);

		wHold = pow(inRange[4 * i + 3], alpha);
		ww += std::abs(pow(wHold, 2));
		w += wHold;
		xw += xHold*wHold;
		yw += yHold*wHold;
		xxw += xHold*xHold*wHold;
		yyw += yHold*yHold*wHold;
		xyw += xHold*yHold*wHold;
		xxxw += xHold*xHold*xHold*wHold;
		xyyw += xHold*yHold*yHold*wHold;
		xxyw += xHold*xHold*yHold*wHold;
		yyyw += yHold*yHold*yHold*wHold;
		xxxxw += xHold*xHold*xHold*xHold*wHold;
		xxyyw += xHold*xHold*yHold*yHold*wHold;
		xxxyw += xHold*yHold*xHold*xHold*wHold;
		yyyyw += yHold*yHold*yHold*yHold*wHold;
		xyyyw += xHold*yHold*yHold*yHold*wHold;
		zw += zHold*wHold;
		xzw += xHold*zHold*wHold;
		yzw += yHold*zHold*wHold;
		xxzw += xHold*xHold*zHold*wHold;
		yyzw += yHold*yHold*zHold*wHold;
		xyzw += xHold*yHold*zHold*wHold;
	}

	sumWSquared = std::abs(pow(w, 2));

	A[0] = w; A[1] = xw; A[2] = yw; A[3] = xxw; A[4] = yyw; A[5] = xyw;
	A[6] = xw;  A[7] = xxw; A[8] = xyw; A[9] = xxxw; A[10] = xyyw; A[11] = xxyw;
	A[12] = yw; A[13] = xyw; A[14] = yyw; A[15] = xxyw; A[16] = yyyw; A[17] = xyyw;
	A[18] = xxw; A[19] = xxxw; A[20] = xxyw; A[21] = xxxxw; A[22] = xxyyw; A[23] = xxxyw;
	A[24] = yyw; A[25] = xyyw; A[26] = yyyw; A[27] = xxyyw; A[28] = yyyyw; A[29] = xyyyw;
	A[30] = xyw; A[31] = xxyw; A[32] = xyyw; A[33] = xxxyw; A[34] = xyyyw; A[35] = xxyyw;
	b[0] = zw; b[1] = xzw; b[2] = yzw; b[3] = xxzw; b[4] = yyzw; b[5] = xyzw;

	coef = Tools::matrixSolver(columnCount, A, b);

	double x = ((coef[2] * coef[5]) - (2 * coef[1] * coef[4])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));
	double y = ((coef[1] * coef[5]) - (2 * coef[2] * coef[3])) / ((4 * coef[3] * coef[4]) - pow(coef[5], 2));

	coordinates.push_back(i_0 + x); //change in dec
	coordinates.push_back(j_0 + y); //change in ra

	if (Tools::getGCDistance(coordinates[0], coordinates[1], i_0, j_0, compPartSetProc.centerDecDeg)*toDeg > psfFWHM/pixelSize)
	{
		Debugger::print("Warn", "Centroiding unsuccessful; forcing original coordinates.");
		coordinates[0] = i_0;
		coordinates[1] = j_0;
	}

	return coordinates;
}


Analysis::~Analysis()
{
}	