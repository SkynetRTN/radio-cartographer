#include "io\Observation.h"

Observation::Observation() {
}

Observation::Observation(std::string filename) {
	this->filename = filename;
}

void Observation::convertFromColumnData(std::vector<std::vector<double>> data) {
	if (data.size() != 11)
		throw std::runtime_error("Expected column data to have 11 columns.");

	int dataPointCount = data[0].size();
	
	// Transfer the vector data
	for (int i = 0; i < dataPointCount; ++i) {
		this->times.emplace_back(data[0][i]);
		this->ras.emplace_back(data[1][i]);
		this->decs.emplace_back(data[2][i]);
		this->azimuths.emplace_back(data[3][i]);
		this->elevations.emplace_back(data[4][i]);

		this->lContinuum.emplace_back(data[5][i]);
		this->rContinuum.emplace_back(data[6][i]);

		this->dataDumps.emplace_back(data[7][i]);
		this->calibrations.emplace_back(data[8][i]);
		this->sweepIndices.emplace_back(data[9][i]);
		this->valids.emplace_back(data[10][i]);
	}
}

/**
 * Returns a new Observation object with data that is contained
 * from the provided 'start' index to the end of the observation.
 *
 * @param start: starting index to copy
 * @return new observation object that is a subset of the original
 *
 */
Observation Observation::splitObservation(int start)
{
	return splitObservation(start, this->times.size());
}

/**
 * Returns a new Observation object with data that is contained
 * from the provided 'start' index to the provided 'end' index.
 *
 * @param start: starting index to copy
 * @param end: ending index to copy
 * @return new observation object that is a subset of the original
 *
 */
Observation Observation::splitObservation(int start, int end) {
	Observation obsNew;

	// Copy shared meta data
	obsNew.filename = this->filename;
	obsNew.telescope = this->telescope;
	obsNew.frequencyStep = this->frequencyStep;
	obsNew.coordinate = this->coordinate;
	obsNew.scanType = this->scanType;

	obsNew.observedFrequencies = this->observedFrequencies;
	obsNew.integratedIndexRange = this->integratedIndexRange;
	obsNew.setDataResolution(this->getDataResolution());

	// Map pattern can change for basket weaves
	obsNew.mapPattern = this->observingModes[start];

	// Allocate memory for vector fields
	obsNew.times.resize(end - start);
	obsNew.ras.resize(end - start);
	obsNew.decs.resize(end - start);
	obsNew.azimuths.resize(end - start);
	obsNew.elevations.resize(end - start);
	obsNew.calibrations.resize(end - start);
	obsNew.dataDumps.resize(end - start);
	obsNew.valids.resize(end - start);
	obsNew.sweepIndices.resize(end - start);
	obsNew.lContinuum.resize(end - start);
	obsNew.rContinuum.resize(end - start);

	// Copy vector data
	for (int i = start; i < end; ++i)
	{
		obsNew.times[i - start] = this->times[i];
		obsNew.ras[i - start] = this->ras[i];
		obsNew.decs[i - start] = this->decs[i];
		obsNew.azimuths[i - start] = this->azimuths[i];
		obsNew.elevations[i - start] = this->elevations[i];
		obsNew.calibrations[i - start] = this->calibrations[i];
		obsNew.dataDumps[i - start] = this->dataDumps[i];
		obsNew.valids[i - start] = this->valids[i];
		obsNew.sweepIndices[i - start] = this->sweepIndices[i];

		obsNew.lContinuum[i - start] = this->lContinuum[i];
		obsNew.rContinuum[i - start] = this->rContinuum[i];
	}

	return obsNew;
}

Observation::Resolution Observation::getDataResolution() {
	return this->resolution;
}

void Observation::setDataResolution(Resolution resolution) {
	this->resolution = resolution;
}

int Observation::getLowerIntegratedIndex() {
	return this->integratedIndexRange[0];
}

int Observation::getHigherIntegratedIndex() {
	return this->integratedIndexRange[1];
}

void Observation::setIntegratedIndexRange(std::vector<int> range) {
	this->integratedIndexRange = range;
}
