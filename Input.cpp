#include "Input.h"

Input::Input() {
}

Input::Input(std::string filename) {
	this->filename = filename;
}

void Input::convertFromColumnData(std::vector<std::vector<double>> data) {
	if (data.size() != 11) {
		throw "Expected column data to have 11 columns.";
	}

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

Input::Resolution Input::getDataResolution() {
	return this->resolution;
}

void Input::setDataResolution(Resolution resolution) {
	this->resolution = resolution;
}

int Input::getLowerIntegratedIndex() {
	return this->integratedIndexRange[0];
}

int Input::getHigherIntegratedIndex() {
	return this->integratedIndexRange[1];
}

void Input::setIntegratedIndexRange(std::vector<int> range) {
	this->integratedIndexRange = range;
}
