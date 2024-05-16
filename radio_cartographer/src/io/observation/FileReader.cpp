#include "io\observation\FileReader.h"
#include "io\observation\FortyFootParser.h"
#include "io\observation\TwentyMeterParser.h"
#include "io\observation\HundredMeterParser.h"
#include "io\observation\LegacyTwentyParser.h"
#include "utils\Tools.h"

#include <fstream>
#include <iostream>


/**
 * Supported file formats include '.md2', '.txt', '.dcr.fits', and '.cyb.fits'.
 * Supported telescopes include the the forty-foot, twenty-meter, and the Green Bank Telescope.
 */
 
FileReader::FileReader(std::string filename) : FileReader(filename, false) {}

FileReader::FileReader(std::string filename, bool logToConsole) {
	if (!exists(filename)) {
		throw std::runtime_error("[ File Reader ] The provided file does not exsit: " + filename);
	}
	setFileType(filename);
	this->filename = filename;
	this->logToConsole = false;
}

/**
 * Constructor for forty-foot input that uses two .md2 files.
 */
FileReader::FileReader(std::vector<std::string> filenames) {
	this->filenames = filenames;

	if (filenames.size() < 1) {
		throw std::runtime_error("[ File Reader] No input files were provided.");
	}

	setFileType(filenames[0]);

	for (const std::string filename : filenames) {
		if (!exists(filename)) {
			throw std::runtime_error("[ File Reader ]" + filename + " does not exist.");
		}
	}

	// If providing one forty-foot input, use other constructor method
	if (fileType == FileType::MD2 && filenames.size() != 2) {
		throw std::runtime_error("[ File Reader ] Expected two input files");
	}
}

/**
 * Reads the input file(s) formatting it for Survey
 */
Observation FileReader::read() {
	if (fileType == FileType::FITS_20M) {
		TwentyMeterParser parser = TwentyMeterParser(filename);
		return parser.parse();
	}
	else if (fileType == FileType::FITS_GBT) {
		HundredMeterParser parser = HundredMeterParser(filename);
		return parser.parse();
	}
	else if (fileType == FileType::TEXT) {
		GBParser parser = GBParser(filename);
		return parser.parseFile();
	}
	else if (fileType == FileType::MD2) {
		if (filename != "") {
			FourtyParser parser = FourtyParser(filename);
			return parser.parseFile();
		}
		else if (filenames.size() == 2) {
			FourtyParser parser = FourtyParser(filenames[0], filenames[1]);
			return parser.parseFile();
		}
	}
}

/**
 * Merges multiple Observation objects to a single object.
 */
Observation FileReader::merge(std::vector<Observation> &observations)
{
	if (observations.size() == 0) {
		throw std::runtime_error("[ File Reader ] Cannot merge an empty vector of observations.");
	}
	else if (observations.size() == 1) {
		return observations[0];
	}

	Observation merged = observations[0];

	for (int i = 1; i < observations.size(); ++i) {
		if (observations[i].getDataResolution() != observations[i - 1].getDataResolution())
			throw std::runtime_error("[ File Reader ] Cannot merge observations with different data resolutions.");

		merged.filename += " + " + observations[i].filename;

		for (int j = 0; j < observations[i].times.size(); ++j)
		{
			merged.times.emplace_back(observations[i].times[j]);
			merged.ras.emplace_back(observations[i].ras[j]);
			merged.decs.emplace_back(observations[i].decs[j]);
			merged.azimuths.emplace_back(observations[i].azimuths[j]);
			merged.elevations.emplace_back(observations[i].elevations[j]);
			merged.dataDumps.emplace_back(observations[i].dataDumps[j]);
			merged.valids.emplace_back(observations[i].valids[j]);
			merged.sweepIndices.emplace_back(observations[i].sweepIndices[j]);

			if (observations[0].getDataResolution() == Observation::Resolution::LOW)
			{
				merged.lSpectra.emplace_back(observations[i].lSpectra[j]);
				merged.rSpectra.emplace_back(observations[i].rSpectra[j]);
			}
			else
			{
				merged.lLowFrequencySpectra.emplace_back(observations[i].lLowFrequencySpectra[j]);
				merged.rLowFrequencySpectra.emplace_back(observations[i].rLowFrequencySpectra[j]);
				merged.lHighFrequencySpectra.emplace_back(observations[i].lHighFrequencySpectra[j]);
				merged.rHighFrequencySpectra.emplace_back(observations[i].rHighFrequencySpectra[j]);
			}
		}
	}
	return merged;
}

/**
 * Determines if the file exists.
 *
 * @return true if the file exists, else false.
 */
bool FileReader::exists(std::string filename) {
	std::string mode = "rb"; // opens non-text file formats

	if (fileType == FileType::TEXT) mode = "r";

	FILE* file = std::fopen(filename.c_str(), mode.c_str());

	if (!file) {
		return false;
	}
	std::fclose(file);
	
	return true;
}

void FileReader::setFileType(std::string filename) {
	if (filename.find(".txt") != std::string::npos) {
		fileType = FileType::TEXT;
	}
	else if (filename.find(".md2") != std::string::npos) {
		fileType = FileType::MD2;
	}
	else if (filename.find(".dcr.fits") != std::string::npos) {
		fileType = FileType::FITS_GBT;
	}
	else if (filename.find(".fits") != std::string::npos) {
		fileType = FileType::FITS_20M;
	}
	else {
		throw std::runtime_error("[ File Reader ] Unsupported file type provided: " + filename);
	}
}

FileReader::FileType FileReader::getFileType() {
	return this->fileType;
}

FileReader::~FileReader() {

}