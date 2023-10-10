#include "FileReader.h"
#include "FourtyParser.h"
#include "FitsParser.h"
#include "GBTParser.h"
#include "GBParser.h"

#include <fstream>
#include <iostream>


/**
 * Supported file formats include '.md2', '.txt', '.dcr.fits', and '.cyb.fits'.
 * Supported telescopes include the the forty foot, twenty meter, and the Green Bank Telescope
 */
 

FileReader::FileReader(std::string filename) {
	this->filename = filename;

	setFileType(filename);

	if (!exists(filename)) {
		throw "No input files were provided.\n";
	}
}

FileReader::FileReader(std::vector<std::string> filenames) {
	this->filenames = filenames;

	if (filenames.size() < 1) {
		throw "No input files were provided.\n";
	}

	setFileType(filenames[0]);

	for (const std::string filename : filenames) {
		if (!exists(filename)) {
			throw filename + " does not exist.\n";
		}
	}

	// If providing one forty foot input, use other constructor method
	if (fileType == FileType::MD2 && filenames.size() != 2) {
		throw "Expected two input files";
	}
}

/**
 * Reads the input file(s) formatting it for Survey
 *
 * @return 2D vector containing relevant data
 */
Input FileReader::read() {
	if (fileType == FileType::FITS_20M) {
		FitsParser parser = FitsParser(filename);
		return parser.parse();
	}
	else if (fileType == FileType::FITS_GBT) {
		GBTParser parser = GBTParser(filename);
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
 * Determines if the file exists.
 *
 * @return true if the file exists, else false.
 */
bool FileReader::exists(std::string filename) {
	std::string mode = "rb"; // opens non-text file formats

	if (fileType == FileType::TEXT) mode = "r";

	FILE* file = std::fopen(filename.c_str(), mode.c_str());

	if (!file) {
		std::cerr << "File " + filename + " does not exist.\n";
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
		throw "Unsupported file type provided: " + filename + "\n";
	}
}

FileReader::FileType FileReader::getFileType() {
	return this->fileType;
}

FileReader::~FileReader() {

}