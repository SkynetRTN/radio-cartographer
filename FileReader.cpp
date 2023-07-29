#include "FileReader.h"
#include "FourtyParser.h"
#include "FitsParser.h"
#include "GBParser.h"

#include <fstream>
#include <iostream>


/**
 * Supported file formats include '.md2', '.txt' and '.fits'.
 */
 

FileReader::FileReader(std::string filename) {
	this->filename = filename;

	setFileType(filename);

	if (!exists(filename)) {
		std::cerr << "No input files were provided.\n";
		std::exit(EXIT_FAILURE);
	}
}

FileReader::FileReader(std::vector<std::string> filenames) {
	this->filenames = filenames;

	if (filenames.size() < 1) {
		std::cerr << "No input files were provided.\n";
		std::exit(EXIT_FAILURE);
	}

	setFileType(filenames[0]);

	for (const std::string filename : filenames) {
		if (!exists(filename)) {
			std::cerr << filename + " does not exist.\n";
			std::exit(EXIT_FAILURE);
		}
	}

	// If providing one forty foot input, use other constructor method
	if (fileType == FileType::MD2 && filenames.size() != 2) {
		std::cerr << "Expected two input files but received ";
		std::cerr << std::to_string(filenames.size()) + ".\n";
		std::exit(EXIT_FAILURE);
	}
}

/**
 * Reads the input file(s) formatting it for Survey
 *
 * @return 2D vector containing relevant data
 */
std::vector<std::vector<double>> FileReader::read() {
	if (fileType == FileType::FITS) {
		FitsParser parser = FitsParser(filename);
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
	else if (filename.find(".fits") != std::string::npos) {
		fileType = FileType::FITS;
	}
	else {
		std::cerr << "Unsupported file type provided: " + filename + "\n";
		std::exit(EXIT_FAILURE);
	}
}

FileReader::~FileReader() {

}