#pragma once
#include <string>
#include <vector>
#include "io\Observation.h"

class FileReader {
public:
	FileReader(std::string);
	FileReader(std::vector<std::string>);
	FileReader(std::string, bool);

	enum class FileType {
		MD2,
		TEXT,
		FITS_20M,
		FITS_GBT
	};

	Observation read();
	static Observation merge(std::vector<Observation>&);

	FileType getFileType();

	~FileReader();

private:
	bool logToConsole;
	std::string filename = "";
	std::vector<std::string> filenames;
	FileType fileType;

	bool exists(std::string);
	void setFileType(std::string);
};