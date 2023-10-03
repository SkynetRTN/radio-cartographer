#pragma once
#include <string>
#include <vector>
#include "Input.h"

class FileReader {
public:
	FileReader(std::string);
	FileReader(std::vector<std::string>);

	enum class FileType {
		MD2,
		TEXT,
		FITS
	};

	Input read();
	FileType getFileType();

	~FileReader();

private:

	std::string filename = "";
	std::vector<std::string> filenames;
	FileType fileType;

	bool exists(std::string);
	void setFileType(std::string);
};