#pragma once
#include <string>
#include <vector>

class FileReader {
public:
	FileReader(std::string);
	FileReader(std::vector<std::string>);

	enum class FileType {
		MD2,
		TEXT,
		FITS
	};

private:

	std::string filename = "";
	std::vector<std::string> filenames;
	FileType fileType;

	std::vector<std::vector<double>> read();

	bool exists(std::string);
	void setFileType(std::string);

protected:
};