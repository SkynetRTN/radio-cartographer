#include "map\Pixel.h"


Pixel::Pixel()
{
	this->SSSProcFlux = 0.0;
	this->SSSScale = 0.0;
	this->SSSWeight = 0.0;
	this->SSSWeight2 = 0.0;
	this->SSSCorrelation = 0.0;
	
	this->LSSFlux = 0.0;
	this->procPath = -1.0;
	this->centroidFlag = false;
}


//getters
double Pixel::getSSSProcFlux()
{
	return SSSProcFlux;
}
double Pixel::getSSSScale()
{
	return SSSScale;
}
double Pixel::getSSSCorrelation()
{
	return SSSCorrelation;
}
double Pixel::getSSSWeight()
{
	return SSSWeight;
}
double Pixel::getSSSWeight2()
{
	return SSSWeight2;
}

double Pixel::getLSSProcFlux()
{
	return LSSProcFlux;
}
double Pixel::getLSSScale()
{
	return LSSScale;
}
double Pixel::getLSSCorrelation()
{
	return LSSCorrelation;
}
double Pixel::getLSSWeight()
{
	return LSSWeight;
}
double Pixel::getLSSWeight2()
{
	return LSSWeight2;
}

//setters
double Pixel::getLSSFlux()
{
	return LSSFlux;
}
double Pixel::getProcPath()
{
	return procPath;
}
bool Pixel::getCentroidFlag()
{
	return centroidFlag;
}
double Pixel::getDec()
{
	return dec;
}
double Pixel::getRa()
{
	return ra;
}
int Pixel::getScanNumber()
{
	return this->scanNumber;
}
int Pixel::getIndexNumber()
{
	return this->indexNumber;
}


void Pixel::setSSSProcFlux(double value)
{
	this->SSSProcFlux = value;
}
void Pixel::setSSSScale(double value)
{
	this->SSSScale = value;
}
void Pixel::setSSSCorrelation(double value)
{
	this->SSSCorrelation = value;
}
void Pixel::setSSSWeight(double value)
{
	this->SSSWeight = value;
}
void Pixel::setSSSWeight2(double value)
{
	this->SSSWeight2 = value;
}

void Pixel::setLSSProcFlux(double value)
{
	this->LSSProcFlux = value;
}
void Pixel::setLSSScale(double value)
{
	this->LSSScale = value;
}
void Pixel::setLSSCorrelation(double value)
{
	this->LSSCorrelation = value;
}
void Pixel::setLSSWeight(double value)
{
	this->LSSWeight = value;
}
void Pixel::setLSSWeight2(double value)
{
	this->LSSWeight2 = value;
}

void Pixel::setLSSFlux(double value)
{
	this->LSSFlux = value;
}
void Pixel::setProcPath(double value)
{
	this->procPath = value;
}
void Pixel::setCentroidFlag(bool val)
{
	this->centroidFlag = val;
}
void Pixel::setDec(double value)
{
	this->dec = value;
}
void Pixel::setRa(double value)
{
	this->ra = value;
}
void Pixel::setScanNumber(int value)
{
	this->scanNumber = value;
}
void Pixel::setIndexNumber(int value)
{
	this->indexNumber = value;
}

Pixel::~Pixel()
{
}