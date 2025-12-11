class Pixel
{
public:
	Pixel();

	//setters
	void setSSSProcFlux(double);
	void setSSSScale(double);
	void setSSSCorrelation(double); //DYLAN
	void setSSSWeight(double);
	void setSSSWeight2(double);

	void setLSSProcFlux(double);
	void setLSSScale(double);
	void setLSSCorrelation(double); //DYLAN
	void setLSSWeight(double);
	void setLSSWeight2(double);


	void setLSSFlux(double);
	void setProcPath(double);
	void setCentroidFlag(bool);
	void setRa(double);
	void setDec(double);
	void setScanNumber(int);
	void setIndexNumber(int);


	//getters
	double getSSSProcFlux();
	double getSSSScale();
	double getSSSCorrelation(); //DYLAN
	double getSSSWeight();
	double getSSSWeight2();

	double getLSSProcFlux();
	double getLSSScale();
	double getLSSCorrelation(); //DYLAN
	double getLSSWeight();
	double getLSSWeight2();

	double getProcPath();
	double getLSSFlux();
	double getRa();
	double getDec();
	int getScanNumber();
	int getIndexNumber();

	bool getCentroidFlag();


	~Pixel();

private:

	double SSSProcFlux;
	double SSSScale;
	double SSSCorrelation;//DYLAN
	double SSSWeight;
	double SSSWeight2;

	double LSSProcFlux;
	double LSSScale;
	double LSSCorrelation;//DYLAN
	double LSSWeight;
	double LSSWeight2;

	double LSSFlux;
	double procPath;
	double ra;
	double dec;

	int scanNumber;
	int indexNumber;
	bool centroidFlag;

};

