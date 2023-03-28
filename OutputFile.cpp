#include "OutputFile.h"
#include <math.h>
#include <iomanip>
#include <time.h>
#include <sstream>
#include <string>

Output::Output()
{

}

void Output::printTelescopeInfo(SurveyParameters sParams)
{
	double psfFWHM;
	std::string telescopeName;

	outputFile.open("Output.txt", std::ofstream::out | std::ofstream::trunc);
	if (sParams.tele == 0)
	{
		telescopeName = "40 foot";
		psfFWHM = 1.22*299792458.0*180.0 / (sParams.frequency* pow(10, 9) * 12.192 * M_PI);
	}
	else if (sParams.tele == 1)
	{
		telescopeName = "20 meter";
		psfFWHM = 1.22*299792458.0*180.0 / (sParams.frequency* pow(10, 9) * 20.0 * M_PI);
	}
	else
	{
		telescopeName = "105 meter";
		psfFWHM = 1.22*299792458.0*180.0 / (sParams.frequency* pow(10, 9) * 105.0 * M_PI);
	}

	outputFile << "TELESCOP," << telescopeName << ",Telescope name\n";
	outputFile << "FREQ," << sParams.frequency << ",Frequency (GHz)\n";
	outputFile << "BEAM," << psfFWHM << ",Beam (deg)\n";
	outputFile.close();

	parameters.push_back(telescopeName);
	parameters.push_back(std::to_string(sParams.frequency));
	parameters.push_back(std::to_string(psfFWHM));
}
void Output::printProcessingHeader()
{
	//outputFile.open("Output.txt", std::ios_base::app);
	//outputFile << "image processing:\n";
	//outputFile << "background subtraction (used = 1), scale (beams) \n";
	//outputFile << "time-delay correction (used = 1), correction (s) \n";
	//outputFile << "rfi cleaning (used = 1), scale (beams) \n";
	//outputFile << "surface modeling minimum scale (beams) \n\n";
	//outputFile.close();
}
void Output::printCalibrationHeader()
{
	
	outputFile.open("Output.txt", std::ios_base::app);

	outputFile << "polarization channel and gain calibration:\n";
	outputFile << "left channel (used = 1), pre-cal (used = 1), t (s), cal (V), cal (Jy) \n";
	outputFile << "left channel (used = 1), post-cal (used = 1), t (s), cal (V), cal (Jy) \n";
	outputFile << "right channel (used = 1), pre-cal (used = 1), t (s), cal (V), cal (Jy) \n";
	outputFile << "right channel (used = 1), post-cal (used = 1), t (s), cal (V), cal (Jy) \n";

	outputFile.close();
}
void Output::printWScale(double wScale)
{    
	outputFile.open("Output.txt", std::ios_base::app);
	outputFile << "SMSCALE," << wScale << ",Minimum surface modeling scale (beams)\n";
	outputFile.close();

	parameters.push_back(std::to_string(wScale));
}

void Output::printPhotometryHolder()
{
	outputFile.open("Output.txt", std::ios_base::app);

	outputFile << "(Photometry=off)\n";
	outputFile << "(Photometry=off)";

	outputFile.close();
}
void Output::printPhotometry(std::vector<double> results)
{
	int num = 12;
	
	outputFile.open("Output.txt", std::ios_base::app);
	outputFile << "RAPIX\t" << "DECPIX\t" << "Mu   \t" << "Sigma\t\t" << "Uncrtnty_in_Mu\t" << "Aper_Pixels\t" << "Aper_Sigma\t" << "Mesred_Phtmtry\t" << "Internal_Err\t" << "Correctd_Phtmtry\t" << "Uncrtnty_in_Ap_w_Corr\t" << "Total_Error_Bar\t" << "\n";

	for (int i = 0; i < results.size() / num; i++)
	{
		outputFile << std::setprecision(6) << (double)results[num * i] << "\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 1] << "\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 2] << "\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 3] << "\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 4] << "\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 5] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 6] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 7] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 8] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 9] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 10] << "\t\t";
		outputFile << std::setprecision(6) << (double)results[num * i + 11] << "\n";
	}
	outputFile.close();
}
void Output::printBgSubtractionInfo(double bgScale)
{
	
	outputFile.open("Output.txt", std::ios_base::app);

	outputFile << "BGSCALE," << bgScale << ",Background subtraction scale (beams)\n";

	outputFile.close();
}
void Output::printTimeShiftInfo(bool ts, double t_int)
{
	outputFile.open("Output.txt", std::ios_base::app);

	if (ts)
	{
		outputFile << "TIMESHIF," << t_int << ",Time delay correction (s) (Used=true)\n";
	}
	else
	{
		outputFile << "TIMESHIF," << t_int << ",Time delay correction (s) (Used=false)\n";
	}

	outputFile.close();
}
void Output::printRfiInfo(double rfiScale)
{
	
	outputFile.open("Output.txt", std::ios_base::app);

	outputFile << "RFISCALE," << rfiScale << ",RFI removal scale (beams)\n";

	outputFile.close();
}


Output::~Output()
{

};
