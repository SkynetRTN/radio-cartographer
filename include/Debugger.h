#ifndef DEBUGGER_H
#define DEBUGGER_H

#include <iostream>
#include <fstream>
#include <time.h>

namespace Debugger
{
	// printing
	void print(std::string, int);

	void print(std::string, double);
	void print(std::string, double, double);	
	
	void print(std::string, std::string);
	void print(std::string, std::string, double);

	// timing
	void clockTime(bool);
}


#endif
