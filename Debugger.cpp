#include "Debugger.h"
#include <iostream>
#include <string>

bool debugging = true;
clock_t t_start;

// printing functions
void Debugger::print(std::string type, int value)
{
	if (debugging || type == "Error" || type == "Warn")
	{
		std::cout << type << ": " << value << std::endl;
	}
}

void Debugger::print(std::string type, double value)
{
	if (debugging || type == "Error" || type == "Warn")
	{
		std::cout << type << ": " << value << std::endl;
	}
}
void Debugger::print(std::string type, double value1, double value2)
{
	if (debugging || type == "Error" || type == "Warn")
	{
		std::cout << type << ": " << value1 << '\t' << value2 << std::endl;
	}
}

void Debugger::print(std::string type, std::string message)
{
	if (debugging || type == "Error" || type == "Warn")
	{
		std::cout << type << ": " << message << std::endl;
	}
}
void Debugger::print(std::string type, std::string message, double value)
{
	if (debugging || type == "Error" || type == "Warn")
	{
		std::cout << type << ": " << message << '\t' << value << std::endl;
	}
}

// timing functions
void Debugger::clockTime(bool start)
{
	if (start)
	{
		t_start = clock();
	}
	else
	{
		clock_t t_end;
		t_end = clock();

		if (debugging)
		{
			std::cout << "Total Time\t" << (double)((t_end - t_start) / CLOCKS_PER_SEC) << std::endl;
		}
	}
	return;
}