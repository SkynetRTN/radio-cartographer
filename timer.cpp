#include "timer.h"
#include <iostream>



timer::timer()
{
}

void timer::clockTime(bool start)
{
	if (start)
	{
		t = clock();
	}
	else
	{
		clock_t tEnd;
		tEnd = clock();
		std::cout << "Total Time\t" << (double)((t - tEnd) / CLOCKS_PER_SEC) << "\n";
	}
	return;
}
timer::~timer()
{
}
