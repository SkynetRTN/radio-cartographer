#pragma once
#include <time.h>
class timer
{
public:
	timer();
	void clockTime(bool);
	~timer();

private:
	clock_t t;

};

