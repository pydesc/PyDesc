/*
 * Copyright 2008-2010, 2017 Pawel Daniluk
 * 
 * This file is part of PyDesc.
 * 
 * PyDesc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PyDesc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.
 */
 
#include "timeout.h"
#include <time.h>
#include <stdio.h>

#ifdef  __APPLE__
#include<mach/mach_time.h>
#include<mach/clock.h>
#else
#include<sys/time.h>
#include<sys/resource.h>
#endif

#include"simple_macros.h"

int start_time=0;
int timeout=0;

long long int timer_mark=0;

void timer_start()
{

	printf("timer_start\n");
#ifdef  __APPLE__
	timer_mark=mach_absolute_time();
#else
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	timer_mark=(usage.ru_utime.tv_sec+usage.ru_stime.tv_sec)*1000+(usage.ru_utime.tv_usec+usage.ru_stime.tv_usec)/1000;
#endif
}

long long int get_timer()
{
#ifdef  __APPLE__
	uint64_t        elapsed;
	uint64_t     elapsedNano;
	// Calculate the duration.

	static mach_timebase_info_data_t    sTimebaseInfo;

	elapsed = mach_absolute_time() - timer_mark;

	// Convert to nanoseconds.

	// Have to do some pointer fun because AbsoluteToNanoseconds
	// works in terms of UnsignedWide, which is a structure rather
	// than a proper 64-bit integer.

	(void) mach_timebase_info(&sTimebaseInfo);

    // Do the maths.  We hope that the multiplication doesn't
    // overflow; the price you pay for working in fixed point.

    elapsedNano = elapsed * sTimebaseInfo.numer / sTimebaseInfo.denom;

	long long int res=elapsedNano/1000000;
#else 

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	long long int res=(usage.ru_utime.tv_sec+usage.ru_stime.tv_sec)*1000+(usage.ru_utime.tv_usec+usage.ru_stime.tv_usec)/1000;

	res-=timer_mark;


#endif

	return res;
}


void set_timeout(int n_sec) 
{
	timeout=n_sec;

	start_time=time(0);
}


int check_timeout()
{
	if(!timeout) return 0;

	int curr_time=time(0);

	if(curr_time>=start_time+timeout) return 1;

	return 0;
}

