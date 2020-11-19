/*
 * mytimer.h
 *
 *  Created on: 29 giu 2017
 *      Author: anonym
 */

#ifndef MYTIMER_H_
#define MYTIMER_H_


#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
class mytimer{

public:
	mytimer() {
		start = timestamp(); // postcondition: elapsed()==0
		pauseSum=0;
  	}
	void restart() {
		start = timestamp(); // post: elapsed()==0
		pauseSum=0;
  	}
	void pause(){
		wpause=timestamp();
	}
	void resume(){
		pauseSum += (timestamp()- wpause);
	}
  	double elapsed(){                  // return elapsed time in seconds
		return timestamp()-start - pauseSum;
	}

private:
	double start;
	double wpause;
	double pauseSum;

	inline double timestamp() {


		//Precise timing
		struct timespec time;
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time);
		return (double)time.tv_sec + (double)time.tv_nsec * 1e-9;


		// Old timing
		//	struct timeval tp;
		//	gettimeofday(&tp, NULL);
		//	return double(tp.tv_sec) + tp.tv_usec / 1000000.;

	}

};







#endif /* MYTIMER_H_ */
