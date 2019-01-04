#ifndef TIMER
#define TIMER

#include <chrono>
#include <iostream>

//Timer for timing seconds, returned as double.

class Timer{
public:
    Timer(void){}
    inline void start(void){
	t1 = std::chrono::high_resolution_clock::now();
    }
    inline void stop(void){
	t2 = std::chrono::high_resolution_clock::now();
    }
    double elapsed(){
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	return time_span.count();
    }
private:
    std::chrono::high_resolution_clock::time_point t1, t2;
};


#endif
