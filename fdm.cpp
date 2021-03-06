#include <iostream>
#include <cstdio>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cassert>
#include <errno.h>
#include <string.h>
#include <cstdint>
#include <complex>
#include <algorithm>
#include <map>
#include "ode.h"
#include "PyMem.h"
#include "timer.h"

/******************************************************************************/
//The level curve class just facilitates sending a vector of doubles to the
//shared memory region for use by matplotlib.

typedef struct
{
    double step;
    double start;
    unsigned long len;
}  dataHeader;

class levelCurve{
private:
    double T; //sample period
    double t0; //initial time
    std::vector<double>state;
public:
    levelCurve(double sample_T=0.01, double start=0.) : T(sample_T), t0(start)
    {	
	state = std::vector<double>();
	state.push_back(0.);
    }
    
    levelCurve(std::vector<double> v, double sample_T=0.01, double start=0.)
	: state(std::move(v)), T(sample_T), t0(start) { }
  
    inline dataHeader makeHeader()
    {
	return dataHeader{ T, t0, state.size() };
    }
   
    inline std::pair<void*, int> makeKey()
    {
	return { &state[0], state.size() * sizeof(state[0]) };
    }
    inline void push( double d )
    {
	state.push_back(d);
    }
    static dataHeader end_header;
};

dataHeader levelCurve::end_header = dataHeader {0., 0., 0};


/******************************************************************************/
/*
struct commonTransforms{
    ZT<double> secondOrder{ {1}, { 1/(T*T), w0*w0 - (2/T*T) + (2*sigma/T), (1/T*T) - (2*sigma/T), 0. } };
    
    ZT<double> hold{ {0}, {1., 0.,} };
    ZT<double> integrator{ {T}, {1., 0.,} };     
    ZT<double> differentiator{ {-1., 1.}, {0,} };
    ZT<double> geometric{ {0}, { 0.8, 0.} };
    ZT<double> exponential{ {T}, { ( 1 - T), 0 } };
    ZT<double> exp_middlediff{};
    ZT<double> average{ {0.333, 0.333, 0.333}, {0.} };
    
    ZT<double> trig{ { sin(T), 0}, { 1., -2*sin(T), 1.0, 0. } };
    ZT<double> doubleIntegral{ {T*T}, {-1., 0.} };
    };*/


/******************************************************************************/
int main( int argc, char **argv )
{

    //parameters for simulation
    double time_max = 10;
    unsigned int steps = 10;     
    int x = 0;
    int v = 1;
    Timer t{};

    

    if ( argc == 3 ){
	time_max = atof( argv[1] );
	steps = atoi( argv[2] );
    }
    double T = ( time_max / steps );

    double zeta = 0.25;
    double w0 = 10;
    double sigma = zeta*w0;

    std::cout << "Sample T is " << T << "\n";
    std::vector<double> den = series::sample( T, {1., 0.} );
    
    ZT<double> test{ {0}, den };
    test.print();

    std::vector<double> in( steps, 0 );
    t.start();

    std::vector<double> results =  odeint<double>( test, in, 1., 0 );
    levelCurve c1 { results };

    dataHeader h_end = {0., 0., 0};


    PyMem p = PyMem();
    dataHeader h1 = c1.makeHeader();
    dataHeader h2 = levelCurve::end_header;
    p.put( {&h1, sizeof(h1)} );
    p.put( c1.makeKey() );
    p.put( {&h2, sizeof(h2)} );
    
    int fd = p.write("shm_out_file");

    t.stop();    

    //This isn't my favorite, but execve takes c strings.
    char c[2];
    char d[9];
    sprintf(c, "%d", fd);
    sprintf(d, "%u", p.size() );
          
    const char* args[] = {"python3", "plotting.py", c, d, (char*)NULL};
    execve ("/usr/bin/python3", const_cast<char**>(args), environ);
    
    printf("exiting with %d: %s\n", errno, strerror(errno));     
    return -errno;
    
    
}

/******************************************************************************/
