#include <iostream>
#include <vector>
#include <algorithm>
#include "ode.h"

void series::print(){
    for ( auto x : taps )
	std::cout << x << '\t';
    std::cout << '\n';
}
    
series series::operator+( series f )
{

    int len;
    series *a = &f;
    series *b = this;
    
    if ( f.taps.size() < taps.size() )
	std::swap(a, b);
    
    len = a->taps.size();
    std::vector<double> res { a->taps };
    
    for ( int ii = 0; ii < b->taps.size(); ii++ )
	res[ii + a->t0 - b->t0] += b->taps[ii];
    
    return series { res, a->t0 };
    
}

series series::operator*( double d )
{
    int len = this->taps.size();
    std::vector<double> v( len );
    
    for ( int ii = 0; ii < len; ii++ )
	v[ii] = this->taps[ii] * d;

    return { v, this->t0 };
}

series series::diff( int order, double T )
{
  
    std::vector<double> u {1.};
    int ii;
    for ( ii = 0; ii < order; ii++ )
    {
	u = diff_backwards(u, T);
    }
    series s { u, ii };
    return s;
}


/******************************************************************************/
std::vector<double> diff_center( std::vector<double> u, double T )
{
    std::vector<double> v (u.size());
    int N = u.size();
    double d = 1/(2*T);
    
    if ( N > 1 )
    {
	//edge cases
	v[1] -= u[0]*d;
	v.insert( v.begin(), d );
	v[N-1] += u[N-1]*d;
	v.push_back( -u[N-1]*d );
	//center difference
	
	for ( int ii = 1; ii < u.size() - 2; ii++ )
	{
	    v[ii-1] += u[ii]*d;
	    v[ii+1] -= u[ii]*d;
	}	
    }    
    else	
    {
	v = {d, 0, -d};
    }
    return v;
}



/******************************************************************************/
std::vector<double> diff_backwards( std::vector<double> u, double T )
{
    std::vector<double> v (u.size());
    int N = u.size();

    if ( N > 1 )
    {

	//i.e. { 1, 2 }
	//to   { 1, -1 + 2, -2 }
	v[N-1] += u[N-1]/T;
	v.push_back( -u[N-1]/T );
	
	for ( int ii = 0; ii < N - 2; ii++ )
	{
	    v[ii]   += u[ii]/T;
	    v[ii+1] -= u[ii]/T;
	}

    }
    else
    {
	v = {1/T, -1/T};
    }
    return v;

}

/******************************************************************************/
//Adding a sequence {0, [0], 0} to {[1]} sensibly produces {0, [1], 0}; the
//z-transform class, however, reinterprets this sequence as causal,
//i.e. {[0], 1, 0}. Hence the check in this method for v[ii] == 0.
/******************************************************************************/
series series::sample( double T, std::vector<double> v, int descending )
{
    if ( descending )
	std::reverse( v.begin(), v.end() );
    
    series s1 { std::vector<double>{v[0]}, 0 };
    for ( int ii = 1; ii < v.size(); ii++ )
    {
	if ( v[ii] != 0 ){
	    series temp = diff(ii, T) * v[ii];
	    std::cout << "with ii = " << ii << ", appending by\n";
	    temp.print();
	    s1 = s1 + temp;
	}
    }
    return s1;    
}



/******************************************************************************/
ZT<double> from_spoly( std::vector<double> num, std::vector<double> den, double T, int descending=0 )
{
    series q = series::sample( T, den, descending );
    series p = series::sample( T, num, descending );

    for ( int ii = 0; ii < q.taps.size(); ii++ )
	q.taps[ii] *= -1;

    if ( q.taps.size() > p.taps.size() )
    {
	series s { std::vector<double>( q.taps.size(), 0 ), q.offset() };
	p = p + s;
    }
    else
    {
	series s { std::vector<double>( p.taps.size(), 0 ), p.offset() };
	q = q + s;
    }

    
    return ZT<double> { p.taps, q.taps };

}


/******************************************************************************/
int main( int argc, char** argv )
{
    ZT<double> test = from_spoly( {1.}, {1., 1.4141, 1.}, 0.01, 0 );
    test.print();
    std::vector<double> forcing( 100, 0. );
    std::vector<double> d = odeint( test, forcing, 1. );

    return 0;
}
