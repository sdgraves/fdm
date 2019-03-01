#include <iostream>
#include <vector>
#include <limits>
#include "ode.h"
#include <cstdlib>

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
    {
	series *temp = a;
	a = b;
	b = temp;
    }
    
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
    double d = 1/(T);
    
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
    else if ( N == 1 )	
    {
	v = {d, 0, -d};
    }
    else
    {
	v = {d};
    }
    return v;
}



/******************************************************************************/
std::vector<double> diff_backwards( std::vector<double> u, double T )
{
    std::vector<double> v (u.size(), 0);
    int N = u.size();

    if ( N > 1 )
    {

	//i.e. { 1, 2 }
	//to   { 1, -1 + 2, -2 }
	v[N-1] = u[N-1]/T;
	v.push_back( -u[N-1]/T );
	
	for ( int ii = 0; ii < N - 1; ii++ )
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
//Adding a sequence {0, [0], 0} to {[1]} might produce {0, [1], 0}; the
//z-transform class, however, reinterprets this sequence as causal,
//i.e. {[0], 1, 0}. Hence the check in this method for v[ii] == 0, We still have
//{ [a], b } +  { [-a], c } = { [0] b+c }.
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
	    s1 = s1 + temp;
	}
    }
    return s1;    
}



/******************************************************************************/
ZT<double> from_spoly( std::vector<double> num, std::vector<double> den, double T, int descending=1 )
{

    series q = series::sample( T, den, descending );
    series p = series::sample( T, num, descending );

    //The ZT class always assumes the 0th register corresponds to the current input/output.
    if ( descending )
    {
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
    }

    return ZT<double> { p.taps, q.taps, T, descending };

}

//gets next time step -- geometric series which limits to eps.
double g( double T, double rate )
{
    double eps = std::numeric_limits<double>::epsilon();
    double x = rate*T + eps;
    return ( x == T ) ? eps : x;
}



//proposed fix to current method -- here the function f is required to keep track for itself of
//needed states, time variance.
template<typename Token, typename Fct>     
std::vector<Token> odeint_noinit( ZT<Token> A, int steps, Fct f, int feedback=0 )
{

    std::vector<Token> vec ( steps );
    Token x,y;
     	  
    for ( int i = 0; i < steps; i++ ){
	x = f(y);
	y = A(feedback*y + x);
	vec[i] = y;
	A.status();
	A.update();
    }
    return vec;
}


template<typename Token, typename Fct>
std::vector<Token> euler( Fct f, double t0, double t1, double T, ZT<Token>p, int feedback=0 )
{
    Token x;
    Token y = p[0];
    std::vector<Token> ret( 1 + (t1 - t0)/T );

    int ii = 0;
    for ( double t = t0; t <= t1; t += T )
    {
	x = f(y,t);
	y = p(feedback*y + x);
	ret[ii++] = y;
	p.update();
//	p.status();
    }
    return ret;
}


//This is a method for setting initial conditions on sequences.
//starting with the case of least generality, assume the difference used is a backwards difference.
template<typename T>
void h( T d, ZT<T> *sys )
{
    //the solution satisfies
    // d = (sys[0] - sys[1])/T
    //or
    // Td - sys[0] = -sys[1]
    // sys[0] - Td = sys[1]
    T val = (*sys)[0] - (sys->rate() * d) ;
    sys->set_out( 1, val );
    
}


//draft of adaptive method
//Note that if the coefficients on the differential equation are allowed to be functions, we have the means
//to handle time-varying and linearized systems. These will probably require additional overhead.
template<typename Token, typename Fct>
void f( Fct forcing, rational r, double t0, double t1, Token threshold, ivp init,
	int feedback=0)
{

    constexpr double min_step = 1E-7;
    static_assert( min_step > 2*std::numeric_limits<double>::epsilon(),
		   "min sample period may reach series limit\n" );    

    //default steps is 100
    double T = (t1 - t0)/100;
    ZT<Token> p = from_spoly( r.num, r.den, T, 0 );
    h(1., &p);
    std::vector<Token> r1 = euler(forcing, t0, t1, T, p, 0);
    double err;
    do {
	std::cout << "sample period: " << T << "\n";
        T = g(0.75, T);

	ZT<Token> q =  from_spoly( r.num, r.den, T, 0 );
	h(1., &q);

	std::vector<Token> r2 = euler(forcing, t0, t1, T, q, 0);
	//the error should be the maximum error, but the systems have to be up/downsampled first.
	err = r2[ r2.size() - 1 ] - r1[ r1.size() - 1 ];
	if ( err < 0 )
	    err = -err;
	std::swap(r1, r2);
    } while ( ( T > min_step ) && ( err > threshold ) );
    
    if ( err > threshold )
	std::cout << "max resolution reached before desired accuracy\n";
    
}


/******************************************************************************/
int main( int argc, char** argv )
{
//    rational s { {1}, {1, 0, 1} };
    rational s { {1}, {1, 1} };

    ivp init { 0., 10. };
    auto x = [](double d, double t) -> double{
	return 0;
    };
    f( x, s, 0, 1, 0.00001, init);
    return 0;
}
