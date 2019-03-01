#ifndef ODEINT
#define ODEINT

#include <vector>
#include <complex>
#include <iostream>
#include <algorithm>


std::vector<double> diff_center( std::vector<double> u, double T=1. );
std::vector<double> diff_backwards( std::vector<double> u, double T=1. );

/******************************************************************************/
//The center difference for recursive use to find difference equations from
//arbitrary constant-coefficient linear differential equations, and a helper
//class to allow adding of vectors corresponding to filter impulse responses
//where t=0 corresponds to different samples in the addend and the augend.

/******************************************************************************/
//assume default time scale is one sample; then this is the backwards difference by default.

class series
{
    int t0;    
public:
    std::vector<double> taps;
    series( std::vector<double> v, int t ) : taps {v}, t0(t) { }
    void print();
    series operator+( series f );
    series operator*( double d );
    static series sample( double T, std::vector<double> v, int descending=1 );
    static series diff( int order, double T=1. );
    inline int offset(){ return t0; }
};
    
/******************************************************************************/
typedef struct rational{
  std::vector<double> num;
  std::vector<double> den;
} rational;

typedef struct ivp{
    double x;
    double v;
} ivp;




/******************************************************************************/
//This is a SISO system class designed for difference equations. The operator()
//returns the value of the current state. Class members only save their own
//time history as far as is needed for the difference equations which govern them
//and handle input history similarly.

//The primary motivation for using a template is to save memory when using
//floats over complex doubles, but should work fine for MIMO systems as long as
//the template T as long as {T} is a field and { T x C } with C the complex numbers
//is a vector space.

//The subscript operator is not guaranteed to fetch the time history of the system
//for subscripts not used by the difference equation, and returns by value.
/******************************************************************************/
template<typename T>
class ZT
{
    using vec = std::vector<T>;
public:
 ZT( vec num=vec{1}, vec den=vec{1}, double t=1, int descending=1 ) : p(den), q(num),
      sample_T(t)
    {

	regs_out = vec(p.size(), 0);
	regs_in = vec(q.size(), 0);

	
	if ( descending )
	{
	    std::reverse(p.begin(), p.end());
	    std::reverse(q.begin(), q.end());
	}

	//set the coefficient on y[k] to 1 to eliminate one multiply
	//this has unintended consequences when either input or output are zero
	if ( p[0] != 1 && p[0] != 0 )
	{
	   T temp = p[0];
	   for ( int ii = 0; ii < p.size(); ii++ )
		p[ii] = p[ii]/temp;
	
	  for ( int ii = 0; ii < q.size(); ii++ )
		q[ii] = q[ii]/temp;
	}
	
    }	
    ZT( ZT<T> f, ZT<T> g )
    {
	p = expand( f.p, g.p );
	q = expand( f.q, g.q );
	regs_out = vec(p.size());
	regs_in = vec(q.size());
    }
    
    //assume difference equations are given in the form
    //
    // a_n x[k-n] + a_n-1 k[k-n+1] + ... + a_0 x[k] = 0;
    //
    //then
    //
    //   s1 = Sum{ a_n x[k-n] + a_n-1 x[k-n+1] + ... + a_1 x [k-1];
    // x[k] = -s1 (for a_0 = 1 ).
    inline T operator()( T x )
    {

	regs_in[0] = x;

	T s1 {};
	T s2 {};
	
	for ( int i = 1;  i < regs_out.size(); i++ )
	    s1 += ( regs_out[i] * p[i] );

	for ( int i = 0; i < regs_in.size(); i++ )
	    s2 += ( regs_in[i] * q[i] );

	regs_out[0] = ( s2 - s1 );
	return ( s2 - s1 );	  
    }

    void status()
    {
	std::cout << "operating with y = (";
	for ( auto ii : regs_out )
	    std::cout << ii << ", ";
	std::cout << ") and x = (";
	for ( auto ii : regs_in )
	    std::cout << ii << ", ";
	std::cout << ")\n";
    }
    inline void update()
    {
	for ( int i = regs_out.size() - 1; i > 0; i-- )
	    regs_out[i] = regs_out[i-1];

	for ( int i = regs_in.size() - 1; i > 0; i-- )
	    regs_in[i] = regs_in[i-1];
	  
    }
    void IVP( T x )
    {
	for ( int i = 0; i < regs_out.size(); i++ )
	{
	    regs_out[i] = x;
	    std::cout << "regs_out " << i << ": " << regs_out[i] << "\n";
	}
    }
    inline T operator[]( int i )
    {
	return ( (0 - i) < regs_out.size() ) ? regs_out[0-i] : T{0} ;
    }
    void set_in( int i, T value )
    {
	regs_in[i] = value;
    }
    void set_out( int i, T value )
    {
	regs_out[i] = value;
    }

    void print()
    {
	std::cout << "(";
	for ( auto x : q )
	    std::cout << ' ' << x << ',';
	std::cout << ") / (";
	for ( auto x : p )
	    std::cout << ' ' << x << ',';
	std::cout << ")\n";
    }
    inline double rate()
    {
	return sample_T;
    }
private:
    vec p;
    vec q;
    vec regs_out; //may want to use std::rotate to maintain regs
    vec regs_in;
    double sample_T;
    vec expand( vec v, vec u )
    {
	
	vec res( v.size() + u.size() - 1 );
	for (int i = (u.size() - 1); i >= 0; i-- )
	{
	    for ( int j = 0; j < v.size(); j++ )
		res[j + i] += v[j]*u[i];	     
	}
	return res;
    }
};

template<typename Token>     
std::vector<Token> odeint( ZT<Token> A, std::vector<Token>x, Token init, int feedback=0 )
{

    std::vector<Token> vec ( x.size() );
    Token y;
    A.IVP( init );
     	  
    for ( int i = 0; i < x.size(); i++ ){
	y = A(feedback*y + x[i]);
	vec[i] = y;
	A.status();
	A.update();
    }
    return vec;
}

#endif
