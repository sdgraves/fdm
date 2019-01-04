#ifndef PYMEM
#define PYMEM
#include <list>
#include <vector>
/******************************************************************************/
/* This is a structure that holds pointers and their corresponding lengths
 * at which data is to be placed in shared memory for use by matplotlib. 
/******************************************************************************/
class PyMem
{
public:
    PyMem();
    void put( std::pair<void*, int>pair );
    template <typename T>
	void put( std::vector<T> v );
    void free(void* data);
    int write( char* str );
    inline unsigned int size(){ return length; }
private:
    std::list<std::pair<void*, unsigned int>> queue;
    unsigned int length;
};


#endif
