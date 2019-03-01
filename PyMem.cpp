#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <vector>
#include <list>
#include "PyMem.h"

/******************************************************************************/
PyMem::PyMem(){
  length = 0;
}

/******************************************************************************/
void PyMem::put(std::pair<void*, int>pair)
{
    queue.push_back( pair );
    length += pair.second;
}

/******************************************************************************/
template<typename T>
void PyMem::put( std::vector<T> v )
{
  unsigned int len = v.size() * sizeof(v[0]);
  queue.push_back( { &v[0], len } );
  length += len;
}		 

/******************************************************************************/
void PyMem::free(void* data)
{
     for ( auto x = queue.begin(); x != queue.end(); x++ )
     {
	  if ( x->first == data )
	  {
	       queue.erase(x);
	       length -= x->second;
	       break;
	  }
     }
}

/******************************************************************************/
int PyMem::write( char* str )
{	
    int fd = shm_open( str, O_RDWR | O_CREAT, S_IRWXU );
    fcntl(fd, F_SETFD, 0 );
    void* buf = mmap( NULL, length, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0 );
    ftruncate(fd, length);
    void* addr = buf;
	
    for ( auto x : queue )
    {	    
	memcpy( addr, x.first, x.second );	    
	addr = (void*)((char*)addr + x.second);
    }
    return fd;
}
