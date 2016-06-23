
#ifndef BL_CUDARENA_H
#define BL_CUDARENA_H

#ifdef BL_CUDA

#include <Arena.H>

//
// A Concrete Class for Dynamic Memory Management
//
// This class derived from Arena handles memory allocations using CUDA's
// Unified Virtual Memory.
//

class CUDArena
    :
    public Arena
{
public:
    //
    // Allocates a dynamic memory arena of size sz.
    // Returns a pointer to this memory.
    //
    virtual void* alloc (std::size_t sz);
    //
    // Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt);
};

#endif /* BL_CUDA */

#endif /*BL_CUDARENA_H*/
