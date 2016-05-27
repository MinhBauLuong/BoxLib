
#ifndef BL_CUDARENA_H
#define BL_CUDARENA_H

#include <Arena.H>
#include <BL_CXX11.H>

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
    virtual void* alloc (std::size_t sz) BL_OVERRIDE;
    //
    // Deletes the arena pointed to by pt.
    //
    virtual void free (void* pt) BL_OVERRIDE;
};

#endif /*BL_CUDARENA_H*/
