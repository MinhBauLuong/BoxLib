#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <BLFort.H>

// This tutorial shows how to translate a Fortran kernel into a C++ kernel. It
// currently works only for regular Boxes, not tiled Boxes (that comes next).

// Prototype for the C++ kernel
void work_cpp(const Box* bx,
              const int ng,
              Real *dataptr);

// Prototype for the Fortran kernel
extern "C"
{
  void work(const int* lo, const int* hi, BL_FORT_FAB_ARG(dfab));
}

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    const int domain_size = 256;
    const int max_box_size = 64;

    // define the lower and upper corner of a 3D domain
    IntVect domain_lo(0, 0, 0);
    IntVect domain_hi(domain_size-1,domain_size-1,domain_size-1);

    // build a box for the domain
    Box domain(domain_lo, domain_hi);

    // build a box array from the 64^3 domain box
    BoxArray ba(domain);
    // break the box array into 32^3 boxes
    ba.maxSize(max_box_size);

    // build a multifab on the box array with 1 component, 2 ghost cells
    MultiFab data(ba, 1, 2);

    // First do the Fortran kernel.
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Calling Fortran kernel ... ";
    }
    for (MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        work(bx.loVect(), bx.hiVect(),
        BL_TO_FORTRAN(data[mfi]));
    }
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "done." << std::endl;
    }

    Real min = data.min(0);
    Real max = data.max(0);
    Real norm0 = data.norm0();
    Real norm1 = data.norm1();
    Real norm2 = data.norm2();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Results from Fortran kernel:" << std::endl;
        std::cout << "min      = " << min  << std::endl;
        std::cout << "max      = " << max  << std::endl;
        std::cout << "max norm = " << norm0 << std::endl;
        std::cout << "L1  norm = " << norm1 << std::endl;
        std::cout << "L2  norm = " << norm2 << std::endl;
    }

    // Now do the C++ kernel. It should produce identical results as the
    // Fortran kernel.
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Calling C++ kernel ... ";
    }
    for (MFIter mfi(data); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Real *dataptr = data[mfi].dataPtr();
        const int ng = data.nGrow();
        work_cpp(&bx, ng, dataptr);
    }
    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "done." << std::endl;
    }

    min = data.min(0);
    max = data.max(0);
    norm0 = data.norm0();
    norm1 = data.norm1();
    norm2 = data.norm2();

    if (ParallelDescriptor::IOProcessor()) {
        std::cout << "Results from C++ kernel:" << std::endl;
        std::cout << "min      = " << min  << std::endl;
        std::cout << "max      = " << max  << std::endl;
        std::cout << "max norm = " << norm0 << std::endl;
        std::cout << "L1  norm = " << norm1 << std::endl;
        std::cout << "L2  norm = " << norm2 << std::endl;
    }


    // print some information for checking
    /* I got
       min      = 7.945268926e-11
       max      = 1
       max norm = 1
       L1  norm = 8680.319857
       L2  norm = 56.24354515
    */

    BoxLib::Finalize();
}

