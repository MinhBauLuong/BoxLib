// --------------------------------------------------------------------------
// GridMoveTest.cpp
// --------------------------------------------------------------------------
//  this file tests the performance for moving grids
// --------------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <VisMF.H>

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

using std::cout;
using std::endl;

const int maxGrid(32);
const int nComp(8);
const int nGhost(1);
const int XDIR(0);
const int nTimes(4);


// --------------------------------------------------------------------------
void UniqueSet(Array<int> &uSet, int setSize, int poolSize) {   // ---- a unique set of random numbers
  if(setSize > poolSize) {
    BoxLib::Abort("**** Error in UniqueSet:  setSize > poolSize.");
  }
  std::set<int> copySet;
  while(copySet.size() < setSize) {
    int r(BoxLib::Random_int(poolSize));
    if(copySet.find(r) == copySet.end()) {
      copySet.insert(r);
      uSet.push_back(r);
    }
  }
  for(int i(0); i < uSet.size(); ++i) {
    std::cout << "uSet[" << i << "]  = " << uSet[i] << std::endl;
  }
}


// --------------------------------------------------------------------------
void SetFabValsToPMap(MultiFab &mf) {
  const Array<int> &newDMA = mf.DistributionMap().ProcessorMap();
  for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
    const int index(mfi.index());
    FArrayBox &fab = mf[mfi];
    for(int i(0); i < fab.nComp(); ++i) {  // ---- setVal to distmap index
      Real val(newDMA[index] + ((static_cast<Real> (i)) / fab.nComp()));
      fab.setVal(val, i);
    }
  }
}


// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);    

    BL_PROFILE_VAR("main()", pmain);
    BL_PROFILE_REGION_START("main");

    int nProcs(ParallelDescriptor::NProcs());
    int myProc(ParallelDescriptor::MyProc());
    bool copyAll(false), copyRandom(false), moveFabs(true);

    // ---- make a box, then a boxarray with maxSize
    Box baseBox(IntVect(0,0,0), IntVect(maxGrid - 1, maxGrid - 1, maxGrid - 1));
    Array<Box> boxes(nProcs);

    for(int p(0); p < boxes.size(); ++p) {
      boxes[p] = baseBox.shift(XDIR, maxGrid);
    }
    BoxArray ba(boxes.dataPtr(), boxes.size());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "ba.size() = " << ba.size() << std::endl;
      std::cout << "ba[0] = " << ba[0] << std::endl;
      std::cout << "ba[1] = " << ba[1] << std::endl;
      std::cout << "ba[last] = " << ba[ba.size() - 1] << std::endl;
    }

    // ---- make a multifab, setval to the index
    MultiFab mf(ba, nComp, nGhost);
    for(int i(0); i < mf.nComp(); ++i) {
      mf.setVal(static_cast<Real> (i), i, 1);
    }
    if(ParallelDescriptor::IOProcessor()) {
      cout << "******** mf.DistributionMap().linkCount = "
           << mf.DistributionMap().linkCount() << endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "dmap = " << mf.DistributionMap() << std::endl;
      std::cout << "dmap.size() = " << mf.DistributionMap().size() << std::endl;
    }

    // ------------------------------------------------------------------
    // ----- very important:  here we are copying a procmap,
    // -----                  but if you just make your own Array<int>
    // -----                  it must have an extra value at the end
    // -----                  set to ParallelDescriptor::MyProc()
    // -----                  see DistributionMapping.H
    // ------------------------------------------------------------------
    const Array<int> &procMap = mf.DistributionMap().ProcessorMap();
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "procMap.size() = " << procMap.size() << std::endl;
    }


    // ------------------ copy entire multifab to nProcs/2) % nProcs
    if(copyAll) {
      // ---- initialize it to (oldmap + (nProcs/2)) % nProcs
      Array<int> newMap(procMap.size());
      for(int i(0); i < procMap.size(); ++i) {
        newMap[i] = (procMap[i] + (nProcs/2)) % nProcs;
      }
      DistributionMapping newDMap(newMap);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "newDMap = " << newDMap << std::endl;
      }

      // ---- make a new multifab with the new map and copy from mf
      MultiFab mfNewMap;
      mfNewMap.define(ba, nComp, nGhost, newDMap, Fab_allocate);
      mfNewMap.setVal(-42.0);

      // ---- now copy from mf
      BL_PROFILE_REGION_START("MFCopy_ProcPlusNPd2");
      BL_PROFILE_VAR("MFCopy_ProcPlusNPd2", mfcppnpd2);
      for(int n(0); n < nTimes; ++n) {
        mfNewMap.copy(mf);
      }
      BL_PROFILE_VAR_STOP(mfcppnpd2);
      BL_PROFILE_REGION_STOP("MFCopy_ProcPlusNPd2");
    }


    // ------------------ copy a random set of grids elsewhere
    if(copyRandom) {
      int nCopies(nProcs/10);
      nCopies = std::min(nCopies, nProcs/2);
      int nRanksInSet(nCopies * 2);

      Array<int> copyArray;
      if(ParallelDescriptor::IOProcessor()) {
        UniqueSet(copyArray, nRanksInSet, nProcs);
      } else {
        copyArray.resize(nRanksInSet);
      }

      ParallelDescriptor::Barrier();
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Starting random copy." << std::endl;
      }
      ParallelDescriptor::Bcast(copyArray.dataPtr(), copyArray.size());
        
      Array<int> newMap(nCopies + 1);
      BoxArray baCopy(nCopies);
      for(int i(0); i < nCopies; ++i) {
        newMap[i] = copyArray[i];
	baCopy.set(i, ba[copyArray[i + nCopies]]);
      }
      newMap[nCopies] = ParallelDescriptor::MyProc();

      DistributionMapping newDMap(newMap);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "newDMap = " << newDMap << std::endl;
      }

      // ---- make a new multifab with the new map and copy from mf
      MultiFab mfNewMap;
      mfNewMap.define(baCopy, nComp, nGhost, newDMap, Fab_allocate);
      mfNewMap.setVal(-42.0);

      // ---- now copy from mf
      BL_PROFILE_REGION_START("MFCopy_RandomSet");
      BL_PROFILE_VAR("MFCopy_RandomSet", mfcprs);
      for(int n(0); n < nTimes; ++n) {
        mfNewMap.copy(mf);
      }
      BL_PROFILE_VAR_STOP(mfcprs);
      BL_PROFILE_REGION_STOP("MFCopy_RandomSet");
    }


    if(moveFabs) {
      VisMF::SetNOutFiles(1);
      Array<int> copyArray;
      int nRanksInSet(8);
      if(nRanksInSet % 2 != 0) {
        BoxLib::Abort("**** Bad nRanksInSet");
      }

      if(ParallelDescriptor::IOProcessor()) {
        UniqueSet(copyArray, nRanksInSet, nProcs);
      } else {
        copyArray.resize(nRanksInSet);
      }
      ParallelDescriptor::Bcast(copyArray.dataPtr(), copyArray.size());

      Array<int> newDistMapArray(mf.DistributionMap().ProcessorMap());
      for(int n(0); n < nRanksInSet; n += 2) {
        newDistMapArray[copyArray[n]] = copyArray[n+1];
      }

      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "newDistMapArray = " << newDistMapArray << std::endl;
      }
      DistributionMapping::CacheStats(std::cout);

      SetFabValsToPMap(mf);
      VisMF::Write(mf, "mfOriginal");

      mf.MoveFabs(newDistMapArray);

      SetFabValsToPMap(mf);
      VisMF::Write(mf, "mfMoves");

      BoxArray ba16(mf.boxArray());
      ba16.maxSize(16);

      Array<int> randomMap(ba16.size());
      if(ParallelDescriptor::IOProcessor()) {
	for(int ir(0); ir < ba16.size(); ++ ir) {
          randomMap[ir] = BoxLib::Random_int(nProcs);
	}
      }
      ParallelDescriptor::Bcast(randomMap.dataPtr(), randomMap.size());
        
      Array<int> newMap16(ba16.size() + 1);
      for(int i(0); i < newMap16.size() - 1; ++i) {
        newMap16[i] = randomMap[i];
      }
      newMap16[newMap16.size() - 1] = myProc;
      DistributionMapping newDMap16(newMap16);
      MultiFab mf16;
      mf16.define(ba16, nComp, nGhost, newDMap16, Fab_allocate);
      mf16.setVal(-1.0);
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "mf16DMA = " << mf16.DistributionMap() << std::endl;
      }
      SetFabValsToPMap(mf16);
      VisMF::Write(mf16, "mf16DMap");
      mf16.copy(mf);
      VisMF::Write(mf16, "mf16");
    }


    BL_PROFILE_REGION_STOP("main");
    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
