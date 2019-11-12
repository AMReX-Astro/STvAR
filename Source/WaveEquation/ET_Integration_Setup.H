#ifndef ET_INTEGRATION_SETUP_H
#define ET_INTEGRATION_SETUP_H

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>

namespace Idx {
    enum ETIndexes {Phi = 0, Pi, NumScalars};
};

#define NUM_GHOST_CELLS 1

#endif