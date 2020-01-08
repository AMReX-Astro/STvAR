#ifndef ET_INTEGRATION_SETUP_K_H 
#define ET_INTEGRATION_SETUP_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {phi, gbar00, gbar01, gbar02, gbar11, gbar12, gbar22, Ksclr, Abar00, Abar01, Abar02, Abar11, Abar12, Abar22, GambarU0, GambarU1, GambarU2, alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif