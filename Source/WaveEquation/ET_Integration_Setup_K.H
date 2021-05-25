#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {PhiReal, PhiImaginary, PiPhiReal, PiPhiImaginary, Chi, Chi0, Chi1, PiChi, PiChi0, PiChi1, A_L_0, A_L_1, A_L_2, AT, PiA_L_0, PiA_L_1, PiA_L_2, PiAT, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {ModPhi, divA, ModChi, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif