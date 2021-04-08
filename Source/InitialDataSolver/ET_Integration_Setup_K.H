#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {u, psi, chi, alphafunc, betafunc, Atilde_LL_00, Atilde_LL_01, Atilde_LL_02, Atilde_LL_11, Atilde_LL_12, Atilde_LL_22, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {Hpsi, Hchi, Hu, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif