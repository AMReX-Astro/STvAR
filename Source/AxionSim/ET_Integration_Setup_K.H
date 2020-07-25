#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {Psi1, Pi1, Psi2, Pi2, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {ModPsi, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif