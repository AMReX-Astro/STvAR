#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {u, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace ADMIdx {
         enum ETAdmIndexes {gamLL00, gamLL01, gamLL02, gamLL11, gamLL12, gamLL22, alpha, betaU0, betaU1, betaU2, NumScalars};
};

namespace Diag { 
         enum DiagnosticIndexes {chi, alpha, NumScalars}; 
};

#define NUM_GHOST_CELLS 4

#endif