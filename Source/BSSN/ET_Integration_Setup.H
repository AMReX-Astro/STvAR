#ifndef ET_INTEGRATION_SETUP_K_H 
#define ET_INTEGRATION_SETUP_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {phi, gambar00, gambar01, gambar02, gambar11, gambar12, gambar22, Ksclr, Abar00, Abar01, Abar02, Abar11, Abar12, Abar22, GambarU0, GambarU1, GambarU2, alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2, NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {H, NumScalars}; 
};

#define NUM_GHOST_CELLS 2

#endif