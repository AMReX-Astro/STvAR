#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {chi, gamtildeLL00, gamtildeLL01, gamtildeLL02, gamtildeLL11, gamtildeLL12, gamtildeLL22, Khat, AtildeLL00, AtildeLL01, AtildeLL02, AtildeLL11, AtildeLL12, AtildeLL22, theta, GamtildeU0, GamtildeU1, GamtildeU2, alpha, betaU0, betaU1, betaU2, rho, epsE, vel30, vel31, vel32, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {u, NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {H, ML0, ML1, ML2, C, Ksclr, NumScalars}; 
};

namespace ConvIdx { 
         enum ConvertedIndexes {Source10, Source11, Source12, Source2, NumScalars}; 
};

namespace SEIdx { 
         enum SEIndexes {rho_fluid, S_fluidL0, S_fluidL1, S_fluidL2, S_fluidLL00, S_fluidLL01, S_fluidLL02, S_fluidLL11, S_fluidLL12, S_fluidLL22, NumScalars}; 
};

#define NUM_GHOST_CELLS 4

#endif