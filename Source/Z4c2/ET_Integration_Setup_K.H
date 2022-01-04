#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {chi, gamtilde_LL_00, gamtilde_LL_01, gamtilde_LL_02, gamtilde_LL_11, gamtilde_LL_12, gamtilde_LL_22, Khat, Atilde_LL_00, Atilde_LL_01, Atilde_LL_02, Atilde_LL_11, Atilde_LL_12, Atilde_LL_22, theta, Gamtilde_U_0, Gamtilde_U_1, Gamtilde_U_2, alpha, beta_U_0, beta_U_1, beta_U_2, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {u, psi, chi, alphafunc, betafunc, Atilde_LL_00, Atilde_LL_01, Atilde_LL_02, Atilde_LL_11, Atilde_LL_12, Atilde_LL_22, NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {gamtilde_LL_00sqrd, gam_LL_00sqrd, Rsclr, sqrtKsqrd, K, H, M, Z, Csqrd, C, NumScalars}; 
         //enum DiagnosticIndexes {Psi4re, Psi4im, ModPsi4, AbsAlpha, Ksclr, NumScalars};
};

#define NUM_GHOST_CELLS 4

#endif