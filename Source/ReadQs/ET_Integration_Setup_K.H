#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {PhiRe, PiPhiRe, PhiIm, PiPhiIm, AT, PiAT, A0, A1, PiA0, PiA1, Psi, alpha, gambar00, gambar01, gambar02, gambar11, gambar12, gambar22, NumScalars};
};

namespace InitIdx { 
         enum ETInitIndexes {NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {NumScalars}; 
};

#define NUM_GHOST_CELLS 4

#endif