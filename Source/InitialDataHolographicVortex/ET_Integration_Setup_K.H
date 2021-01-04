#ifndef ET_INTEGRATION_SETUP_K_H
#define ET_INTEGRATION_SETUP_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

namespace Idx { 
         enum ETIndexes {Phi, AL3, NumScalars}; 
};

namespace InitIdx { 
         enum ETInitIndexes {Q6, Q7, NumScalars}; 
};

namespace Diag { 
         enum DiagnosticIndexes {NumScalars}; 
};

namespace ConvIdx { 
         enum ConvertedIndexes {NumScalars}; 
};

namespace SEIdx { 
         enum SEIndexes {NumScalars}; 
};

#define NUM_GHOST_CELLS 4

#endif