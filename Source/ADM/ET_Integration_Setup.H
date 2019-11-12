#ifndef ET_INTEGRATION_SETUP_H
#define ET_INTEGRATION_SETUP_H

namespace Idx {
    enum ETIndexes {alpha, betax, betay, betaz, 
                    gxx, gxy, gxz, gyx, gyy, gyz, gzx, gzy, gzz,
                    ginvxx, ginvxy, ginvxz, ginvyx, ginvyy, ginvyz, ginvzx, ginvzy, ginvzz,
                    Kxx, Kxy, Kxz, Kyx, Kyy, Kyz, Kzx, Kzy, Kzz, NumScalars};
};

#define NUM_GHOST_CELLS 1

#endif