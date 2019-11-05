#ifndef ET_INTEGRATION_K_H
#define ET_INTEGRATION_K_H

#if defined(WAVEEQ)
#include "ET_Integration_K_WaveEq.H"
#elif defined(MAXWELL)
#include "ET_Integration_K_Maxwell.H"
#else
#include "ET_Integration_K_Default.H"
#endif

#endif
