#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <ET_Integration_Setup.H>

AMREX_GPU_DEVICE
inline
void
state_init(int i, int j, int k,
           amrex::Array4<amrex::Real> const& state_fab,
           amrex::Real time,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept
{
    // init phi
    amrex::Real x = (i+0.5) * dx[0];
    state_fab(i, j, k, Idx::Phi) = sin(2.0*M_PI*(x - time));
    // init pi
    state_fab(i, j, k, Idx::Pi) = -2.0*M_PI*cos(2.0*M_PI*(x-time));
}
#endif
