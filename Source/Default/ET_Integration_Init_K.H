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
    state_fab(i, j, k, Idx::Phi) = 0.0;
}
#endif
