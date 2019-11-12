#ifndef ET_INTEGRATION_RHS_K_H
#define ET_INTEGRATION_RHS_K_H

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <ET_Integration_Setup.H>

AMREX_GPU_DEVICE
inline
void
state_rhs(int i, int j, int k,
          amrex::Array4<amrex::Real> const& rhs_fab, 
          amrex::Array4<amrex::Real const> const& state_fab,
          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept
{
    rhs_fab(i, j, k, Idx::Phi) = 1.0;
}
#endif
