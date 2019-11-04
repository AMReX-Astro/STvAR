#ifndef ET_INTEGRATION_K_H
#define ET_INTEGRATION_K_H

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>

namespace Idx {
    enum ETIndexes {Phi = 0, Pi};
};

AMREX_GPU_DEVICE
inline
void
state_rhs(int i, int j, int k,
          amrex::Array4<amrex::Real> const& rhs_fab, 
          amrex::Array4<amrex::Real const> const& state_fab,
          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept
{
    rhs_fab(i, j, k, Idx::Phi) = state_fab(i, j, k, Idx::Pi);
    rhs_fab(i, j, k, Idx::Pi)  = (state_fab(i, j + 1, k, Idx::Phi) + state_fab(i, j - 1, k, Idx::Phi) - 2*state_fab(i, j, k, Idx::Phi))/(dx[1]*dx[1]) + (state_fab(i + 1, j, k, Idx::Phi) + state_fab(i - 1, j, k, Idx::Phi) - 2*state_fab(i, j, k, Idx::Phi))/(dx[0]*dx[0]);  
}

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
