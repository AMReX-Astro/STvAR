#ifndef ET_INTEGRATION_K_DEFAULT_H
#define ET_INTEGRATION_K_DEFAULT_H

#include <AMReX_REAL.H>
#include <AMReX_Array4.H>

namespace Idx {
    enum ETIndexes {Phi = 0, NumScalars};
};

#define NUM_GHOST_CELLS 0

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
