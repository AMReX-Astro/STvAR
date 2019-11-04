#ifndef ET_INTEGRATION_K_H
#define ET_INTEGRATION_K_H
#include "ET_Integration_Indexes.H"

AMREX_GPU_DEVICE
inline
void
state_rhs(int i, int j, int k,
          amrex::Array4<amrex::Real> const& rhs_fab, 
          amrex::Array4<amrex::Real const> const& state_fab,
          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept
{
    rhs_fab(i, j, k, IPHI) = state_fab(i, j, k, IPI);
    rhs_fab(i, j, k, IPI)  = (state_fab(i, j + 1, k, IPHI) + state_fab(i, j - 1, k, IPHI) - 2*state_fab(i, j, k, IPHI))/(dx[1]*dx[1]) + (state_fab(i + 1, j, k, IPHI) + state_fab(i - 1, j, k, IPHI) - 2*state_fab(i, j, k, IPHI))/(dx[0]*dx[0]);  
}
#endif
