#ifndef ET_INTEGRATION_K_H
#define ET_INTEGRATION_K_H

AMREX_GPU_DEVICE
inline
void
phi_rhs(int i, int j, int k, int n, 
        amrex::Array4<amrex::Real> const& phi_rhs_fab, 
        amrex::Array4<amrex::Real const> const& phi_old_fab,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept
{
    // Set phi_rhs_fab = 1.0
    phi_rhs_fab(i, j, k, n) = 1.0;
}
#endif