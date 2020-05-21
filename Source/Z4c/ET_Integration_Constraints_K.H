#ifndef ET_INTEGRATION_CONSTRAINTS_K_H 
#define ET_INTEGRATION_CONSTRAINTS_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab, 
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0];
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1];
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];

       
        diag(i, j, k, Diag::H) = 0;

}
#endif