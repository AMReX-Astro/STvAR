#ifndef ET_INTEGRATION_POST_UPDATE_K_H
#define ET_INTEGRATION_POST_UPDATE_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
}
#endif