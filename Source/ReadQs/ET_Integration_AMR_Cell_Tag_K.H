#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H
#define ET_INTEGRATION_AMR_CELL_TAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k, 
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{

        return state_fab(i, j, k, Idx::Q6) < error_threshold;
}
#endif