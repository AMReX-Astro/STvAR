#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H
#define ET_INTEGRATION_AMR_CELL_TAG_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k, 
        int lev,
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        //amrex::Real alpha = state_fab(i, j, k, Idx::alpha);
        
        //amrex::Real v = 2.0;
        
        amrex::Real PhiReal = state_fab(i, j, k, Idx::PhiReal);
        amrex::Real PhiImaginary = state_fab(i, j, k, Idx::PhiImaginary);
        
        amrex::Real ModPhi = std::sqrt(std::pow(PhiReal,2)+std::pow(PhiImaginary,2));
        
        return ModPhi < error_threshold; 

        //return false; //std::fabs(alpha) < error_threshold;
}
#endif