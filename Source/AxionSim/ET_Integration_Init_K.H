#ifndef ET_INTEGRATION_INIT_K_H
#define ET_INTEGRATION_INIT_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2]; 
        state_fab(i, j, k, Idx::Psi1) = std::exp(-std::pow(x0, 2) - std::pow(x1, 2) - std::pow(x2, 2));

        state_fab(i, j, k, Idx::Pi1) = 0;

        state_fab(i, j, k, Idx::Psi2) = -std::exp(-std::pow(x0, 2) - std::pow(x1, 2) - std::pow(x2, 2));

        state_fab(i, j, k, Idx::Pi2) = 0;

}
#endif