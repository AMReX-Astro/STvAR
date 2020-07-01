#ifndef ET_INTEGRATION_POST_UPDATE_K_H 
#define ET_INTEGRATION_POST_UPDATE_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        const amrex::GeometryData& geom) noexcept 
{
}
#endif