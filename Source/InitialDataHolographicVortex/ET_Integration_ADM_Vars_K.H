#ifndef ET_INTEGRATION_ADM_VARS_K_H
#define ET_INTEGRATION_ADM_VARS_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
fill_adm_vars_fab(int i, int j, int k, 
        amrex::Array4<const amrex::Real> const& state_fab,
        amrex::Array4<amrex::Real> const& adm_vars_array) noexcept 
{
        //ADM variables are written from Z4c variables. 

}
#endif