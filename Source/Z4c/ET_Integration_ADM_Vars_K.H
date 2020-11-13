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

        adm_vars_array(i, j, k, ADMIdx::gamLL00) = std::pow(state_fab(i, j, k, Idx::chi), -1)*state_fab(i, j, k, Idx::gamtildeLL00);
        adm_vars_array(i, j, k, ADMIdx::gamLL01) = std::pow(state_fab(i, j, k, Idx::chi), -1)*state_fab(i, j, k, Idx::gamtildeLL01);
        adm_vars_array(i, j, k, ADMIdx::gamLL02) = std::pow(state_fab(i, j, k, Idx::chi), -1)*state_fab(i, j, k, Idx::gamtildeLL02);
        adm_vars_array(i, j, k, ADMIdx::gamLL11) = std::pow(state_fab(i, j, k, Idx::chi),  -1)*state_fab(i, j, k, Idx::gamtildeLL11);
        adm_vars_array(i, j, k, ADMIdx::gamLL12) = std::pow(state_fab(i, j, k, Idx::chi), -1)*state_fab(i, j, k, Idx::gamtildeLL12);
        adm_vars_array(i, j, k, ADMIdx::gamLL22) = std::pow(state_fab(i, j, k, Idx::chi), -1)*state_fab(i, j, k, Idx::gamtildeLL22);
        
        adm_vars_array(i, j, k, ADMIdx::alpha) = state_fab(i, j, k, Idx::alpha);
        adm_vars_array(i, j, k, ADMIdx::betaU0) = state_fab(i, j, k, Idx::betaU0);
        adm_vars_array(i, j, k, ADMIdx::betaU1) = state_fab(i, j, k, Idx::betaU1);
        adm_vars_array(i, j, k, ADMIdx::betaU2) = state_fab(i, j, k, Idx::betaU2);

}
#endif