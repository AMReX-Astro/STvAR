

#ifndef ET_INTEGRATION_STRESS_ENERGY_K_H
#define ET_INTEGRATION_STRESS_ENERGY_K_H

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
fill_stress_energy_fab(int i, int j, int k,
        amrex::Array4<amrex::Real> const& stress_energy_arr,
        amrex::Array4<const amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept
{
        const auto domain_xlo = geom.ProbLo();
        
        amrex::Real x0 = (i + 0.5)*geom.CellSize(0) + domain_xlo[0]; 
        amrex::Real x1 = (j + 0.5)*geom.CellSize(1) + domain_xlo[1]; 
        amrex::Real x2 = (k + 0.5)*geom.CellSize(2) + domain_xlo[2];


        amrex::Real chi = state_fab(i, j, k, Idx::chi);

        amrex::Real gamtildeLL00 = state_fab(i, j, k, Idx::gamtildeLL00);
        amrex::Real gamtildeLL01 = state_fab(i, j, k, Idx::gamtildeLL01);
        amrex::Real gamtildeLL02 = state_fab(i, j, k, Idx::gamtildeLL02);
        amrex::Real gamtildeLL10 = state_fab(i, j, k, Idx::gamtildeLL01);
        amrex::Real gamtildeLL11 = state_fab(i, j, k, Idx::gamtildeLL11);
        amrex::Real gamtildeLL12 = state_fab(i, j, k, Idx::gamtildeLL12);
        amrex::Real gamtildeLL20 = state_fab(i, j, k, Idx::gamtildeLL02);
        amrex::Real gamtildeLL21 = state_fab(i, j, k, Idx::gamtildeLL12);
        amrex::Real gamtildeLL22 = state_fab(i, j, k, Idx::gamtildeLL22);
        amrex::Real alpha = state_fab(i, j, k, Idx::alpha);

        amrex::Real betaU0 = state_fab(i, j, k, Idx::betaU0);
        amrex::Real betaU1 = state_fab(i, j, k, Idx::betaU1);
        amrex::Real betaU2 = state_fab(i, j, k, Idx::betaU2);
        amrex::Real gamtildeUU00 = (gamtildeLL11*gamtildeLL22 - gamtildeLL12*gamtildeLL21)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU01 = (-gamtildeLL01*gamtildeLL22 + gamtildeLL02*gamtildeLL21)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU02 = (gamtildeLL01*gamtildeLL12 - gamtildeLL02*gamtildeLL11)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU10 = (-gamtildeLL10*gamtildeLL22 + gamtildeLL12*gamtildeLL20)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU11 = gamtildeLL00*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU12 = -gamtildeLL00*(gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU20 = (gamtildeLL10*gamtildeLL21 - gamtildeLL11*gamtildeLL20)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU21 = -gamtildeLL00*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU22 = gamtildeLL00*(gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamLL00 = gamtildeLL00/chi;
        amrex::Real gamLL01 = gamtildeLL01/chi;
        amrex::Real gamLL02 = gamtildeLL02/chi;
        amrex::Real gamLL10 = gamtildeLL10/chi;
        amrex::Real gamLL11 = gamtildeLL11/chi;
        amrex::Real gamLL12 = gamtildeLL12/chi;
        amrex::Real gamLL20 = gamtildeLL20/chi;
        amrex::Real gamLL21 = gamtildeLL21/chi;
        amrex::Real gamLL22 = gamtildeLL22/chi;
        amrex::Real gamUU00 = chi*gamtildeUU00;
        amrex::Real gamUU01 = chi*gamtildeUU01;
        amrex::Real gamUU02 = chi*gamtildeUU02;
        amrex::Real gamUU10 = chi*gamtildeUU10;
        amrex::Real gamUU11 = chi*gamtildeUU11;
        amrex::Real gamUU12 = chi*gamtildeUU12;
        amrex::Real gamUU20 = chi*gamtildeUU20;
        amrex::Real gamUU21 = chi*gamtildeUU21;
        amrex::Real gamUU22 = chi*gamtildeUU22;
        amrex::Real P_fluid = 0;

        amrex::Real iota = 0.10000000000000001;

        amrex::Real h_fluid = 0;

        amrex::Real W_fluid = 1;

        amrex::Real u3L0 = 0;
        amrex::Real u3L1 = 0;
        amrex::Real u3L2 = 0;
        stress_energy_arr(i, j, k, SEIdx::rho_fluid) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidL0) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidL1) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidL2) = 0;

        amrex::Real S_fluidL0 = 0;
        amrex::Real S_fluidL1 = 0;
        amrex::Real S_fluidL2 = 0;
        stress_energy_arr(i, j, k, SEIdx::S_fluidLL00) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidLL01) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidLL02) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidLL11) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidLL12) = 0;

        stress_energy_arr(i, j, k, SEIdx::S_fluidLL22) = 0;

}
#endif