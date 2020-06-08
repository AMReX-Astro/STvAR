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
        amrex::Real gamtildeLL00 = state_fab(i, j, k, Idx::gamtildeLL00);
        amrex::Real gamtildeLL01 = state_fab(i, j, k, Idx::gamtildeLL01);
        amrex::Real gamtildeLL02 = state_fab(i, j, k, Idx::gamtildeLL02);
        amrex::Real gamtildeLL10 = state_fab(i, j, k, Idx::gamtildeLL01);
        amrex::Real gamtildeLL11 = state_fab(i, j, k, Idx::gamtildeLL11);
        amrex::Real gamtildeLL12 = state_fab(i, j, k, Idx::gamtildeLL12);
        amrex::Real gamtildeLL20 = state_fab(i, j, k, Idx::gamtildeLL02);
        amrex::Real gamtildeLL21 = state_fab(i, j, k, Idx::gamtildeLL12);
        amrex::Real gamtildeLL22 = state_fab(i, j, k, Idx::gamtildeLL22);
        amrex::Real detgamtilde = gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20;
        amrex::Real scalefactor = std::pow(detgamtilde, -0.33333333333333331);
        state_fab(i, j, k, Idx::gamtildeLL00) = gamtildeLL00*scalefactor;

        state_fab(i, j, k, Idx::gamtildeLL01) = gamtildeLL01*scalefactor;

        state_fab(i, j, k, Idx::gamtildeLL02) = gamtildeLL02*scalefactor;

        state_fab(i, j, k, Idx::gamtildeLL11) = gamtildeLL11*scalefactor;

        state_fab(i, j, k, Idx::gamtildeLL12) = gamtildeLL12*scalefactor;

        state_fab(i, j, k, Idx::gamtildeLL22) = gamtildeLL22*scalefactor;

        amrex::Real AtildeLL00 = state_fab(i, j, k, Idx::AtildeLL00);
        amrex::Real AtildeLL01 = state_fab(i, j, k, Idx::AtildeLL01);
        amrex::Real AtildeLL02 = state_fab(i, j, k, Idx::AtildeLL02);
        amrex::Real AtildeLL10 = state_fab(i, j, k, Idx::AtildeLL01);
        amrex::Real AtildeLL11 = state_fab(i, j, k, Idx::AtildeLL11);
        amrex::Real AtildeLL12 = state_fab(i, j, k, Idx::AtildeLL12);
        amrex::Real AtildeLL20 = state_fab(i, j, k, Idx::AtildeLL02);
        amrex::Real AtildeLL21 = state_fab(i, j, k, Idx::AtildeLL12);
        amrex::Real AtildeLL22 = state_fab(i, j, k, Idx::AtildeLL22);
        amrex::Real gamtildeUU00 = (gamtildeLL11*gamtildeLL22 - gamtildeLL12*gamtildeLL21)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU01 = (-gamtildeLL01*gamtildeLL22 + gamtildeLL02*gamtildeLL21)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU02 = (gamtildeLL01*gamtildeLL12 - gamtildeLL02*gamtildeLL11)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU10 = (-gamtildeLL10*gamtildeLL22 + gamtildeLL12*gamtildeLL20)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU11 = gamtildeLL00*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU12 = -gamtildeLL00*(gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU20 = (gamtildeLL10*gamtildeLL21 - gamtildeLL11*gamtildeLL20)/(gamtildeLL00*gamtildeLL11*gamtildeLL22 - gamtildeLL00*gamtildeLL12*gamtildeLL21 - gamtildeLL01*gamtildeLL10*gamtildeLL22 + gamtildeLL01*gamtildeLL12*gamtildeLL20 + gamtildeLL02*gamtildeLL10*gamtildeLL21 - gamtildeLL02*gamtildeLL11*gamtildeLL20);
        amrex::Real gamtildeUU21 = -gamtildeLL00*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real gamtildeUU22 = gamtildeLL00*(gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)/((gamtildeLL00*gamtildeLL11 - gamtildeLL01*gamtildeLL10)*(gamtildeLL00*gamtildeLL22 - gamtildeLL02*gamtildeLL20) - (gamtildeLL00*gamtildeLL12 - gamtildeLL02*gamtildeLL10)*(gamtildeLL00*gamtildeLL21 - gamtildeLL01*gamtildeLL20));
        amrex::Real TrAtilde = 0 + AtildeLL00*gamtildeUU00 + AtildeLL01*gamtildeUU01 + AtildeLL02*gamtildeUU02 + AtildeLL10*gamtildeUU10 + AtildeLL11*gamtildeUU11 + AtildeLL12*gamtildeUU12 + AtildeLL20*gamtildeUU20 + AtildeLL21*gamtildeUU21 + AtildeLL22*gamtildeUU22;
        state_fab(i, j, k, Idx::AtildeLL00) = AtildeLL00 - 0.33333333333333331*TrAtilde*gamtildeLL00;

        state_fab(i, j, k, Idx::AtildeLL01) = AtildeLL01 - 0.33333333333333331*TrAtilde*gamtildeLL01;

        state_fab(i, j, k, Idx::AtildeLL02) = AtildeLL02 - 0.33333333333333331*TrAtilde*gamtildeLL02;

        state_fab(i, j, k, Idx::AtildeLL11) = AtildeLL11 - 0.33333333333333331*TrAtilde*gamtildeLL11;

        state_fab(i, j, k, Idx::AtildeLL12) = AtildeLL12 - 0.33333333333333331*TrAtilde*gamtildeLL12;

        state_fab(i, j, k, Idx::AtildeLL22) = AtildeLL22 - 0.33333333333333331*TrAtilde*gamtildeLL22;

}
#endif