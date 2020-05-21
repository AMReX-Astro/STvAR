#ifndef ET_INTEGRATION_INIT_K_H 
#define ET_INTEGRATION_INIT_K_H 

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup.H> 

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

        amrex::Real r = std::sqrt(std::pow(x0, 2) + std::pow(x1, 2) + std::pow(x2, 2));
        amrex::Real c1 = 1.168642873;
        amrex::Real c2 = -1.168642873;
        amrex::Real r1 = std::sqrt(std::pow(x1, 2) + std::pow(x2, 2) + std::pow(-c1 + x0, 2));
        amrex::Real r2 = std::sqrt(std::pow(x1, 2) + std::pow(x2, 2) + std::pow(-c2 + x0, 2));
        amrex::Real M1 = 0.453;
        amrex::Real M2 = 0.453;
        
        amrex::Real eps1;
        amrex::Real eps2;
        
        if(r1<= 0*M1)
          eps1 = 0.1*M1;
        else
          eps1 = 0;

        if(r2<= 0*M2)
          eps2 = 0.1*M2;
        else
          eps2 = 0;
        
        state_fab(i, j, k, Idx::chi) = std::pow((1.0/2.0)*M1/(r1+eps1) + (1.0/2.0)*M2/(r2+eps2) + 1,-4);

        state_fab(i, j, k, Idx::gamtildeLL00) = 1;

        state_fab(i, j, k, Idx::gamtildeLL01) = 0;

        state_fab(i, j, k, Idx::gamtildeLL02) = 0;

        state_fab(i, j, k, Idx::gamtildeLL11) = 1;

        state_fab(i, j, k, Idx::gamtildeLL12) = 0;

        state_fab(i, j, k, Idx::gamtildeLL22) = 1;

        state_fab(i, j, k, Idx::Khat) = 0;

        state_fab(i, j, k, Idx::AtildeLL00) = 0;

        state_fab(i, j, k, Idx::AtildeLL01) = 0;

        state_fab(i, j, k, Idx::AtildeLL02) = 0;

        state_fab(i, j, k, Idx::AtildeLL11) = 0;

        state_fab(i, j, k, Idx::AtildeLL12) = 0;

        state_fab(i, j, k, Idx::AtildeLL22) = 0;

        state_fab(i, j, k, Idx::theta) = 0;

        state_fab(i, j, k, Idx::GamtildeU0) = 0;

        state_fab(i, j, k, Idx::GamtildeU1) = 0;

        state_fab(i, j, k, Idx::GamtildeU2) = 0;

        state_fab(i, j, k, Idx::alpha) = std::pow((1.0/2.0)*M1/r1 + (1.0/2.0)*M2/r2 + 1, -2);

        state_fab(i, j, k, Idx::betaU0) = 0;

        state_fab(i, j, k, Idx::betaU1) = 0;

        state_fab(i, j, k, Idx::betaU2) = 0;

}
#endif