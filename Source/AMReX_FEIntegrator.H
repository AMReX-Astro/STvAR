#ifndef AMREX_FE_INTEGRATOR_H
#define AMREX_FE_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include "AMReX_IntegratorBase.H"

class FEIntegrator : public IntegratorBase
{
private:
    std::unique_ptr<amrex::MultiFab> F_tmp_ptr;

public:
    FEIntegrator(amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external, amrex::Real time);
    amrex::Real advance(const amrex::Real timestep);
};

#endif
