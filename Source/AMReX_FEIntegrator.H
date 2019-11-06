#ifndef AMREX_FE_INTEGRATOR_H
#define AMREX_FE_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class FEIntegrator
{
private:
    std::unique_ptr<amrex::MultiFab> F_tmp_ptr;
    amrex::MultiFab& S_old;
    amrex::MultiFab& S_new;

public:
    FEIntegrator(amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external);
    amrex::Real advance(std::function<void(amrex::MultiFab&, amrex::MultiFab&, amrex::Real)> f, const amrex::Real time, const amrex::Real timestep);
};

#endif
