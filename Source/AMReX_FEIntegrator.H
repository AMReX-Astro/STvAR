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
    std::function<void(amrex::MultiFab&, amrex::MultiFab&, amrex::Real)> Fun;
    amrex::Real time;

public:
    FEIntegrator(std::function<void(amrex::MultiFab&, amrex::MultiFab&, amrex::Real)> F, amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external, amrex::Real time);
    amrex::Real advance(const amrex::Real timestep);
    amrex::MultiFab& get_new_data();
    amrex::MultiFab& get_old_data();
    amrex::Real get_time();
};

#endif
