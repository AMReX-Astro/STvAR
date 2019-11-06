#ifndef AMREX_INTEGRATOR_BASE_H
#define AMREX_INTEGRATOR_BASE_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class IntegratorBase
{
private:
    std::function<void(amrex::MultiFab&, const amrex::MultiFab&, const amrex::Real)> Fun;

protected:
    amrex::MultiFab& S_old;
    amrex::MultiFab& S_new;
    amrex::Real time;

public:
    IntegratorBase(std::function<void(amrex::MultiFab&, const amrex::MultiFab&, const amrex::Real)> F, amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external, amrex::Real time);
    virtual amrex::Real advance(const amrex::Real timestep) { return timestep; }
    void rhs(amrex::MultiFab& S_rhs, const amrex::MultiFab& S_data, const amrex::Real time);
    amrex::MultiFab& get_new_data();
    amrex::MultiFab& get_old_data();
    amrex::Real get_time();
};

#endif
