#ifndef AMREX_TIME_INTEGRATOR_H
#define AMREX_TIME_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class TimeIntegrator
{
private:
    int number_nodes;
    amrex::Vector<amrex::MultiFab> F_nodes;
    std::unique_pointer<amrex::MultiFab> S_tmp_ptr;
    amrex::MultiFab& S_old;
    amrex::MultiFab& S_new;
    amrex::MultiFab& S_tmp;
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;
    amrex::Vector<amrex::Real> weights;
    amrex::Vector<amrex::Real> nodes;

public:
    TimeIntegrator(MultiFab const& S_old_external, MultiFab const& S_new_external);
    Real advance(std::function<void(MultiFab&, const MultiFab&, const Real)> F, const Real time, const Real timestep);
};

#endif
