#ifndef AMREX_TIME_INTEGRATOR_H
#define AMREX_TIME_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class TimeIntegrator
{
private:
    enum StateTimes {Old = 0, Tmp, New, NumStates};
    int number_nodes;
    amrex::Vector<amrex::MultiFab> F_nodes;
    amrex::Vector<amrex::MultiFab> S_val;
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;
    amrex::Vector<amrex::Real> weights;
    amrex::Vector<amrex::Real> nodes; 

public:
    TimeIntegrator(const BoxArray& ba, const DistributionMapping& dm, const int Ncomp, const int Nghost);
    amrex::MultiFab& get_old();
    amrex::MultiFab& get_new();
    amrex::MultiFab& get_tmp();
    Real advance(std::function<void(MultiFab&, const MultiFab&, const Real)> F, const Real time, const Real timestep);
};

#endif