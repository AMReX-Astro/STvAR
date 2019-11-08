#ifndef AMREX_TIME_INTEGRATOR_H
#define AMREX_TIME_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "AMReX_IntegratorBase.H"
#include "AMReX_FEIntegrator.H"
#include "AMReX_RKIntegrator.H"

namespace IntegratorTypes {
    enum TimeIntegratorTypes {ForwardEuler = 0,
                              ExplicitRungeKutta};
};

class TimeIntegrator
{
private:
    amrex::Real time;
    amrex::Real timestep;
    int step_number;
    std::unique_ptr<IntegratorBase> integrator_ptr;
    std::function<void ()> post_timestep;

public:
    TimeIntegrator(amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external, amrex::Real time);
    void set_post_timestep(std::function<void ()> F);
    void set_post_update(std::function<void (amrex::MultiFab&)> F);
    void set_rhs(std::function<void(amrex::MultiFab&, const amrex::MultiFab&, const amrex::Real)> F);
    void integrate(const amrex::Real start_timestep, const amrex::Real end_time, const int nsteps);
    amrex::MultiFab& get_new_data();
    amrex::MultiFab& get_old_data();
    amrex::Real get_time();
    amrex::Real get_timestep();
    int get_step_number();
};

#endif
