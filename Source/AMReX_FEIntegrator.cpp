#include "AMReX_FEIntegrator.H"

using namespace amrex;

FEIntegrator::FEIntegrator(std::function<void(amrex::MultiFab&, amrex::MultiFab&, amrex::Real)> F,
                           amrex::MultiFab& S_old_external, 
                           amrex::MultiFab& S_new_external, 
                           amrex::Real initial_time) : Fun(F),
                                                       time(initial_time),
                                                       S_old(S_old_external),
                                                       S_new(S_new_external)
{
    // Create temporary MultiFab
    F_tmp_ptr = std::make_unique<MultiFab>(S_old.boxArray(), S_old.DistributionMap(), S_old.nComp(), S_old.nGrow());
}

Real FEIntegrator::advance(const Real timestep)
{
    // Assume before advance() that S_new is valid data at the current time ("time" argument)
    // So we update S_old by copying the current state.
    MultiFab::Copy(S_old, S_new, 0, 0, S_old.nComp(), S_old.nGrow());

    // F_tmp = RHS(S_old, t_old)
    MultiFab& F_tmp = *F_tmp_ptr;
    Fun(F_tmp, S_old, time);

    // S_new += timestep * dS/dt
    MultiFab::Saxpy(S_new, timestep, F_tmp, 0, 0, S_new.nComp(), S_new.nGrow());

    // Update time
    time += timestep;

    // Return timestep
    return timestep;
}

amrex::MultiFab& FEIntegrator::get_new_data()
{
    return S_new;
}

amrex::MultiFab& FEIntegrator::get_old_data()
{
    return S_old;
}

amrex::Real FEIntegrator::get_time()
{
    return time;
}