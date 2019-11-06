#include "AMReX_FEIntegrator.H"

using namespace amrex;

FEIntegrator::FEIntegrator(MultiFab& S_old_external, MultiFab& S_new_external) : S_old(S_old_external), S_new(S_new_external)
{
    // Create temporary MultiFab and store a reference to it
    F_tmp_ptr = std::make_unique<MultiFab>(S_old.boxArray(), S_old.DistributionMap(), S_old.nComp(), S_old.nGrow());
}

Real FEIntegrator::advance(std::function<void(MultiFab&, MultiFab&, Real)> f, const Real time, const Real timestep)
{
    // Assume before advance() that S_new is valid data at the current time ("time" argument)
    // So we update S_old by copying the current state.
    MultiFab::Copy(S_old, S_new, 0, 0, S_old.nComp(), S_old.nGrow());

    // F_tmp = RHS(S_old, t_old)
    MultiFab& F_tmp = *F_tmp_ptr;
    f(F_tmp, S_old, time);

    MultiFab::Saxpy(S_new, timestep, F_tmp, 0, 0, S_new.nComp(), S_new.nGrow());

    // Return timestep
    return timestep;
}