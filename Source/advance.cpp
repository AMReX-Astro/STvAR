#include "ET_Integration.H"
#include "AMReX_FEIntegrator.H"

using namespace amrex;

void advance(MultiFab& state_new_mf, MultiFab& state_old_mf, Real time, Real dt, const Geometry& geom)
{
    int ncomp = state_new_mf.nComp();

    // Make a time integrator
    FEIntegrator integrator(state_old_mf, state_new_mf);

    // Create a RHS source function
    auto source_fun = [&](MultiFab& rhs, MultiFab& state, const Real time){
      fill_state_rhs(rhs, state, geom);
    };

    // Call the time integrator advance
    integrator.advance(source_fun, time, dt);

    // Fill ghost cells for each grid from valid regions of another grid
    state_new_mf.FillBoundary(geom.periodicity());
}
