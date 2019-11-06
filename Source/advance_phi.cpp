#include "ET_Integration.H"
#include "AMReX_FEIntegrator.H"

using namespace amrex;

void advance_phi (MultiFab& phi_new_mf, MultiFab& phi_old_mf, Real time, Real dt, const Geometry& geom)
{
    int ncomp = phi_new_mf.nComp();

    // Fill ghost cells for each grid from valid regions of another grid
    phi_old_mf.FillBoundary();

    // Make a time integrator
    FEIntegrator integrator(phi_old_mf, phi_new_mf);

    // Create a RHS source function
    auto source_fun = [&](MultiFab& rhs, MultiFab& state, const Real time){
      fill_phi_rhs(rhs, state, geom);
    };

    // Call the time integrator advance
    integrator.advance(source_fun, time, dt);

    // Fill ghost cells for each grid from valid regions of another grid
    phi_new_mf.FillBoundary();
}
