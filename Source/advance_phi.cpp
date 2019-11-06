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

    /*
    // Create a MultiFab containing the time integration RHS
    MultiFab rhs_mf(phi_new_mf.boxArray(), phi_new_mf.DistributionMap(), ncomp, 0);
    fill_phi_rhs(rhs_mf, phi_old_mf, geom);

    // Loop over grids to do a forward euler integration in time
    for ( MFIter mfi(phi_new_mf); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();

      const auto& phi_new_fab = phi_new_mf.array(mfi);
      const auto& phi_old_fab = phi_old_mf.array(mfi);
      const auto&     rhs_fab =     rhs_mf.array(mfi);

      // For each grid, loop over all the valid points
      AMREX_FOR_4D(bx, ncomp, i, j, k, n,
      {
         // Right now rhs_fab(i,j,k,n) = 1 so this adds dt to every value in every time step
         phi_new_fab(i,j,k,n) = phi_old_fab(i,j,k,n) + dt * rhs_fab(i,j,k,n);
      });
    }
    */

    // Fill ghost cells for each grid from valid regions of another grid
    phi_new_mf.FillBoundary();
}
