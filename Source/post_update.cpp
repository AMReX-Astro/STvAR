#include "ET_Integration.H"

using namespace amrex;

void rescale_state (MultiFab& state_mf)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
  for ( MFIter mfi(state_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.tilebox();
    const auto ncomp = state_mf.nComp();

    const auto& state_fab = state_mf.array(mfi);

    // For each grid, loop over all the valid points
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
	// for example:
	amrex::Real det = state_fab(i, j, k, Idx::gbar00)*state_fab(i, j, k, Idx::gbar11)*state_fab(i, j, k, Idx::gbar22) - state_fab(i, j, k, Idx::gbar00)*std::pow(state_fab(i, j, k, Idx::gbar12), 2) - std::pow(state_fab(i, j, k, Idx::gbar01), 2)*state_fab(i, j, k, Idx::gbar22) + 2*state_fab(i, j, k, Idx::gbar01)*state_fab(i, j, k, Idx::gbar02)*state_fab(i, j, k, Idx::gbar12) - std::pow(state_fab(i, j, k, Idx::gbar02), 2)*state_fab(i, j, k, Idx::gbar11);
	state_fab(i, j, k, Idx::gbar00) = state_fab(i, j, k, Idx::gbar00) * std::pow(det,-1/3);
	state_fab(i, j, k, Idx::gbar01) = state_fab(i, j, k, Idx::gbar01) * std::pow(det,-1/3);
	state_fab(i, j, k, Idx::gbar02) = state_fab(i, j, k, Idx::gbar02) * std::pow(det,-1/3);
	state_fab(i, j, k, Idx::gbar11) = state_fab(i, j, k, Idx::gbar11) * std::pow(det,-1/3);
	state_fab(i, j, k, Idx::gbar12) = state_fab(i, j, k, Idx::gbar12) * std::pow(det,-1/3);
	state_fab(i, j, k, Idx::gbar22) = state_fab(i, j, k, Idx::gbar22) * std::pow(det,-1/3);
    });
  }
}
