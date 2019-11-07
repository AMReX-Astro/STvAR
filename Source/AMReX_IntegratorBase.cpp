#include "AMReX_IntegratorBase.H"

using namespace amrex;

IntegratorBase::IntegratorBase(amrex::MultiFab& S_old_external,
                               amrex::MultiFab& S_new_external,
                               amrex::Real initial_time) : time(initial_time),
                                                           S_old(S_old_external),
                                                           S_new(S_new_external) {}

void IntegratorBase::rhs(amrex::MultiFab& S_rhs, const amrex::MultiFab& S_data, const amrex::Real time)
{
    Fun(S_rhs, S_data, time);
}

amrex::MultiFab& IntegratorBase::get_new_data()
{
    return S_new;
}

amrex::MultiFab& IntegratorBase::get_old_data()
{
    return S_old;
}

amrex::Real IntegratorBase::get_time()
{
    return time;
}
