#ifndef ET_INTEGRATION_H
#define ET_INTEGRATION_H
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

void main_main();
void init_phi(amrex::MultiFab& phi_new, amrex::Real time, const amrex::Geometry& geom);
void advance_phi(amrex::MultiFab& phi_new, amrex::MultiFab& phi_old, amrex::Real time, amrex::Real dt, const amrex::Geometry& geom);
void fill_state_rhs(amrex::MultiFab& phi_rhs_mf, amrex::MultiFab& phi_old_mf, const amrex::Geometry& geom);
#endif
