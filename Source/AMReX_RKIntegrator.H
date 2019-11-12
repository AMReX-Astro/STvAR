#ifndef AMREX_RK_INTEGRATOR_H
#define AMREX_RK_INTEGRATOR_H
#include <functional>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include "AMReX_IntegratorBase.H"

namespace ButcherTableauTypes {
    enum ButcherTypes {User = 0, ForwardEuler, Trapezoid, SSPRK3, RK4, NumTypes};
};

class RKIntegrator : public IntegratorBase
{
private:
    int tableau_type;
    int number_nodes;
    bool use_adaptive_timestep;
    amrex::Vector<amrex::MultiFab> F_nodes;
    std::unique_ptr<amrex::MultiFab> S_tmp_ptr;
    amrex::Vector<amrex::Vector<amrex::Real> > tableau;
    amrex::Vector<amrex::Real> weights;
    amrex::Vector<amrex::Real> extended_weights;
    amrex::Vector<amrex::Real> nodes;
    void initialize_preset_tableau();
    void initialize_parameters();
    void initialize_stages();

public:
    RKIntegrator(amrex::MultiFab& S_old_external, amrex::MultiFab& S_new_external, amrex::Real time);
    amrex::Real advance(const amrex::Real timestep);
};

#endif
