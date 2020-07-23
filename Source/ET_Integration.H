#ifndef ET_INTEGRATION_H
#define ET_INTEGRATION_H
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <ET_Integration_Setup_K.H>
#include <ET_Integration_Init_K.H>
#include <ET_Integration_Rhs_K.H>
#include <ET_Integration_Diag_K.H>
#include <ET_Integration_Post_Update_K.H>
#include <ET_Integration_AMR_Cell_Tag_K.H>

#ifdef PROBLEM_LOADS_INITDATA
#include <ET_Integration_LoadInitData_K.H>
#endif

namespace Variable
{
    extern amrex::Vector<std::string> names;
    void Initialize();
}

namespace Diagnostics
{
    extern amrex::Vector<std::string> names;
    void Initialize();
}
#endif
