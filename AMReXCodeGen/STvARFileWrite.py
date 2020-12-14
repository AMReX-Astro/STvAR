from STvAR import *

ijkString = ['i','j','k']

Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_rhs(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& rhs_fab, 
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Array4<amrex::Real> const& stress_energy_arr,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
""" 

Diag_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_diagnostics(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& diag, 
        amrex::Array4<amrex::Real const> const& state_fab,
        const amrex::Real time_lev,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
"""

LoadInit_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
initialize_from_data(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Array4<amrex::Real const> const& initial_data,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
"""

PostUpdate_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_post_update(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
"""

InitFromScratch_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
state_init(int i, int j, int k, 
        amrex::Array4<amrex::Real> const& state_fab, 
        amrex::Real time, const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo(); 
        
"""

ConvertedVariables_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
void 
fill_conv_vars_fab(int i, int j, int k, 
        amrex::Array4<const amrex::Real> const& state_fab,
        amrex::Array4<amrex::Real> const& rhs_fab,
        amrex::Array4<amrex::Real> const& conv_vars_arr,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx, 
        const amrex::GeometryData& geom) noexcept 
{
        const auto domain_xlo = geom.ProbLo();
        
"""

AMRtagging_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 
#include <ET_Integration_Setup_K.H> 

AMREX_GPU_DEVICE 
inline 
bool
state_is_tagged(int i, int j, int k, 
        amrex::Array4<amrex::Real const> const& state_fab,
        amrex::Real error_threshold,
        const amrex::Real time,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        const amrex::GeometryData& geom) noexcept 
{

        const auto domain_xlo = geom.ProbLo(); 

"""

Setup_Header_String = """

#include <AMReX_REAL.H> 
#include <AMReX_Array4.H> 

"""

def RHS_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_RHS_K_H" + "\n" + "#define ET_INTEGRATION_RHS_K_H" 
    Header = IFNDEF_DEFINE+Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def Diag_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_DIAG_K_H" + "\n" + "#define ET_INTEGRATION_DIAG_K_H" 
    Header = IFNDEF_DEFINE+Diag_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def LoadInit_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_LOADINITDATA_K_H" + "\n" + "#define ET_INTEGRATION_LOADINITDATA_K_H"
    Header = IFNDEF_DEFINE+LoadInit_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def PostUpdate_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_POST_UPDATE_K_H" + "\n" + "#define ET_INTEGRATION_POST_UPDATE_K_H"
    Header = IFNDEF_DEFINE+PostUpdate_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def InitFromScratch_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_INIT_K_H" + "\n" + "#define ET_INTEGRATION_INIT_K_H"
    Header = IFNDEF_DEFINE+InitFromScratch_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def ConvertedVariables_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_CONVERTED_VARS_K_H" + "\n" + "#define ET_INTEGRATION_CONVERTED_VARS_K_H"
    Header = IFNDEF_DEFINE+ConvertedVariables_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header
    
    
def AMRtagging_Header(version = "standard", dim = 3):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_AMR_CELL_TAG_K_H" + "\n" + "#define ET_INTEGRATION_AMR_CELL_TAG_K_H"
    Header = IFNDEF_DEFINE+AMRtagging_Header_String
    coord_string = ""
    for i in range(dim):
        coord_string += "        amrex::Real x" + str(i) + " = (" + ijkString[i] + " + 0.5)*geom.CellSize(" + str(i) + ") + domain_xlo[" + str(i) + "]; \n"
    Header += coord_string
    return Header

def Write_Setup_File(Statenames, Initnames, Diagnames, Convnames, Nghostcells = 2, version = "standard"):
    IFNDEF_DEFINE = "#ifndef ET_INTEGRATION_SETUP_K_H" + "\n" + "#define ET_INTEGRATION_SETUP_K_H"
    FileString = IFNDEF_DEFINE+Setup_Header_String
    
    IdxNamespaceStr = "namespace Idx { \n"
    IdxNamespaceStr += "         enum ETIndexes {"
    
    Idx_string = ""
    for itr in Statenames:
        Idx_string += itr+", "
    Idx_string += "NumScalars"
    
    IdxNamespaceStr += Idx_string
    IdxNamespaceStr += "}; \n};\n\n"
    
    FileString += IdxNamespaceStr
    
    InitIdxNamespaceStr = "namespace InitIdx { \n"
    InitIdxNamespaceStr += "         enum ETInitIndexes {"
    
    InitIdx_string = ""
    for itr in Initnames:
        InitIdx_string += itr+", "
    InitIdx_string += "NumScalars"
    
    InitIdxNamespaceStr += InitIdx_string
    InitIdxNamespaceStr += "}; \n};\n\n"
    
    FileString += InitIdxNamespaceStr
    
    DiagNamespaceStr = "namespace Diag { \n"
    DiagNamespaceStr += "         enum DiagnosticIndexes {"
    
    Diag_string = ""
    for itr in Diagnames:
        Diag_string += itr+", "
    Diag_string += "NumScalars"
    
    DiagNamespaceStr += Diag_string
    DiagNamespaceStr += "}; \n};\n\n"
    
    FileString += DiagNamespaceStr
    
    ConvNamespaceStr = "namespace ConvIdx { \n"
    ConvNamespaceStr += "         enum ConvertedIndexes {"
    
    Conv_string = ""
    for itr in Convnames:
        Conv_string += itr+", "
    Conv_string += "NumScalars"
    
    ConvNamespaceStr += Conv_string
    ConvNamespaceStr += "}; \n};\n\n"
    
    FileString += ConvNamespaceStr
    
    return FileString

def Write_Setup_File_Custom(Customnames, customIndex):
    CustomNamespaceStr = "namespace " + customIndex +"Idx { \n"
    CustomNamespaceStr += "         enum " + customIndex + "Indexes {"
    
    Custom_string = ""
    for itr in Customnames:
        Custom_string += itr+", "
    Custom_string += "NumScalars"
    
    CustomNamespaceStr += Custom_string
    CustomNamespaceStr += "}; \n};\n\n"
    
    FileString = CustomNamespaceStr
        
    return FileString

def Write_Setup_File_NGhost(Nghostcells):
    FileString = "#define NUM_GHOST_CELLS "+str(Nghostcells)+"\n\n"
    FileString += "#endif"
    return FileString

def Closer(version = "standard"):
    CloserString = "}\n"
    CloserString += "#endif"
    return CloserString

def VarString(Statenames, version = "standard"):
    varstring = "names = {"
    if len(Statenames) > 0:
        for itr in range(len(Statenames)-1):
            varstring += "\""+Statenames[itr]+"\", "
        varstring += "\""+Statenames[len(Statenames)-1]+"\""
    varstring += "};"
    return varstring





