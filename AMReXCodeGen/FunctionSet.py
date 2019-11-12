import NRPy_param_funcs as par
# The indexedexp module defines various functions for defining and managing indexed quantities like tensors and pseudotensors
import indexedexp as ixp
# The grid module defines various parameters related to a numerical grid or the dimensionality of indexed expressions
# For example, it declares the parameter DIM, which specifies the dimensionality of the indexed expression
import grid as gri
from outputC import *
import sympy
from sympy import symbols, IndexedBase, Indexed, Idx, preorder_traversal
import numpy as np


Nx, Ny, Nz, Nn= symbols('Nx Ny Nz Nn', integer=True)
i = Idx('i', Nx)
j = Idx('j', Ny)
k = Idx('k', Nz)
n = Idx('n', Nn)
dx, dy, dz = symbols('dx dy dz')

directions = ['x','y','z']

def shift(E, idx_shift):
    # This function takes a generic Sympy expression and
    # returns a new Sympy expression where every Sympy Indexed
    # object in E has been shifted by idx_shift.
    # - idx_shift should be of length D, the dimension of E
    
    def shift_indexed(S, idx_shift):
        # This function returns a new IndexedBase object with shifted indices
        # - S should be a Sympy Indexed object
        # - idx_shift should be a tuple or list of index offsets to apply
        # - idx_shift should be of length D, the dimension of S
        base = S.base
        indices = [si + di for si, di in zip(S.indices, idx_shift)]
        return base[indices]

    return E.replace(lambda expr: type(expr) == Indexed, lambda expr: shift_indexed(expr, idx_shift))

def Dc(E, direction):
    assert(direction == 'x' or direction == 'y' or direction == 'z')
    if direction == 'x':
        shift_hi = (1, 0, 0, 0)
        shift_lo = (-1, 0, 0, 0)
        delta = dx
    elif direction == 'y':
        shift_hi = (0, 1, 0, 0)
        shift_lo = (0, -1, 0, 0)
        delta = dy
    elif direction == 'z':
        shift_hi = (0, 0, 1, 0)
        shift_lo = (0, 0, -1, 0)
        delta = dz
    return (shift(E, shift_hi) - shift(E, shift_lo))/(2 * delta)

def Dc2(E, direction):
    assert(direction == 'x' or direction == 'y' or direction == 'z')
    if direction == 'x':
        shift_hi = (1, 0, 0, 0)
        shift_lo = (-1, 0, 0, 0)
        delta = dx*dx
    elif direction == 'y':
        shift_hi = (0, 1, 0, 0)
        shift_lo = (0, -1, 0, 0)
        delta = dy*dy
    elif direction == 'z':
        shift_hi = (0, 0, 1, 0)
        shift_lo = (0, 0, -1, 0)
        delta = dz*dz
    return (shift(E, shift_hi) - 2*E + shift(E, shift_lo))/(delta)

def DcTen(E, direction):
    retE = ixp.zerorank1(len(E))
    for itr in range(len(E)):
        retE[itr] = Dc(E[itr], direction)
    return retE

def Dc2Ten(E, direction):
    retE = ixp.zerorank1(len(E))
    for itr in range(len(E)):
        retE[itr] = Dc2(E[itr], direction)
    return retE

def grad(phi):
    retGradPhi = ixp.zerorank1(len(directions))
    for itr in range(len(directions)):
        retGradPhi[itr] = Dc(phi, directions[itr])
    return retGradPhi

def Lap(phi):
    return Dc2(phi,'x')+Dc2(phi,'y')+Dc2(phi,'z') 
    
def div(E):
    div = 0
    for itr in range(len(E)):
        div += Dc(E[itr],directions[itr])
    return div

def LapTen(E):
    retLapE = ixp.zerorank1()
    for itr in range(len(directions)):
        retLapE += np.array(Dc2Ten(E,directions[itr]))
    return retLapE

def AMReXcode(expr, varnames, declare_rhs = False, rhsname = "", declare_state = False, statename = ""):
    str_expr = str(expr)
    
    str_expr = str_expr.replace("[","(").replace("]",")")
    str_expr = str_expr.replace("dx","dx[0]")
    str_expr = str_expr.replace("dy","dx[1]")
    str_expr = str_expr.replace("dz","dx[2]")
    str_expr = str_expr.replace("dx[0]**2","(dx[0]*dx[0])")
    str_expr = str_expr.replace("dx[1]**2","(dx[1]*dx[1])")
    str_expr = str_expr.replace("dx[2]**2","(dx[2]*dx[2])")
    str_expr = str_expr.replace("pi","M_PI")
    str_expr = str_expr+";"
    for name in varnames:
        str_expr = str_expr.replace('state_fab'+name,'state_fab')
    for name in varnames:
        str_expr = str_expr.replace(name,"Idx::"+name)
    
    if declare_rhs == True:
        str_expr = "rhs_fab(i, j, k, Idx::"+rhsname+ ") = " + str_expr
        
    if declare_state == True:
        str_expr = "state_fab(i, j, k, Idx::"+statename+ ") = " + str_expr
        
    return str_expr
    
def createSETUP(name, varnames, nghostcells):
    fileSETUP = open(name, "w+")
    fileSETUP.write("#ifndef ET_INTEGRATION_SETUP_K_H \n")
    fileSETUP.write("#define ET_INTEGRATION_SETUP_K_H \n\n")

    fileSETUP.write("#include <AMReX_REAL.H> \n")
    fileSETUP.write("#include <AMReX_Array4.H> \n\n")
    
    fileSETUP.write("namespace Idx { \n")
    fileSETUP.write("         enum ETIndexes {")
    
    Idx_string = ""
    for itr in varnames:
        Idx_string += itr+", "
    Idx_string += "NumScalars"
    
    fileSETUP.write(Idx_string)
    fileSETUP.write("}; \n};\n\n")
    fileSETUP.write("#define NUM_GHOST_CELLS "+str(nghostcells)+"\n\n")
    fileSETUP.write("#endif")

    fileSETUP.close()

def createRHS(name):
    fileRHS = open(name, "w+")
    fileRHS.write("#ifndef ET_INTEGRATION_RHS_K_H \n")
    fileRHS.write("#define ET_INTEGRATION_RHS_K_H \n\n")

    fileRHS.write("#include <AMReX_REAL.H> \n")
    fileRHS.write("#include <AMReX_Array4.H> \n")
    fileRHS.write("#include <ET_Integration_Setup.H> \n\n")

    fileRHS.write("AMREX_GPU_DEVICE \ninline \nvoid \n")
    fileRHS.write("state_rhs(int i, int j, int k, \n")
    fileRHS.write("        amrex::Array4<amrex::Real> const& rhs_fab, \n")
    fileRHS.write("        amrex::Array4<amrex::Real const> const& state_fab, \n")
    fileRHS.write("        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept \n{\n")
    fileRHS.close()

def addRHS(name,RHS):
    fileRHS = open(name,"a+")
    fileRHS.write("         "+RHS+"\n\n")
    fileRHS.close()
    
def finishRHS(name):
    fileRHS = open(name, "a+")
    fileRHS.write("}\n")
    fileRHS.write("#endif")
    fileRHS.close()
    
    
def createINIT(name):
    fileINIT = open(name, "w+")
    fileINIT.write("#ifndef ET_INTEGRATION_INIT_K_H \n")
    fileINIT.write("#define ET_INTEGRATION_INIT_K_H \n\n")

    fileINIT.write("#include <AMReX_REAL.H> \n")
    fileINIT.write("#include <AMReX_Array4.H> \n")
    fileINIT.write("#include <ET_Integration_Setup.H> \n\n")

    fileINIT.write("AMREX_GPU_DEVICE \ninline \nvoid \n")
    fileINIT.write("state_init(int i, int j, int k, \n")
    fileINIT.write("        amrex::Array4<amrex::Real> const& state_fab, \n")
    fileINIT.write("        amrex::Real time, \n")
    fileINIT.write("        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx) noexcept \n{\n")
    fileINIT.write("        amrex::Real x = (i + 0.5)*dx[0];\n")
    fileINIT.write("        amrex::Real y = (j + 0.5)*dx[1];\n")
    fileINIT.write("        amrex::Real z = (k + 0.5)*dx[2];\n\n")
    fileINIT.close()

def addINIT(name,INIT):
    fileINIT = open(name,"a+")
    fileINIT.write("        "+INIT+"\n\n")
    fileINIT.close()
    
def finishINIT(name):
    fileINIT = open(name, "a+")
    fileINIT.write("}\n")
    fileINIT.write("#endif")
    fileINIT.close()
    
    
    
    
    
