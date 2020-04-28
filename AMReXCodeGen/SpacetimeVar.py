import indexedexp as ixp
import numpy as np

import sympy as sp
from sympy import symbols, IndexedBase, Indexed, Idx #, preorder_traversal

import finite_difference as fin
from sympy.printing.cxxcode import *
from sympy.printing.fcode import FCodePrinter

class CustomCXX17Printer(CXX17CodePrinter):
    def _print_Indexed(self, expr):
        return FCodePrinter._print_Indexed(self, expr)
    
printer = CustomCXX17Printer()

def AMReXcode(expr, varnames= "", declare_rhs = False, rhsname = "", declare_state = False, statename = "", declare_diag = False, diagname = ""):
    str_expr = str(printer.doprint(expr))

    str_expr = str_expr.replace("pi","M_PI")
    str_expr = str_expr+";"
    for name in varnames:
        str_expr = str_expr.replace(name,"Idx::"+name)
    
    if declare_rhs == True:
        str_expr = "rhs_fab(i, j, k, Idx::"+rhsname+ ") = " + str_expr
        
    if declare_state == True:
        str_expr = "state_fab(i, j, k, Idx::"+statename+ ") = " + str_expr
        
    if declare_diag == True:
        str_expr = "diag(i, j, k, Diag::"+diagname+ ") = " + str_expr
        
    return str_expr
    
Nx, Ny, Nz, Nn= symbols('Nx Ny Nz Nn', integer=True)
i = Idx('i', Nx)
j = Idx('j', Ny)
k = Idx('k', Nz)
n = Idx('n', Nn)

DX = ixp.zerorank1(DIM=3)
DX[0], DX[1], DX[2] = symbols('dx[0] dx[1] dx[2]')

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

class stVar:
    decl = []
    declDiag = []
    def __init__(self, symb, var = 0, dim = 3, declare = False, declareDiag = False):
        self.symb = symbols(symb)
        self.var = symbols('0')
        if declare == True: 
            self.var = IndexedBase('state_fab')[i,j,k,str(symb)]
            stVar.decl.append(str(symb))
        if declareDiag == True:
            stVar.declDiag.append(str(symb))
            
    def AMReXDeclare(self):
        return "        amrex::Real " + str(self.symb) + " = " + AMReXcode(self.var,stVar.decl)+'\n'
    
    def AMReXReal(self):
        return "        amrex::Real " + str(self.symb) + " = " + AMReXcode(self.var)+'\n'
    
    def AMReXRHS(self):
        return "        "+AMReXcode(self.var,declare_rhs = True, rhsname = str(self.symb))+"\n\n"
    
    def AMReXInit(self):
        return "        "+AMReXcode(self.var,declare_state = True, statename = str(self.symb))+"\n\n"
    
    def AMReXDiag(self):
        return "        "+AMReXcode(self.var,declare_diag = True, diagname = str(self.symb))+"\n\n"
            
    
class stVarRank1(stVar):
    def __init__(self, symb, dim = 3, declare = False):
        self.symb = ixp.declarerank1(str(symb),DIM=dim)
        self.var = ixp.zerorank1(DIM=dim)
        self.dim = dim
        if declare == True: 
            for itri in range(dim):
                self.var[itri] = IndexedBase('state_fab')[i,j,k,str(self.symb[itri])]
                stVar.decl.append(str(self.symb[itri]))
                
        self.symb = np.array(self.symb)
        
    def AMReXDeclare(self):
        expr = ""
        for i in range(self.dim):
            expr += "        amrex::Real " + str(self.symb[i]) + " = " + AMReXcode(self.var[i],stVar.decl)+'\n'
        return expr
    
    def AMReXReal(self):
        expr = ""
        for i in range(self.dim):
            expr += "        amrex::Real " + str(self.symb[i]) + " = " + AMReXcode(self.var[i])+'\n'
        return expr
            
    def AMReXRHS(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.var[i],declare_rhs = True, rhsname = str(self.symb[i]))+"\n\n"
        return expr
    
    def AMReXInit(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.var[i],declare_state = True, statename = str(self.symb[i]))+"\n\n"
        return expr
        
class stVarRank2(stVar):
    def __init__(self, symb, sym = 'nosym', dim = 3, declare = False, resetsym = True):
        self.symb = ixp.declarerank2(str(symb), sym, DIM=dim)
        self.var = ixp.zerorank2(DIM=dim)
        self.dim = dim
        self.sym = sym
        if declare == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    self.var[itri][itrj] = IndexedBase('state_fab')[i,j,k,str(self.symb[itri][itrj])]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.decl.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.decl.append(str(self.symb[itri][itrj]))
        if resetsym == True:
                self.symb = ixp.declarerank2(str(symb),'nosym', DIM=dim)
                
        self.symb = np.array(self.symb)
        
    def AMReXDeclare(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                expr += "        amrex::Real " + str(self.symb[i][j]) + " = " + AMReXcode(self.var[i][j],stVar.decl)+'\n'
        return expr
    
    def AMReXReal(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                expr += "        amrex::Real " + str(self.symb[i][j]) + " = " + AMReXcode(self.var[i][j])+'\n'
        return expr
    
    def AMReXRHS(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.var[i][j],declare_rhs = True, rhsname = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i,self.dim):
                    expr += "        "+AMReXcode(self.var[i][j],declare_rhs = True, rhsname = str(self.symb[i][j]))+"\n\n"
        return expr
    
    def AMReXInit(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.var[i][j],declare_state = True, statename = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i,self.dim):
                    expr += "        "+AMReXcode(self.var[i][j],declare_state = True, statename = str(self.symb[i][j]))+"\n\n"
        return expr

class stVarRank3(stVar):
    def __init__(self, symb, sym = 'nosym', dim = 3):
        self.symb = ixp.declarerank3(str(symb), sym, DIM=dim)
        self.var = ixp.zerorank3(DIM=dim)
        self.dim = dim
        self.symb = np.array(self.symb)
        
    def AMReXDeclare(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        amrex::Real " + str(self.symb[i][j][k]) + " = " + AMReXcode(self.var[i][j][k],stVar.decl)+'\n'
        return expr
    
    def AMReXReal(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        amrex::Real " + str(self.symb[i][j][k]) + " = " + AMReXcode(self.var[i][j][k])+'\n'
        return expr
    
    def AMReXRHS(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.var[i][j][k],declare_rhs = True, rhsname = str(self.symb[i][j][k]))+"\n\n"
        return expr
    
    def AMReXInit(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.var[i][j][k],declare_state = True, statename = str(self.symb[i][j][k]))+"\n\n"
        return expr
    
    
        
class stVarRank4(stVar):
    def __init__(self, symb, sym = 'none', dim = 3):
        self.symb = ixp.declarerank4(str(symb),sym, DIM=dim)
        self.var = ixp.zerorank4(DIM=dim)
        self.dim = dim
        
        self.symb = np.array(self.symb)
        
    def AMReXDeclare(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        amrex::Real " + str(self.symb[i][j][k][l]) + " = " + AMReXcode(self.var[i][j][k][l],stVar.decl)+'\n'
        return expr
    
    def AMReXReal(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        amrex::Real " + str(self.symb[i][j][k][l]) + " = " + AMReXcode(self.var[i][j][k][l])+'\n'
        return expr
    
    def AMReXRHS(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.var[i][j][k][l],declare_rhs = True, rhsname = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
    
    def AMReXInit(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.var[i][j][k][l],declare_state = True, statename = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
        
def Dsymb(var, difftype = '1'):
    if difftype == 'KO':
        dvar = symbols('KO'+str(var))
    else:
        if len(str(difftype)) == 1:
            dvar = symbols('d'+str(var)+str(difftype))
        elif len(difftype) == 2:
            dvar = symbols('dd'+str(var)+str(difftype))
    return dvar

def Dvar(E, directions = '0', CnUpDn = 'cn', order = 2):
    difftype = 'd'
    if CnUpDn == 'up':
        difftype += 'up'
    elif CnUpDn == 'dn':
        difftype += 'dn'
    
    delta = DX[int(directions[0])]
    
    if len(directions) == 1:
        difftype += 'D' + directions
    elif len(directions) == 2:
        difftype += 'DD' + directions
        delta *= DX[int(directions[1])]
    
    fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl(difftype,FDORDER=order)
    shiftE = 0
    for i in range(len(fdcoeffs)):
        shiftE += fdcoeffs[i]*shift(E,fdstencl[i])
    shiftE = shiftE/delta
    return shiftE

def Diffup1(E, direction, order):
    fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl('dD'+direction,FDORDER=order)
    delta = DX[int(direction)]
        
    shiftE = 0
    for i in range(len(fdcoeffs)):
        shiftE += fdcoeffs[i]*shift(E,fdstencl[i])
    shiftE = shiftE/delta
    return shiftE

def Diffdn1(E, direction, order):
    fdcoeffs, fdstencl = fin.compute_fdcoeffs_fdstencl('ddnD'+direction,FDORDER=order)
    delta = DX[int(direction)]
        
    shiftE = 0
    for i in range(len(fdcoeffs)):
        shiftE += fdcoeffs[i]*shift(E,fdstencl[i])
    shiftE = shiftE/delta
    return shiftE

def KOdiss(E, direction, order, sigma=0.1):
    delta = DX[int(direction)]
    r = int((2+order)/2)
    for i in range(r):
        E = sp.simplify(Diffdn1(E,direction,1))
    for i in range(r):
        E = sp.simplify(Diffup1(E,direction,1))        
    E = (-1)**(r+1)/(2**(2*r))*delta**(2*r-1)*sigma*E
    return E

def FullKOdiss(E, order, sigma = 0.1, dim = 3):
    KOD = 0
    for i in range(dim):
        KOD += KOdiss(E, str(i), order, sigma)
    return KOD
        

def DstVar(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, dim = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVar('stvarOut')
        stvarOut.var = FullKOdiss(stVarIn.var, orderD, sigma, dim)
        stvarOut.symb = symbols('KO'+str(stVarIn.symb))
    elif difftype == 1:
        stvarOut = stVarRank1('stvarOut')
        for i in range(dim):
            stvarOut.var[i] = Dvar(stVarIn.var, directions = str(i), CnUpDn = CnUpDnRank1, order = orderD)
            stvarOut.symb[i] = Dsymb(stVarIn.symb, difftype = str(i))
    elif difftype == 2:
        stvarOut = stVarRank2('stvarOut')
        for i in range(dim):
            for j in range(dim):
                stvarOut.var[i][j] = Dvar(stVarIn.var, directions = str(i)+str(j), CnUpDn = CnUpDnRank1, order = orderD)
                stvarOut.symb[i][j] = Dsymb(stVarIn.symb, difftype = str(i)+str(j))
    return stvarOut

def DstVarRank1(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, dim = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVarRank1('stvarout')
        for i in range(dim):
            stvarOut.var[i] = FullKOdiss(stVarIn.var[i], orderD, sigma, dim)
            stvarOut.symb[i] = symbols('KO'+str(stVarIn.symb[i]))
    elif difftype == 1:
        stvarOut = stVarRank2('stvarOut')
        for i in range(dim):
            for j in range(dim):
                stvarOut.var[i][j] = Dvar(stVarIn.var[i], directions = str(j), CnUpDn = CnUpDnRank1, order = orderD)
                stvarOut.symb[i][j] = Dsymb(stVarIn.symb[i], difftype = str(j))
                
    elif difftype == 2:
        stvarOut = stVarRank3('stvarOut')
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    stvarOut.var[i][j][k] = Dvar(stVarIn.var[i], directions = str(j)+str(k), CnUpDn = CnUpDnRank1, order = orderD)
                    stvarOut.symb[i][j][k] = Dsymb(stVarIn.symb[i], difftype = str(j)+str(k))
    return stvarOut

def DstVarRank2(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, dim = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVarRank2('stvarout')
        for i in range(dim):
            for j in range(dim):
                stvarOut.var[i][j] = FullKOdiss(stVarIn.var[i][j], orderD, sigma, dim)
                stvarOut.symb[i][j] = symbols('KO'+str(stVarIn.symb[i][j]))
    elif difftype == 1:
        stvarOut = stVarRank3('stvarOut')
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    stvarOut.var[i][j][k] = Dvar(stVarIn.var[i][j], directions = str(k), CnUpDn = CnUpDnRank1, order = orderD)
                    stvarOut.symb[i][j][k] = Dsymb(stVarIn.symb[i][j], difftype = str(k))
                
    elif difftype == 2:
        stvarOut = stVarRank4('stvarOut')
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        stvarOut.var[i][j][k][l] = Dvar(stVarIn.var[i][j], directions = str(k)+str(l), CnUpDn = CnUpDnRank1, order = orderD)
                        stvarOut.symb[i][j][k][l] = Dsymb(stVarIn.symb[i][j], difftype = str(k)+str(l))
    return stvarOut
