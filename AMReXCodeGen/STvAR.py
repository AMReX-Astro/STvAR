import indexedexp as ixp
import numpy as np

import sympy as sp
from sympy import symbols, IndexedBase, Indexed, Idx

import finite_difference as fin
from sympy.printing.cxxcode import *
from sympy.printing.fcode import FCodePrinter

class CustomCXX17Printer(CXX17CodePrinter):
    def _print_Indexed(self, expr):
        return FCodePrinter._print_Indexed(self, expr)
    
printer = CustomCXX17Printer()

def AMReXcode(expr, varnames = "", diagnames = "", initnames = "", convnames = "", customnames = "", declare_rhs = False, rhsname = "", declare_state = False, statename = "", declare_diag = False, diagname = "", declare_init = False, initname = "", declare_conv = False, convname = "", declare_custom = False, customname = "", customprefix = "", customIndex = ""):
    str_expr = str(printer.doprint(expr))

    str_expr = str_expr.replace("pi","M_PI")
    str_expr = str_expr+";"
    for name in varnames:
        str_expr = str_expr.replace("Idx" + name,"Idx::"+name)
        
    for name in diagnames:
        str_expr = str_expr.replace("Diag" + name,"Diag::"+name)
        
    for name in initnames:
        str_expr = str_expr.replace("InitIdx" + name,"InitIdx::"+name)
     
    for name in convnames:
        str_expr = str_expr.replace("ConvIdx" + name, "ConvIdx::"+name)
        
    for name in customnames:
        str_expr = str_expr.replace(customprefix+"Idx" + name, customprefix+"Idx::"+name)
    
    if declare_rhs == True:
        str_expr = "rhs_fab(i, j, k, Idx::"+rhsname+ ") = " + str_expr
        
    if declare_state == True:
        str_expr = "state_fab(i, j, k, Idx::"+statename+ ") = " + str_expr
        
    if declare_diag == True:
        str_expr = "diag(i, j, k, Diag::"+diagname+ ") = " + str_expr
        
    if declare_init == True:
        str_expr = "initial_data(i, j, k, InitIdx::"+initname+ ") = " + str_expr
        
    if declare_conv == True:
        str_expr = "conv_vars_arr(i, j, k, ConvIdx::"+convname+ ") = " + str_expr
        
    if declare_custom == True:
        str_expr = customprefix + "(i, j, k, " + customIndex + "Idx::" + customname+ ") = " + str_expr
        
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
    declState = []
    declDiag = []
    declInit = []
    declConv = []
    declCustom = []
    def __init__(self, symb, expr = 0, Isymb = 0, state = False, rhs = False, diag = False, initial = False, conv = False, custom = False, custpre = "", custIdx = ""):
        self.symb = symbols(symb)
        self.expr = expr
        self.Isymb = Isymb
        self.custpre = custpre
        self.custIdx = custIdx
        if state == True:
            Idxsymb = symbols("Idx" + str(symb))
            self.Isymb = IndexedBase('state_fab')[i,j,k,Idxsymb]
            self.expr = self.Isymb
            stVar.declState.append(str(symb))
        if rhs == True:
            Idxsymb = symbols("Idx" + str(symb))
            self.Isymb = IndexedBase('rhs_fab')[i,j,k,Idxsymb]
        if diag == True:
            Idxsymb = symbols("Diag"+str(symb))
            self.Isymb = IndexedBase('diag')[i,j,k,Idxsymb]
            self.expr = self.Isymb
            stVar.declDiag.append(str(symb))
        if initial == True:
            Idxsymb = symbols("InitIdx"+str(symb))
            self.Isymb = IndexedBase('initial_data')[i,j,k,Idxsymb]
            self.expr = self.Isymb
            stVar.declInit.append(str(symb))
        if conv == True:
            Idxsymb = symbols("ConvIdx"+str(symb))
            self.Isymb = IndexedBase('conv_vars_arr')[i,j,k,Idxsymb]
            self.expr = self.Isymb
            stVar.declConv.append(str(symb))
        if custom == True:
            Idxsymb = symbols(custIdx +"Idx" + str(symb))
            self.Isymb = IndexedBase(custpre)[i,j,k,Idxsymb]
            self.expr = self.Isymb
            stVar.declCustom.append(str(symb))
            
    
    def AMReXSymb2Expr(self):
        return "        amrex::Real " + str(self.symb) + " = " + AMReXcode(self.expr,stVar.declState,stVar.declDiag,stVar.declInit,stVar.declConv,stVar.declCustom)+"\n\n"
    
    def AMReXSymb2State(self):
        return "        amrex::Real " + str(self.symb) + " = " + AMReXcode(self.Isymb,stVar.declState,stVar.declDiag,stVar.declInit,stVar.declConv,stVar.declCustom)+"\n\n"
    
    def AMReXSetRHS(self):
        return "        "+AMReXcode(self.expr,stVar.declState,stVar.declDiag, stVar.declInit,stVar.declConv, stVar.declCustom, declare_rhs = True, rhsname = str(self.symb))+"\n\n"
    
    def AMReXSetState(self):
        return "        "+AMReXcode(self.expr,stVar.declState,stVar.declDiag, stVar.declInit,stVar.declConv, stVar.declCustom, declare_state = True, statename = str(self.symb))+"\n\n"
    
    def AMReXSetDiag(self):
        return "        "+AMReXcode(self.expr,stVar.declState, stVar.declDiag, stVar.declInit,stVar.declConv, stVar.declCustom, declare_diag = True, diagname = str(self.symb))+"\n\n"
    
    def AMReXSetConv(self):
        return "        "+AMReXcode(self.expr,stVar.declState, stVar.declDiag, stVar.declInit,stVar.declConv, stVar.declCustom, declare_conv = True, convname = str(self.symb))+"\n\n"
    
    def AMReXSetCustom(self):
        return "        "+AMReXcode(self.expr,stVar.declState, stVar.declDiag, stVar.declInit,stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb), customprefix = self.custpre, customIndex = self.custIdx)+"\n\n"
            
    
class stVarRank1(stVar):
    def __init__(self, symb, dim = 3, state = False, rhs = False, diag = False, initial = False, conv = False, custom = False, custpre = "", custIdx = ""):
        self.symb = ixp.declarerank1(str(symb),DIM=dim)
        self.expr = ixp.zerorank1(DIM=dim)
        self.Isymb = ixp.zerorank1(DIM=dim)
        self.custpre = custpre
        self.custIdx = custIdx
        self.dim = dim
        if state == True: 
            for itri in range(dim):
                Idxsymb = symbols("Idx" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase('state_fab')[i,j,k,Idxsymb]
                self.expr[itri] = self.Isymb[itri]
                stVar.declState.append(str(self.symb[itri]))
        if rhs == True:
            for itri in range(dim):
                Idxsymb = symbols("Idx" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase('rhs_fab')[i,j,k,Idxsymb]
        if diag == True:
            for itri in range(dim):
                Idxsymb = symbols("Diag" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase('diag')[i,j,k,Idxsymb]
                self.expr[itri] = self.Isymb[itri]
                stVar.declDiag.append(str(self.symb[itri]))
        if initial == True: 
            for itri in range(dim):
                Idxsymb = symbols("InitIdx" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase('initial_data')[i,j,k,Idxsymb]
                self.expr[itri] = self.Isymb[itri]
                stVar.declInit.append(str(self.symb[itri]))
        if conv == True: 
            for itri in range(dim):
                Idxsymb = symbols("ConvIdx" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase('conv_vars_arr')[i,j,k,Idxsymb]
                self.expr[itri] = self.Isymb[itri]
                stVar.declConv.append(str(self.symb[itri]))
        if custom == True: 
            for itri in range(dim):
                Idxsymb = symbols(custIdx+"Idx" + str(self.symb[itri]))
                self.Isymb[itri] = IndexedBase(custpre)[i,j,k,Idxsymb]
                self.expr[itri] = self.Isymb[itri]
                stVar.declCustom.append(str(self.symb[itri]))
                
        self.symb = np.array(self.symb)
        
    def AMReXSymb2Expr(self):
        expr = ""
        for i in range(self.dim):
            expr += "        amrex::Real " + str(self.symb[i]) + " = " + AMReXcode(self.expr[i],stVar.declState, stVar.declDiag,stVar.declInit,stVar.declConv,stVar.declCustom)+'\n'
        return expr
    
    def AMReXSymb2State(self):
        expr = ""
        for i in range(self.dim):
            expr += "        amrex::Real " + str(self.symb[i]) + " = " + AMReXcode(self.Isymb[i],stVar.declState, stVar.declDiag,stVar.declInit,stVar.declConv,stVar.declCustom)+'\n'
        return expr
            
    def AMReXSetRHS(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.expr[i], stVar.declState, stVar.declDiag, stVar.declInit,stVar.declConv,stVar.declCustom, declare_rhs = True, rhsname = str(self.symb[i]))+"\n\n"
        return expr
    
    def AMReXSetState(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.expr[i], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv,stVar.declCustom, declare_state = True, statename = str(self.symb[i]))+"\n\n"
        return expr
    
    def AMReXSetDiag(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.expr[i], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv,stVar.declCustom, declare_diag = True, diagname = str(self.symb[i]))+"\n\n"
        return expr
    
    def AMReXSetConv(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.expr[i], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv,stVar.declCustom, declare_conv = True, convname = str(self.symb[i]))+"\n\n"
        return expr
    
    def AMReXSetCustom(self):
        expr = ""
        for i in range(self.dim):
            expr += "        "+AMReXcode(self.expr[i], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb[i]), customprefix = self.custpre, customIndex = self.custIdx)+"\n\n"
        return expr
    
    
class stVarRank2(stVar):
    def __init__(self, symb, sym = 'nosym', dim = 3, state = False, rhs = False, diag = False, initial = False, conv = False, custom = False, custpre = "", custIdx = "", resetsym = True):
        self.symb = ixp.declarerank2(str(symb), sym, DIM=dim)
        self.expr = ixp.zerorank2(DIM=dim)
        self.Isymb = ixp.zerorank2(DIM=dim)
        self.custpre = custpre
        self.custIdx = custIdx
        self.dim = dim
        self.sym = sym
        if state == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols("Idx" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase('state_fab')[i,j,k,Idxsymb]
                    self.expr[itri][itrj] = self.Isymb[itri][itrj]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.declState.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.declState.append(str(self.symb[itri][itrj]))
        if rhs == True:
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols("Idx" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase('rhs_fab')[i,j,k,Idxsymb]
                        
        if diag == True:
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols("Diag" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase('diag')[i,j,k,Idxsymb]
                    self.expr[itri][itrj] = self.Isymb[itri][itrj]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.declDiag.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.declDiag.append(str(self.symb[itri][itrj]))
        
        if initial == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols("InitIdx" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase('initial_data')[i,j,k,Idxsymb]
                    self.expr[itri][itrj] = self.Isymb[itri][itrj]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.declInit.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.declInit.append(str(self.symb[itri][itrj]))
                        
        if conv == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols("ConvIdx" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase('conv_vars_arr')[i,j,k,Idxsymb]
                    self.expr[itri][itrj] = self.Isymb[itri][itrj]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.declConv.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.declConv.append(str(self.symb[itri][itrj]))
                        
        if custom == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    Idxsymb = symbols(custIdx + "Idx" + str(self.symb[itri][itrj]))
                    self.Isymb[itri][itrj] = IndexedBase(custpre)[i,j,k,Idxsymb]
                    self.expr[itri][itrj] = self.Isymb[itri][itrj]
                    
            if sym == 'sym01':
                for itri in range(dim):
                    for itrj in range(itri,dim):
                        stVar.declCustom.append(str(self.symb[itri][itrj]))
            else:
                for itri in range(dim):
                    for itrj in range(dim):
                        stVar.declCustom.append(str(self.symb[itri][itrj]))
        
        if resetsym == True:
                self.symb = ixp.declarerank2(str(symb),'nosym', DIM=dim)
                
        self.symb = np.array(self.symb)    
    
    def AMReXSymb2Expr(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                expr += "        amrex::Real " + str(self.symb[i][j]) + " = " + AMReXcode(self.expr[i][j],stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
    
    def AMReXSymb2State(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                expr += "        amrex::Real " + str(self.symb[i][j]) + " = " + AMReXcode(self.Isymb[i][j],stVar.declState, stVar.declDiag,stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
            
    def AMReXSetRHS(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_rhs = True, rhsname = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i,self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, declare_rhs = True, rhsname = str(self.symb[i][j]))+"\n\n"
        return expr
    
    def AMReXSetState(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_state = True, statename = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i, self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_state = True, statename = str(self.symb[i][j]))+"\n\n"
        return expr
    
    def AMReXSetDiag(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_diag = True, diagname = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i, self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_diag = True, diagname = str(self.symb[i][j]))+"\n\n"
        return expr
    
    def AMReXSetConv(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_conv = True, convname = str(self.symb[i][j]))+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i, self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_conv = True, convname = str(self.symb[i][j]))+"\n\n"
        return expr
    
    def AMReXSetCustom(self):
        expr = ""
        if self.sym == 'nosym':
            for i in range(self.dim):
                for j in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb[i][j]), customprefix = custpre, customIndex = custIdx)+"\n\n"
        elif self.sym == 'sym01':
            for i in range(self.dim):
                for j in range(i, self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb[i][j]), customprefix = self.custpre, customIndex = self.custIdx)+"\n\n"
        return expr
    

class stVarRank3(stVar):
    def __init__(self, symb, sym = 'nosym', dim = 3, state = False, rhs = False, diag = False, initial = False, conv = False, custom = False, custpre = "", custIdx = ""):
        self.symb = ixp.declarerank3(str(symb), sym, DIM=dim)
        self.expr = ixp.zerorank3(DIM=dim)
        self.Isymb = ixp.zerorank3(DIM=dim)
        self.custpre = custpre
        self.custIdx = custIdx
        self.dim = dim
        if state == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols("Idx" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase('state_fab')[i,j,k,Idxsymb]
                        self.expr[itri][itrj][itrk] = self.Isymb[itri][itrj][itrk]
                        stVar.declState.append(str(self.symb[itri][itrj][itrk]))
        if rhs == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols("Idx" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase('rhs_fab')[i,j,k,Idxsymb]
        if diag == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols("Diag" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase('diag')[i,j,k,Idxsymb]
                        self.expr[itri][itrj][itrk] = self.Isymb[itri][itrj][itrk]
                        stVar.declDiag.append(str(self.symb[itri][itrj][itrk]))
        
        if initial == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols("InitIdx" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase('initial_data')[i,j,k,Idxsymb]
                        self.expr[itri][itrj][itrk] = self.Isymb[itri][itrj][itrk]
                        stVar.declInit.append(str(self.symb[itri][itrj][itrk]))
                        
        if conv == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols("ConvIdx" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase('conv_vars_arr')[i,j,k,Idxsymb]
                        self.expr[itri][itrj][itrk] = self.Isymb[itri][itrj][itrk]
                        stVar.declConv.append(str(self.symb[itri][itrj][itrk]))
                      
        if custom == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        Idxsymb = symbols(custIdx + "Idx" + str(self.symb[itri][itrj][itrk]))
                        self.Isymb[itri][itrj][itrk] = IndexedBase(custpre)[i,j,k,Idxsymb]
                        self.expr[itri][itrj][itrk] = self.Isymb[itri][itrj][itrk]
                        stVar.declCustom.append(str(self.symb[itri][itrj][itrk]))
                
        self.symb = np.array(self.symb)
        
    def AMReXSymb2Expr(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        amrex::Real " + str(self.symb[i][j][k]) + " = " + AMReXcode(self.expr[i][j][k],stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
    
    def AMReXSymb2State(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        amrex::Real " + str(self.symb[i][j][k]) + " = " + AMReXcode(self.Isymb[i][j][k],stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
            
    def AMReXSetRHS(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j][k], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_rhs = True, rhsname = str(self.symb[i][j][k]))+"\n\n"
        return expr
    
    def AMReXSetState(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j][k], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_state = True, statename = str(self.symb[i][j][k]))+"\n\n"
        return expr
    
    def AMReXSetDiag(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j][k], stVar.declState, stVar.declDiag, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_diag = True, diagname = str(self.symb[i][j][k]))+"\n\n"
        return expr 
    
    def AMReXSetConv(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j][k], stVar.declState, stVar.declDiag, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_conv = True, convname = str(self.symb[i][j][k]))+"\n\n"
        return expr
    
    def AMReXSetCustom(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    expr += "        "+AMReXcode(self.expr[i][j][k], stVar.declState, stVar.declDiag, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb[i][j][k]), customprefix = self.custpre, customindex = self.custIdx)+"\n\n"
        return expr
    
    
class stVarRank4(stVar):
    def __init__(self, symb, sym = 'none', dim = 3, state = False, rhs = False, diag = False, initial = False, conv = False, custom = False, custpre = "", custIdx = ""):
        self.symb = ixp.declarerank4(str(symb),sym,DIM=dim)
        self.expr = ixp.zerorank4(DIM=dim)
        self.Isymb = ixp.zerorank4(DIM=dim)
        self.custpre = custpre
        self.custIdx = custIdx
        self.dim = dim
        if state == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols("Idx" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase('state_fab')[i,j,k,Idxsymb]
                            self.expr[itri][itrj][itrk][itrl] = self.Isymb[itri][itrj][itrk][itrl]
                            stVar.declState.append(str(self.symb[itri][itrj][itrk][itrl]))
        if rhs == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols("Idx" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase('rhs_fab')[i,j,k,Idxsymb]
        if diag == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols("Diag" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase('diag')[i,j,k,Idxsymb]
                            self.expr[itri][itrj][itrk][itrl] = self.Isymb[itri][itrj][itrk][itrl]
                            stVar.declDiag.append(str(self.symb[itri][itrj][itrk][itrl]))
        
        if initial == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols("InitIdx" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase('initial_data')[i,j,k,Idxsymb]
                            self.expr[itri][itrj][itrk][itrl] = self.Isymb[itri][itrj][itrk][itrl]
                            stVar.declInit.append(str(self.symb[itri][itrj][itrk][itrl]))
                            
        if conv == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols("ConvIdx" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase('conv_vars_arr')[i,j,k,Idxsymb]
                            self.expr[itri][itrj][itrk][itrl] = self.Isymb[itri][itrj][itrk][itrl]
                            stVar.declConv.append(str(self.symb[itri][itrj][itrk][itrl]))

                            
        if custom == True: 
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            Idxsymb = symbols(custIdx+ "Idx" + str(self.symb[itri][itrj][itrk][itrl]))
                            self.Isymb[itri][itrj][itrk][itrl] = IndexedBase(custpre)[i,j,k,Idxsymb]
                            self.expr[itri][itrj][itrk][itrl] = self.Isymb[itri][itrj][itrk][itrl]
                            stVar.declCustom.append(str(self.symb[itri][itrj][itrk][itrl]))
                            
        self.symb = np.array(self.symb)
        
    def AMReXSymb2Expr(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        amrex::Real " + str(self.symb[i][j][k][l]) + " = " + AMReXcode(self.expr[i][j][k][l],stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
    
    def AMReXSymb2State(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        amrex::Real " + str(self.symb[i][j][k][l]) + " = " + AMReXcode(self.Isymb[i][j][k][l],stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom)+'\n'
        return expr
            
    def AMReXSetRHS(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.expr[i][j][k][l], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_rhs = True, rhsname = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
    
    def AMReXSetState(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.expr[i][j][k][l], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_state = True, statename = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
    
    def AMReXSetDiag(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.expr[i][j][k][l], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_diag = True, diagname = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
    
    def AMReXSetConv(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.expr[i][j][k][l], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_conv = True, convname = str(self.symb[i][j][k][l]))+"\n\n"
        return expr
    
    def AMReXSetCustom(self):
        expr = ""
        for i in range(self.dim):
            for j in range(self.dim):
                for k in range(self.dim):
                    for l in range(self.dim):
                        expr += "        "+AMReXcode(self.expr[i][j][k][l], stVar.declState, stVar.declDiag, stVar.declInit, stVar.declConv, stVar.declCustom, declare_custom = True, customname = str(self.symb[i][j][k][l]), customprefix = self.custpre, customIndex = self.custIdx)+"\n\n"
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

def Ds(var, difftype = '1'):
    return Dsymb(var.symb, difftype)

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

def Dv(E, directions = '0', CnUpDn = 'cn', order = 2):
    return Dvar(E.var, directions, CnUpDn, order)

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
    
    
def DstVar(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, DIM = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVar('stvarOut')
        stvarOut.expr = FullKOdiss(stVarIn.expr, orderD, sigma, DIM)
        stvarOut.Isymb = stvarOut.expr
        stvarOut.symb = symbols('KO'+str(stVarIn.symb))
    elif difftype == 1:
        stvarOut = stVarRank1('stvarOut', dim = DIM)
        for i in range(DIM):
            stvarOut.expr[i] = Dvar(stVarIn.expr, directions = str(i), CnUpDn = CnUpDnRank1, order = orderD)
            stvarOut.Isymb[i] = stvarOut.expr[i]
            stvarOut.symb[i] = Dsymb(stVarIn.symb, difftype = str(i))
    elif difftype == 2:
        stvarOut = stVarRank2('stvarOut', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                stvarOut.expr[i][j] = Dvar(stVarIn.expr, directions = str(i)+str(j), CnUpDn = CnUpDnRank1, order = orderD)
                stvarOut.Isymb[i][j] = stvarOut.expr[i][j]
                stvarOut.symb[i][j] = Dsymb(stVarIn.symb, difftype = str(i)+str(j))
    return stvarOut

def DstVarRank1(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, DIM = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVarRank1('stvarout', dim = DIM)
        for i in range(DIM):
            stvarOut.expr[i] = FullKOdiss(stVarIn.expr[i], orderD, sigma, DIM)
            stvarOut.Isymb[i] = stvarOut.expr[i]
            stvarOut.symb[i] = symbols('KO'+str(stVarIn.symb[i]))
    elif difftype == 1:
        stvarOut = stVarRank2('stvarOut', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                stvarOut.expr[i][j] = Dvar(stVarIn.expr[i], directions = str(j), CnUpDn = CnUpDnRank1, order = orderD)
                stvarOut.Isymb[i][j] = stvarOut.expr[i][j]
                stvarOut.symb[i][j] = Dsymb(stVarIn.symb[i], difftype = str(j))
                
    elif difftype == 2:
        stvarOut = stVarRank3('stvarOut', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    stvarOut.expr[i][j][k] = Dvar(stVarIn.expr[i], directions = str(j)+str(k), CnUpDn = CnUpDnRank1, order = orderD)
                    stvarOut.Isymb[i][j][k] = stvarOut.expr[i][j][k]
                    stvarOut.symb[i][j][k] = Dsymb(stVarIn.symb[i], difftype = str(j)+str(k))
    return stvarOut

def DstVarRank2(stVarIn, difftype = 1, CnUpDnRank1 = 'cn', orderD = 2, DIM = 3, sigma = 0.1):
    if difftype == 'KO':
        stvarOut = stVarRank2('stvarout', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                stvarOut.expr[i][j] = FullKOdiss(stVarIn.expr[i][j], orderD, sigma, DIM)
                stvarOut.Isymb[i][j] = stvarOut.expr[i][j]
                stvarOut.symb[i][j] = symbols('KO'+str(stVarIn.symb[i][j]))
    elif difftype == 1:
        stvarOut = stVarRank3('stvarOut', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    stvarOut.expr[i][j][k] = Dvar(stVarIn.expr[i][j], directions = str(k), CnUpDn = CnUpDnRank1, order = orderD)
                    stvarOut.Isymb[i][j][k] = stvarOut.expr[i][j][k]
                    stvarOut.symb[i][j][k] = Dsymb(stVarIn.symb[i][j], difftype = str(k))
                
    elif difftype == 2:
        stvarOut = stVarRank4('stvarOut', dim = DIM)
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    for l in range(DIM):
                        stvarOut.expr[i][j][k][l] = Dvar(stVarIn.expr[i][j], directions = str(k)+str(l), CnUpDn = CnUpDnRank1, order = orderD)
                        stvarOut.Isymb[i][j][k][l] = stvarOut.expr[i][j][k][l]
                        stvarOut.symb[i][j][k][l] = Dsymb(stVarIn.symb[i][j], difftype = str(k)+str(l))
    return stvarOut    
    
    
def TagCondition(TagVar, conditional = '>', ExprOrSymb = 'expr'):
    String =  stVar('String')
    if ExprOrSymb == 'expr':
        String.expr = str(TagVar.expr) + ' ' + conditional + ' ' + 'error_threshold'
    elif ExprOrSymb == 'symb':
        String.expr = str(TagVar.symb) + ' ' + conditional + ' ' + 'error_threshold'
    return '        return ' + AMReXcode(String.expr,stVar.declState) + '\n'    
    
    
    

