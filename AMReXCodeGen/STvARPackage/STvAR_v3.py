import numpy as np

import sympy as sp
from sympy import symbols, IndexedBase, Indexed, Idx

from sympy.printing.cxx import *
from sympy.printing.fortran import FCodePrinter

class CustomCXX17Printer(CXX17CodePrinter):
    def _print_Indexed(self, expr):
        return FCodePrinter._print_Indexed(self, expr)
    
printer = CustomCXX17Printer()

import STvARIndexing as idx
import STvARSymbolParsing as ssp
import STvARFiniteDifference as fd

def amrexcode(expr, varnames=[]):
    str_expr = str(printer.doprint(expr))
    str_expr = str_expr.replace("pi", "M_PI")
    for var in varnames:
        str_expr = str_expr.replace(var.vartype+var.varname,var.vartype + "Idx::" + var.varname)

    return str_expr

Nx, Ny, Nz, Nn= symbols('Nx Ny Nz Nn', integer=True)
i = Idx('i', Nx)
j = Idx('j', Ny)
k = Idx('k', Nz)
n = Idx('n', Nn)

class barevar:
    def __init__(self, varname, vartype, varprefix):
        self.varname = varname
        self.vartype = vartype
        self.varprefix = varprefix
        
class stvar:
    gridvars = []
    vartypes = []
    def __init__(self, symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", isreal = False):
        self.symb = symbols(symb, real = isreal)
        self.expr = expr
        self.isymb = isymb
        self.vartype = vartype
        self.rootname = symb
        self.rank = ''
        
        if(varprefix == ""):
            self.varprefix = vartype
            
        else:
            self.varprefix = varprefix
        
        if gridvar == True:
            var = barevar(str(self.symb), self.vartype, self.varprefix)
            Idxsymb = symbols(self.vartype+str(self.symb))
            self.isymb = IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
            if addtolist:
                stvar.gridvars.append(var)
                if self.vartype not in stvar.vartypes:
                    stvar.vartypes.append(self.vartype)
            
    def symb2expr(self):
        return "        amrex::Real " + str(self.symb) + " = " + amrexcode(self.expr,stvar.gridvars)+";\n\n"
    
    def symb2isymb(self):
        return "        amrex::Real " + str(self.symb) + " = " + amrexcode(self.isymb,stvar.gridvars)+";\n\n"
    
    def setisymb(self, prefix = ""):
        if prefix == "":
            rstring = "        " + amrexcode(self.isymb,stvar.gridvars) + " = " + amrexcode(self.expr,stvar.gridvars)+";\n\n"
        else:
            Idxsymb = symbols(self.vartype+str(self.symb))
            tempisymb = IndexedBase(prefix)[i,j,k,Idxsymb]
            rstring = "        " + amrexcode(tempisymb, stvar.gridvars) + " = " + amrexcode(self.expr,stvar.gridvars)+ ";\n\n"
        
        return rstring
    
    def diff(self, dstring, Accuracy = 2, dim = 3):
        dsymbrootname, dsymbrank = fd.dsymb_name_rank(self, dstring)
        if 'DD' in dstring:
            stvarout = stvarrank2(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = 'sym01', dim = dim)
            for i in range(dim):
                for j in range(dim):
                    stvarout.isymb[i][j] = fd.disymb(self.isymb, dstring + str(i) + str(j), Accuracy)
        elif 'dKODFull' in dstring:
            stvarout = stvar(ssp.unparse_name_rank(dsymbrootname, dsymbrank))
            for i in range(dim):
                stvarout.isymb += fd.disymb(self.isymb, 'dKOD'+str(i), Accuracy)
        else:
            stvarout = stvarrank1(ssp.unparse_name_rank(dsymbrootname, dsymbrank), dim = dim)
            for i in range(dim):
                stvarout.isymb[i] = fd.disymb(self.isymb, dstring +str(i), Accuracy)
        stvarout.expr = stvarout.isymb
        return stvarout
    
class stvarrank1:
    def __init__(self, symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", isreal = False, dim = 3):
        self.symb = idx.declarerank1(str(symb),DIM=dim)
        self.expr = idx.zerorank1(DIM=dim)
        self.isymb = idx.zerorank1(DIM=dim)
        self.vartype = vartype
        self.rootname, self.rank = ssp.parse_name_rank(symb)
        if len(self.rank) != 1:
            print("The rank must 1.")
        self.dim = dim
        
        if(varprefix == ""):
            self.varprefix = vartype
            
        else:
            self.varprefix = varprefix
            
        if gridvar == True:
            for itri in range(dim):
                var = barevar(str(self.symb[itri]), self.vartype, self.varprefix)
                Idxsymb = symbols(self.vartype+str(self.symb[itri]))
                self.isymb[itri] = IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                if addtolist:
                    stvar.gridvars.append(var)
                    if self.vartype not in stvar.vartypes:
                        stvar.vartypes.append(self.vartype)
                        
        self.symb = np.array(self.symb)
    
    def symb2expr(self):
        expr = ""
        for itri in range(self.dim):
            if self.symb[itri] != 0:
                expr += "        amrex::Real " + str(self.symb[itri]) + " = " + amrexcode(self.expr[itri],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def symb2isymb(self):
        expr = ""
        for itri in range(self.dim):
            if self.symb[itri] != 0:
                expr += "        amrex::Real " + str(self.symb[itri]) + " = " + amrexcode(self.isymb[itri],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def setisymb(self, prefix = ""):
        expr = ""
        for itri in range(self.dim):
            if self.isymb[itri] != 0: 
                if prefix == "":
                    expr += "        " + amrexcode(self.isymb[itri],stvar.gridvars) + " = " + amrexcode(self.expr[itri],stvar.gridvars)+";\n"
                else:
                    Idxsymb = symbols(self.vartype+str(self.symb[itri]))
                    tempisymb = IndexedBase(prefix)[i,j,k,Idxsymb]
                    expr += "        " + amrexcode(tempisymb, stvar.gridvars) + " = " + amrexcode(self.expr[itri],stvar.gridvars)+ ";\n"
        expr += '\n'
        return expr
    
    def diff(self, dstring, Accuracy = 2, dim = 3, makesparse = True):
        dsymbrootname, dsymbrank = fd.dsymb_name_rank(self, dstring)
        if 'DD' in dstring:
            stvarout = stvarrank3(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = 'sym12', dim = dim)
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        stvarout.isymb[i][j][k] = fd.disymb(self.isymb[i], dstring + str(j) + str(k), Accuracy)
                        if (makesparse == True) and (self.symb[i] == 0):
                            stvarout.symb[i][j][k] = 0
                            stvarout.isymb[i][j][k] = 0
        elif 'dKODFull' in dstring:
            stvarout = stvarrank1(ssp.unparse_name_rank(dsymbrootname, dsymbrank), dim = dim)
            for i in range(dim):
                for j in range(dim):
                    stvarout.isymb[i] += fd.disymb(self.isymb[i], 'dKOD'+str(j), Accuracy)
                    if (makesparse == True) and (self.symb[i] == 0):
                        stvarout.symb[i] = 0
                        stvarout.isymb[i] = 0
        else:
            stvarout = stvarrank2(ssp.unparse_name_rank(dsymbrootname, dsymbrank), dim = dim)
            for i in range(dim):
                for j in range(dim):
                    stvarout.isymb[i][j] = fd.disymb(self.isymb[i], dstring + str(j), Accuracy)
                    if (makesparse == True) and (self.symb[i] == 0):
                        stvarout.symb[i][j] = 0
                        stvarout.isymb[i][j] = 0
        stvarout.expr = stvarout.isymb
        return stvarout
    
class stvarrank2:
    def __init__(self, symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", isreal = False, dim = 3, sym = 'nosym', resetsym = True):
        self.symb = idx.declarerank2(str(symb), sym, DIM=dim)
        self.expr = idx.zerorank2(DIM=dim)
        self.isymb = idx.zerorank2(DIM=dim)
        self.vartype = vartype
        self.rootname, self.rank = ssp.parse_name_rank(symb)
        if len(self.rank) != 2:
            print("The rank must 2.")
        self.dim = dim
        self.sym = sym
        
        if(varprefix == ""):
            self.varprefix = vartype
            
        else:
            self.varprefix = varprefix
            
        if gridvar == True:
            for itri in range(dim):
                for itrj in range(dim):
                    var = barevar(str(self.symb[itri][itrj]), self.vartype, self.varprefix)
                    if str(self.symb[itri][itrj])[0] == '-':
                        Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj])[1:])
                        self.isymb[itri][itrj] = -IndexedBase(self.varprefix,real = isreal)[i,j,k,Idxsymb]
                    elif self.symb[itri][itrj] == 0:
                        self.isymb[itri][itrj] = 0
                    else:
                        Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj]))
                        self.isymb[itri][itrj] = IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                    if addtolist:
                        if idx.reducebysymmetry2(itri,itrj,self.sym):
                            stvar.gridvars.append(var)
                        if self.vartype not in stvar.vartypes:
                            stvar.vartypes.append(self.vartype)

        #if resetsym == True:
         #   self.symb = ixp.declarerank2(str(symb), 'nosym', DIM=dim)
            
        self.symb = np.array(self.symb)
        #Should I do the same for expr and isymb?  Test this later on.
        
    def symb2expr(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                if self.symb[itri][itrj] != 0:
                    if idx.reducebysymmetry2(itri,itrj,self.sym):
                        expr += "        amrex::Real " + str(self.symb[itri][itrj]) + " = " + amrexcode(self.expr[itri][itrj],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def symb2isymb(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                if self.symb[itri][itrj] != 0:
                    if idx.reducebysymmetry2(itri,itrj,self.sym):
                        expr += "        amrex::Real " + str(self.symb[itri][itrj]) + " = " + amrexcode(self.isymb[itri][itrj],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def setisymb(self, prefix = ""):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                if self.isymb[itri][itrj] != 0:
                    if prefix == "":
                        if idx.reducebysymmetry2(itri,itrj,self.sym):
                            expr += "        " + amrexcode(self.isymb[itri][itrj],stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj],stvar.gridvars)+";\n"
                    else:
                        Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj]))
                        tempisymb = IndexedBase(prefix)[i,j,k,Idxsymb]
                        if idx.reducebysymmetry2(itri,itrj,self.sym):
                            expr += "        " + amrexcode(tempisymb, stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj],stvar.gridvars)+ ";\n"
        expr += '\n'
        return expr
    
    def diff(self, dstring, Accuracy = 2, dim = 3, makesparse = True):
        dsymbrootname, dsymbrank = fd.dsymb_name_rank(self, dstring)
        if 'DD' in dstring:
            if self.sym != 'nosym':
                stvarout = stvarrank4(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = self.sym + 'sym23', dim = dim)
            else:
                stvarout = stvarrank4(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = 'sym23', dim = dim)
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for l in range(dim):
                            stvarout.isymb[i][j][k][l] = fd.disymb(self.isymb[i][j], dstring + str(k) + str(l), Accuracy)
                            if (makesparse == True) and (self.symb[i][j] == 0):
                                stvarout.symb[i][j][k][l] = 0
                                stvarout.isymb[i][j][k][l] = 0
        elif 'dKODFull' in dstring:
            stvarout = stvarrank2(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = self.sym, dim = dim)
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        stvarout.isymb[i][j] += fd.disymb(self.isymb[i][j], 'dKOD'+str(k), Accuracy)
                        if (makesparse == True) and (self.symb[i][j] == 0):
                            stvarout.symb[i][j] = 0
                            stvarout.isymb[i][j] = 0
        else:
            stvarout = stvarrank3(ssp.unparse_name_rank(dsymbrootname, dsymbrank), sym = self.sym, dim = dim)
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        stvarout.isymb[i][j][k] = fd.disymb(self.isymb[i][j], dstring + str(k), Accuracy)
                        if (makesparse == True) and (self.symb[i][j] == 0):
                            stvarout.symb[i][j][k] = 0
                            stvarout.isymb[i][j][k] = 0
        stvarout.expr = stvarout.isymb
        return stvarout
    
class stvarrank3:
    def __init__(self, symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", isreal = False, dim = 3, sym = 'nosym', resetsym = True):
        self.symb = idx.declarerank3(str(symb), symmetry_option = sym, DIM=dim)
        self.expr = idx.zerorank3(DIM=dim)
        self.isymb = idx.zerorank3(DIM=dim)
        self.vartype = vartype
        self.rootname, self.rank = ssp.parse_name_rank(symb)
        if len(self.rank) != 3:
            print("The rank must 3.")
        self.dim = dim
        self.sym = sym
        
        if(varprefix == ""):
            self.varprefix = vartype
            
        else:
            self.varprefix = varprefix
            
        if gridvar == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        var = barevar(str(self.symb[itri][itrj][itrk]), self.vartype, self.varprefix)
                        if str(self.symb[itri][itrj][itrk])[0] == '-':
                            Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk])[1:])
                            self.isymb[itri][itrj][itrk] = -IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                        elif self.symb[itri][itrj][itrk] == 0:
                            self.isymb[itri][itrj][itrk] = 0
                        else:
                            Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk]))
                            self.isymb[itri][itrj][itrk] = IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                        if addtolist:
                            if idx.reducebysymmetry3(itri, itrj, itrk, self.sym):
                                stvar.gridvars.append(var)
                            if self.vartype not in stvar.vartypes:
                                stvar.vartypes.append(self.vartype)
                                

        #if resetsym == True:
         #   self.symb = ixp.declarerank2(str(symb), 'nosym', DIM=dim)
            
        self.symb = np.array(self.symb)
        #Should I do the same for expr and isymb?  Test this later on.
        
    def symb2expr(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    if self.symb[itri][itrj][itrk] != 0:
                        if idx.reducebysymmetry3(itri, itrj, itrk, self.sym):
                            expr += "        amrex::Real " + str(self.symb[itri][itrj][itrk]) + " = " + amrexcode(self.expr[itri][itrj][itrk],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def symb2isymb(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    if self.symb[itri][itrj][itrk] != 0:
                        if idx.reducebysymmetry3(itri, itrj, itrk, self.sym):
                            expr += "        amrex::Real " + str(self.symb[itri][itrj][itrk]) + " = " + amrexcode(self.isymb[itri][itrj][itrk],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def setisymb(self, prefix = ""):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    if self.isymb[itri][itrj][itrk] != 0:
                        if prefix == "":
                            if idx.reducebysymmetry3(itri, itrj, itrk, self.sym):
                                expr += "        " + amrexcode(self.isymb[itri][itrj][itrk],stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj][itrk],stvar.gridvars)+";\n"
                        else:
                            Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk]))
                            tempisymb = IndexedBase(prefix)[i,j,k,Idxsymb]
                            if idx.reducebysymmetry3(itri, itrj, itrk, self.sym):
                                expr += "        " + amrexcode(tempisymb, stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj][itrk],stvar.gridvars)+ ";\n"
        expr += '\n'
        return expr
    
class stvarrank4:
    def __init__(self, symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", isreal = False, dim = 3, sym = 'nosym', resetsym = True):
        self.symb = idx.declarerank4(str(symb), symmetry_option = sym, DIM=dim)
        self.expr = idx.zerorank4(DIM=dim)
        self.isymb = idx.zerorank4(DIM=dim)
        self.vartype = vartype
        self.rootname, self.rank = ssp.parse_name_rank(symb)
        if len(self.rank) != 4:
            print("The rank must 4.")
        self.dim = dim
        self.sym = sym
        
        if(varprefix == ""):
            self.varprefix = vartype
            
        else:
            self.varprefix = varprefix
            
        if gridvar == True:
            for itri in range(dim):
                for itrj in range(dim):
                    for itrk in range(dim):
                        for itrl in range(dim):
                            var = barevar(str(self.symb[itri][itrj][itrk][itrl]), self.vartype, self.varprefix)
                            if str(self.symb[itri][itrj][itrk])[0] == '-':
                                Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk][itrl])[1:])
                                self.isymb[itri][itrj][itrk][itrl] = -IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                            elif self.symb[itri][itrj][itrk][itrl] == 0:
                                self.isymb[itri][itrj][itrk][itrl] = 0
                            else:
                                Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk][itrl]))
                                self.isymb[itri][itrj][itrk][itrl] = IndexedBase(self.varprefix, real = isreal)[i,j,k,Idxsymb]
                            if addtolist:
                                if idx.reducebysymmetry4(itri, itrj, itrk, itrl, self.sym):
                                    stvar.gridvars.append(var)
                                if self.vartype not in stvar.vartypes:
                                    stvar.vartypes.append(self.vartype)

        #if resetsym == True:
         #   self.symb = ixp.declarerank2(str(symb), 'nosym', DIM=dim)
            
        self.symb = np.array(self.symb)
        #Should I do the same for expr and isymb?  Test this later on.
        
    def symb2expr(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    for itrl in range(self.dim):
                        if self.symb[itri][itrj][itrk][itrl] != 0:
                            if idx.reducebysymmetry4(itri, itrj, itrk, itrl, self.sym):
                                expr += "        amrex::Real " + str(self.symb[itri][itrj][itrk][itrl]) + " = " + amrexcode(self.expr[itri][itrj][itrk][itrl],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def symb2isymb(self):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    for itrl in range(self.dim):
                        if self.symb[itri][itrj][itrk][itrl] != 0:
                            if idx.reducebysymmetry4(itri, itrj, itrk, itrl, self.sym):
                                expr += "        amrex::Real " + str(self.symb[itri][itrj][itrk][itrl]) + " = " + amrexcode(self.isymb[itri][itrj][itrk][itrl],stvar.gridvars)+";\n"
        expr += '\n'
        return expr
    
    def setisymb(self, prefix = ""):
        expr = ""
        for itri in range(self.dim):
            for itrj in range(self.dim):
                for itrk in range(self.dim):
                    for itrl in range(self.dim):
                        if self.isymb[itri][itrj][itrk][itrl] != 0:
                            if prefix == "":
                                if idx.reducebysymmetry4(itri, itrj, itrk, itrl, self.sym):
                                    expr += "        " + amrexcode(self.isymb[itri][itrj][itrk][itrl],stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj][itrk][itrl],stvar.gridvars)+";\n"
                            else:
                                Idxsymb = symbols(self.vartype+str(self.symb[itri][itrj][itrk][itrl]))
                                tempisymb = IndexedBase(prefix)[i,j,k,Idxsymb]
                                if idx.reducebysymmetry4(itri, itrj, itrk, itrl, self.sym):
                                    expr += "        " + amrexcode(tempisymb, stvar.gridvars) + " = " + amrexcode(self.expr[itri][itrj][itrk][itrl],stvar.gridvars)+ ";\n"
        expr += '\n'
        return expr
    
def makestvar(symb, expr=0, isymb=0, gridvar = False, addtolist = True, vartype = "", varprefix = "", dim = 3, sym = 'nosym', resetsym = True):
    if symb.find('_') == -1:
        return stvar(symb, expr, isymb, gridvar, addtolist, vartype, varprefix)
    else:
        index_type_string = symb[symb.find('_')+1:]
        
        if len(index_type_string) == 1:
            return stvarrank1(symb, expr, isymb, gridvar, addtolist, vartype, varprefix, dim)
        elif len(index_type_string) == 2:
            return stvarrank2(symb, expr, isymb, gridvar, addtolist, vartype, varprefix, dim, sym, resetsym)
        elif len(index_type_string) == 3:
            return stvarrank2(symb, expr, isymb, gridvar, addtolist, vartype, varprefix, dim, sym, resetsym)
        elif len(index_type_string) == 4:
            return stvarrank2(symb, expr, isymb, gridvar, addtolist, vartype, varprefix, dim, sym, resetsym)
        else:
            print("Only tensors up to rank 4 are supported")
    
def partitionbytype():
    partitionedlist = []
    for i in range(len(stvar.vartypes)):
        partitionedlist.append([x for x in stvar.gridvars if x.vartype == stvar.vartypes[i]])
        
    return partitionedlist
    
def tagcondition(stvar1, conditional, stvar2):
    String = stvar('String')
    String.symb = stvar1.expr
    String.expr = stvar2.expr
    
    return "        return " + amrexcode(stvar1.expr,stvar.gridvars) + ' ' + conditional + ' ' + amrexcode(stvar2.expr,stvar.gridvars)+";\n\n"