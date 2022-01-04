from STvAR_v3 import *
from STvARFiniteDifference import *
import STvARIndexing as sidx
import STvARSymbolParsing as ssp

def inversemetric(metric):
    Dim = metric.dim

    if metric.rank == 'LL':
        inverserank = 'UU'
    elif metric.rank == 'UU':
        inverserank = 'LL'
    else:
        print("Rank type nsupported")
        
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, inverserank), sym = metric.sym, dim = Dim)
    
    if(Dim == 3):
        inversemetric.expr = sidx.symm_matrix_inverter3x3(metric.symb)
    elif (Dim == 4):
        inversemetric.expr = sidx.symm_matrix_inverter4x4(metric.symb)
    else:
        for i in range(Dim):
            for j in range(Dim):
                inversemetric.expr[i][j] += sp.simplify(sp.Matrix(metric.symb).inv()[i,j])
    return inversemetric

def detmetric(metric):
    Dim = metric.dim
    
    detmetric = stvar("det" + metric.rootname)
    
    if(Dim == 3):
        detmetric.expr = sidx.symm_matrix_det3x3(metric.symb)
    elif (Dim == 4):
        detmetric.expr = sidx.symm_matrix_det4x4(metric.symb)
    else:
        for i in range(Dim):
            for j in range(Dim):
                inversemetric.expr[i][j] += sp.simplify(sp.Matrix(metric.symb).det())
    return detmetric
    
def ChristoffelLLL(metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    rank = 'LLL'
    if metric.rank != 'LL':
        print("Rank must be LL type")
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, rank), sym = 'sym12', dim=Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                ChristoffelLLL.expr[i][j][k] += 1/2*(dsymb(metric.symb[i][j], dstring+str(k)) + dsymb(metric.symb[i][k], dstring+str(j))- dsymb(metric.symb[j][k], dstring+str(i)))
    
    return ChristoffelLLL

def ChristoffelULL(metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    rank = 'ULL'
    if metric.rank != 'LL':
        print("Rank must be LL type")
        
    inverserank = 'UU'
        
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, inverserank), sym = metric.sym, dim = Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, rank), sym = 'sym12', dim=Dim)
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'LLL'), sym = 'sym12', dim = Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    ChristoffelULL.expr[i][j][k] += inversemetric.symb[i][l]*ChristoffelLLL.symb[l][j][k]
                    
    return ChristoffelULL
                                 
def ConformalChristoffelULL(chi, metric, dstring = 'dD'):
    Dim = metric.dim
    metricrootname = metric.rootname
    chirootname = chi.rootname
    rank = 'ULL'
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metricrootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + metricrootname, rank), sym = 'sym12', dim=Dim)
    
    ConformalChristoffelULL = stvarrank3(ssp.unparse_name_rank('ConformalChristoffel' + chirootname + metricrootname, rank), sym = 'sym12', dim=Dim)
        
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                ConformalChristoffelULL.expr[i][j][k] += ChristoffelULL.symb[i][j][k] - 1/2*((sp.eye(3)[i,j]/chi.symb)*dsymb(chi.symb, dstring+str(k))+(sp.eye(3)[i,k]/chi.symb)*dsymb(chi.symb, dstring + str(j)))
                for l in range(Dim):
                    ConformalChristoffelULL.expr[i][j][k] += 1/2*metric.symb[j][k]*inversemetric.symb[i][l]*dsymb(chi.symb, dstring+ str(l))/chi.symb
                    
    return ConformalChristoffelULL

def ConformalChristoffelLLL(chi, metric, dstring = 'dD'):
    Dim = metric.dim
    metricrootname = metric.rootname
    chirootname = chi.rootname
    rank = 'LLL'
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    ConformalChristoffelULL = stvarrank3(ssp.unparse_name_rank('ConformalChristoffel'+ chirootname + metricrootname, 'ULL'), sym = 'sym12', dim=Dim)
    ConformalChristoffelLLL = stvarrank3(ssp.unparse_name_rank('ConformalChristoffel'+ chirootname + metricrootname, 'LLL'), sym = 'sym12', dim=Dim)
    
        
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    ConformalChristoffelLLL.expr[i][j][k] += (metric.symb[i][l]/(chi.symb))*ConformalChristoffelULL.symb[l][j][k]
                    
    return ConformalChristoffelLLL

def partialinversemetricUUL(metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    dsymb_name, dsymb_rank = dsymb_name_rank(inversemetric,dstring)
    partialinversemetric = stvarrank3(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym01')
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    partialinversemetric.expr[i][j][k] += -ChristoffelULL.symb[i][k][l]*inversemetric.symb[l][j]-ChristoffelULL.symb[j][k][l]*inversemetric.symb[l][i]
    
    return partialinversemetric

def ChristoffelDU(metric):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    
    ChristoffelDU = stvarrank1(ssp.unparse_name_rank('ChristoffelD' + rootname, 'U'), dim = Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                ChristoffelDU.expr[i] += inversemetric.symb[j][k]*ChristoffelULL.symb[i][j][k]
                
    return ChristoffelDU

def partialChristoffelLLLL(metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname

    if metric.rank != 'LL':
        print("Rank must be LL type")
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'LLL'), sym = 'sym12', dim=Dim)
    dsymb_name, dsymb_rank = dsymb_name_rank(ChristoffelLLL, dstring)
    partialChristoffelLLLL = stvarrank4(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym12')
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    partialChristoffelLLLL.expr[i][j][k][l] += 1/2*(dsymb(metric.symb[i][j], dstring +'D'+str(k)+str(l))+dsymb(metric.symb[i][k], dstring +'D'+str(j)+str(l))-dsymb(metric.symb[j][k], dstring +'D'+str(i)+str(l)))
                    
    return partialChristoffelLLLL

def partialChristoffelULLL(metric, dstring= 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'LLL'), sym = 'sym12', dim=Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    
    dsymb_name, dsymb_rank = dsymb_name_rank(inversemetric,dstring)
    partialinversemetric = stvarrank3(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym01')
    
    dsymb_name, dsymb_rank = dsymb_name_rank(ChristoffelLLL, dstring)
    partialChristoffelLLLL = stvarrank4(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym12')
    
    dsymb_name, dsymb_rank = dsymb_name_rank(ChristoffelULL, dstring)
    partialChristoffelULLL = stvarrank4(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym12')
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    for m in range(Dim):
                        partialChristoffelULLL.expr[i][j][k][l] += partialinversemetric.symb[i][m][l]*ChristoffelLLL.symb[m][j][k]
                        partialChristoffelULLL.expr[i][j][k][l] += inversemetric.symb[i][m]*partialChristoffelLLLL.symb[m][j][k][l]
    
    return partialChristoffelULLL

def RiemannTensorLLLL(metric, dstring = 'dDD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
        
    RiemannTensorLLLL = stvarrank4(ssp.unparse_name_rank("Riemann" + metric.rootname, 'LLLL'), sym = 'asym01asym23', dim = Dim)
    
    for i in range(Dim):
        for k in range(Dim):
            for l in range(Dim):
                for m in range(Dim):
                    for m in range(Dim):
                        RiemannTensorLLLL.expr[i][j][k][l] += metric.symb[i][m]*RiemannTensorULLL.symb[m][j][k][l]
                        
    return RiemannTensorLLLL

def RiemannTensorULLL(metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
        
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    
    dsymb_name, dsymb_rank = dsymb_name_rank(ChristoffelULL, dstring)
    partialChristoffelULLL = stvarrank4(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym12')
    
    RiemannTensorULLL = stvarrank4(ssp.unparse_name_rank("Riemann" + metric.rootname, 'ULLL'), sym = 'asym23', dim = Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    RiemannTensorULLL.expr[i][j][k][l] += partialChristoffelULLL.symb[i][j][l][k]-partialChristoffelULLL.symb[i][j][k][l]
                    for m in range(Dim):
                        RiemannTensorULLL.expr[i][j][k][l] += ChristoffelULL.symb[i][m][k]*ChristoffelULL.symb[m][j][l]-ChristoffelULL.symb[i][m][l]*ChristoffelULL.symb[m][j][k]
    
    return RiemannTensorULLL
    

def RicciTensorLL(metric, dstring = 'dDD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'LLL'), sym = 'sym12', dim=Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    RicciTensorLL = stvarrank2(ssp.unparse_name_rank("Ricci" + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                for l in range(Dim):
                    RicciTensorLL.expr[i][j] += 1/2*inversemetric.symb[k][l]*(dsymb(metric.symb[k][j], dstring + str(i) + str(l)))
                    RicciTensorLL.expr[i][j] += 1/2*inversemetric.symb[k][l]*(dsymb(metric.symb[i][l], dstring + str(k) + str(j)))
                    RicciTensorLL.expr[i][j] += -1/2*inversemetric.symb[k][l]*(dsymb(metric.symb[k][l], dstring + str(i) + str(j)))
                    RicciTensorLL.expr[i][j] += -1/2*inversemetric.symb[k][l]*(dsymb(metric.symb[i][j], dstring + str(k) + str(l)))
                    
                    for m in range(Dim):
                        RicciTensorLL.expr[i][j] += inversemetric.symb[k][l]*(ChristoffelULL.symb[m][i][l]*ChristoffelLLL.symb[m][k][j]-ChristoffelULL.symb[m][i][j]*ChristoffelLLL.symb[m][k][l])
                        
    return RicciTensorLL

def RicciScalar(metric):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
        
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)    
    RicciTensorLL = stvarrank2(ssp.unparse_name_rank("Ricci" + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    RicciScalar = stvar("RicciScalar" + metric.rootname)
    
    for i in range(Dim):
        for j in range(Dim):
            RicciScalar.expr += inversemetric.symb[i][j]*RicciTensorLL.symb[i][j]
        
    return RicciScalar

def Z4cRicciTildeTensorLL(metric, GamU, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    ChristoffelLLL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'LLL'), sym = 'sym12', dim=Dim)
    ChristoffelULL = stvarrank3(ssp.unparse_name_rank('Christoffel' + rootname, 'ULL'), sym = 'sym12', dim=Dim)
    ChristoffelDU = stvarrank1(ssp.unparse_name_rank('ChristoffelD' + rootname, 'U'), dim = Dim)
    Z4cRicciTildeTensorLL = stvarrank2(ssp.unparse_name_rank("Z4cRicciTilde" + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    
    
    for i in range(Dim):
        for j in range(Dim):
            for k in range(Dim):
                Z4cRicciTildeTensorLL.expr[i][j] += 1/2*(metric.symb[k][i]*dsymb(GamU.symb[k], dstring + str(j))+metric.symb[k][j]*dsymb(GamU.symb[k],dstring + str(i)))
                Z4cRicciTildeTensorLL.expr[i][j] += 1/2*ChristoffelDU.symb[k]*(ChristoffelLLL.symb[i][j][k]+ChristoffelLLL.symb[j][i][k])
                for l in range(Dim):
                    Z4cRicciTildeTensorLL.expr[i][j] += -1/2*inversemetric.symb[k][l]*dsymb(metric.symb[i][j], dstring + 'D'+str(k)+str(l))
                    for m in range(Dim):
                        Z4cRicciTildeTensorLL.expr[i][j] += inversemetric.symb[l][m]*(ChristoffelULL.symb[k][l][i]*ChristoffelLLL.symb[j][k][m]+ChristoffelULL.symb[k][l][j]*ChristoffelLLL.symb[i][k][m])
                        Z4cRicciTildeTensorLL.expr[i][j] += inversemetric.symb[l][m]*ChristoffelULL.symb[k][i][m]*ChristoffelLLL.symb[k][l][j]
                    
    return Z4cRicciTildeTensorLL

def CovariantD(tensor, ChristoffelULL, dstring = 'dD'):
    Dim = ChristoffelULL.dim
    dsymbrootname, dsymbrank = dsymb_name_rank(tensor, dstring)
    if (np.ndim(tensor.symb) == 0):
        if 'DD' in dstring:
            first_order_dstring = "".join(list(dstring)[0:-1])
            Outtensor = stvarrank2(ssp.unparse_name_rank('Cov' + dsymbrootname, dsymbrank), sym = 'sym01', dim = Dim)
            
            for i in range(Dim):
                for j in range(Dim):
                    Outtensor.expr[i][j] += dsymb(tensor.symb, dstring + str(i)+str(j))
                    for k in range(Dim):
                        Outtensor.expr[i][j] += -ChristoffelULL.symb[k][i][j]*dsymb(tensor.symb, first_order_dstring + str(k))
        else:
            Outtensor = stvarrank1(ssp.unparse_name_rank('Cov' + dsymbrootname, dsymbrank), dim = Dim)
            for i in range(Dim):
                Outtensor.expr[i] += dsymb(tensor.symb, dstring + str(i))
                
    elif (np.ndim(tensor.symb) == 1):
        if 'DD' in dstring:
            first_order_dstring = "".join(list(dstring)[0:-1])
            Outtensor = stvarrank3(ssp.unparse_name_rank('Cov' + dsymbrootname, dsymbrank), sym = 'sym12', dim = Dim)
            dsymb_name, dsymb_rank = dsymb_name_rank(ChristoffelULL, first_order_dstring)
            partialChristoffelULLL = stvarrank4(ssp.unparse_name_rank(dsymb_name,dsymb_rank), sym = 'sym12')
            for k in range(Dim):
                for i in range(Dim):
                    for j in range(Dim):
                        Outtensor.expr[k][i][j] += dsymb(tensor.symb[k], dstring + str(i) + str(j))
                        for l in range(Dim):
                            if tensor.rank == 'U':
                                Outtensor.expr[k][i][j] += partialChristoffelULLL.symb[k][j][l][i]*tensor.symb[l] + ChristoffelULL.symb[k][j][l]*dsymb(tensor.symb[l],"".join(list(dstring)[0:-1]) + str(i))
                                Outtensor.expr[k][i][j] += -ChristoffelULL.symb[l][i][j]*dsymb(tensor.symb[k],"".join(list(dstring)[0:-1]) + str(l))
                                Outtensor.expr[k][i][j] += ChristoffelULL.symb[k][i][l]*dsymb(tensor.symb[l],"".join(list(dstring)[0:-1]) + str(j))
                                for m in range(Dim):
                                    Outtensor.expr[k][i][j] += -ChristoffelULL.symb[m][i][j]*ChristoffelULL.symb[k][m][l]*tensor.symb[l]
                                    Outtensor.expr[k][i][j] += ChristoffelULL.symb[k][i][m]*ChristoffelULL.symb[m][j][l]*tensor.symb[l]
                            elif tensor.rank == 'L':
                                Outtensor.expr[k][i][j] += -partialChristoffelULLL.symb[l][j][k][i]*tensor.symb[l] - ChristoffelULL.symb[l][j][k]*dsymb(tensor.symb[l],"".join(list(dstring)[0:-1]) + str(i))
                                Outtensor.expr[k][i][j] += -ChristoffelULL.symb[l][i][j]*dsymb(tensor.symb[k],"".join(list(dstring)[0:-1]) + str(l))
                                Outtensor.expr[k][i][j] += -ChristoffelULL.symb[l][i][k]*dsymb(tensor.symb[l],"".join(list(dstring)[0:-1]) + str(j))
                                for m in range(Dim):
                                    Outtensor.expr[k][i][j] += ChristoffelULL.symb[m][i][j]*ChristoffelULL.symb[l][m][k]*tensor.symb[l]
                                    Outtensor.expr[k][i][j] += ChristoffelULL.symb[m][i][k]*ChristoffelULL.symb[l][j][m]*tensor.symb[l]
                                
                
        else:
            Outtensor = stvarrank2(ssp.unparse_name_rank('Cov' + dsymbrootname, dsymbrank), sym = 'nosym', dim = Dim)
            for i in range(Dim):
                for j in range(Dim):
                    Outtensor.expr[i][j] += dsymb(tensor.symb[i], dstring + str(j))
                    for k in range(Dim):
                        if tensor.rank == 'U':
                            Outtensor.expr[i][j] += ChristoffelULL.symb[i][j][k]*tensor.symb[k]
                        elif tensor.rank == 'L':
                            Outtensor.expr[i][j] += -ChristoffelULL.symb[k][i][j]*tensor.symb[k]
                
    return Outtensor

def Z4cRicciTensorFromScalar(chi, metric, dstring = 'dD'):
    Dim = metric.dim
    dsymbrootname_order1, dsymbrank_order1 = dsymb_name_rank(chi, dstring)
    dsymbrootname_order2, dsymbrank_order2 = dsymb_name_rank(chi, dstring+'D')
    
    
    CovDchi = stvarrank1(ssp.unparse_name_rank('Cov' + dsymbrootname_order1, dsymbrank_order1), dim = Dim)
    CovDDchi = stvarrank2(ssp.unparse_name_rank('Cov' + dsymbrootname_order2, dsymbrank_order2), sym = 'sym01', dim = Dim)
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    
    Z4cRicciTensorFromScalar = stvarrank2(ssp.unparse_name_rank("Z4cRicciFrom" + chi.rootname + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    
    for i in range(Dim):
        for j in range(Dim):
            Z4cRicciTensorFromScalar.expr[i][j] += 1/(2*chi.symb)*CovDDchi.symb[i][j] - 1/(4*chi.symb**2)*CovDchi.symb[i]*CovDchi.symb[j]
            for k in range(Dim):
                for l in range(Dim):
                    Z4cRicciTensorFromScalar.expr[i][j] += metric.symb[i][j]*(1/(2*chi.symb)*inversemetric.symb[k][l]*CovDDchi.symb[k][l]-3/(4*chi.symb**2)*inversemetric.symb[k][l]*CovDchi.symb[k]*CovDchi.symb[l])
                    
    return Z4cRicciTensorFromScalar
    
def Z4cRicciTensorLL(chi, metric, dstring = 'dD'):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    Z4cRicciTildeTensorLL = stvarrank2(ssp.unparse_name_rank("Z4cRicciTilde" + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    Z4cRicciTensorFromScalar = stvarrank2(ssp.unparse_name_rank("Z4cRicciFrom" + chi.rootname + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    Z4cRicciTensor = stvarrank2(ssp.unparse_name_rank("Z4cRicci" + chi.rootname + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
                                         
    
    Z4cRicciTensorLL.expr = Z4cRicciTildeTensorLL.symb + Z4cRicciTensorFromScalar.symb
    
    return Z4cRicciTensorLL
                        
def Z4cRicciScalar(chi, metric):
    Dim = metric.dim
    rootname = metric.rootname
    
    if metric.rank != 'LL':
        print("Rank must be LL type")
    
    inversemetric = stvarrank2(ssp.unparse_name_rank("inv" + metric.rootname, 'UU'), sym = metric.sym, dim = Dim)
    Z4cRicciTensor = stvarrank2(ssp.unparse_name_rank("Z4cRicciFrom" + chi.rootname + metric.rootname, 'LL'), sym = 'sym01', dim = Dim)
    Z4cRicciScalar = stvar("Z4cRicciScalar" + chi.rootname + metric.rootname)
    
    for i in range(Dim):
        for j in range(Dim):
            Z4cRicciScalar.expr += chi.symb*inversemetric.symb[i][j]*Z4cRicciTensor.symb[i][j]
     
    return Z4cRicciScalar
    
    
    
    