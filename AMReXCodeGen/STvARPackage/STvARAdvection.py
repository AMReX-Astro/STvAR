from STvAR_v3 import *
from STvARFiniteDifference import *


def AdvectiveD(tensor, V, upwindlevel= '1'):
    Dim = V.dim
    if (np.ndim(tensor.symb) == 0):
        Outtensor = stvar(ssp.unparse_name_rank('AdvD' + V.rootname + tensor.rootname, tensor.rank))
        for i in range(Dim):
            Outtensor.expr += V.symb[i]*((V.symb[i]>0)*dsymb(tensor.symb,'dup' + upwindlevel + 'D'+ str(i))+(V.symb[i]<0)*dsymb(tensor.symb,'ddn' + upwindlevel + 'D'+str(i)))
    elif (np.ndim(tensor.symb) == 1):
        Outtensor = stvarrank1(ssp.unparse_name_rank('AdvD' + V.rootname + tensor.rootname, tensor.rank))
        for i in range(Dim):
            for j in range(Dim):
                Outtensor.expr[i] += V.symb[j]*((V.symb[j]>0)*dsymb(tensor.symb[i],'dup' + upwindlevel + 'D'+ str(j))+(V.symb[j]< 0)*dsymb(tensor.symb[i],'ddn' + upwindlevel + 'D'+str(j)))
        
    elif (np.ndim(tensor.symb) == 2):
        Outtensor = stvarrank2(ssp.unparse_name_rank('AdvD' + V.rootname + tensor.rootname, tensor.rank), sym = tensor.sym)
        for i in range(Dim):
            for j in range(Dim):
                for k in range(Dim):
                    Outtensor.expr[i][j] += V.symb[k]*((V.symb[k]>0)*dsymb(tensor.symb[i][j],'dup' + upwindlevel + 'D'+ str(k))+(V.symb[k] < 0)*dsymb(tensor.symb[i][j],'ddn' + upwindlevel + 'D'+str(k)))
        
    
    return Outtensor
    
    