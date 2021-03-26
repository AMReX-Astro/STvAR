import sympy as sp
import STvARIndexing as idx
import STvARSymbolParsing as ssp
from sympy import symbols, IndexedBase, Indexed, Idx

def FiniteDiff(dstring, Accuracy = 2):
    StencilSize = Accuracy + 1
    Shift = 0
    if "dup" in dstring:
        if "dupD" in dstring:
            Shift = 1
        else:
            Shift = int(dstring[3])
    elif "ddn" in dstring:
        if "ddnD" in dstring:
            Shift = -1
        else:
            Shift = -int(dstring[3])
        
    M = sp.zeros(StencilSize,StencilSize)
    MinIdx = - sp.Rational(Accuracy,2)+Shift 
    for i in range(StencilSize):
        for j in range(StencilSize):
            M[(i,j)] = (MinIdx + j)**(i)
    Minv = sp.zeros(StencilSize,StencilSize)
    Minv = M**(-1)
    matrixcol = 1
    difftype = "First"
    if "DD" in dstring:
        if dstring[len(dstring)-1] == dstring[len(dstring)-2]:
            difftype = "Second"
            matrixcol = 2
        else:
            difftype = "Mixed"
    elif "dKOD" in dstring:
        difftype = "KreissOliger"
        matrixcol = StencilSize - 1
    else:
        pass
    
    FDCoeffs = []
    FDStencil = []
    if difftype != "Mixed":
        for i in range(StencilSize):
            Idx4 = [0,0,0,0]
            FDCoeff = sp.factorial(matrixcol)*Minv[(i,matrixcol)]
            if FDCoeff != 0:
                FDCoeffs.append(FDCoeff)
                if difftype == "KreissOliger":
                    FDCoeffs[i] *= (-1)**(sp.Rational((StencilSize+1),2))/2**matrixcol
                gridpos = i + int(MinIdx)
                if gridpos != 0:
                    direction = int(dstring[len(dstring)-1])
                    Idx4[direction] = gridpos
                FDStencil.append(Idx4)
    else:
        for i in range(StencilSize):
            for j in range(StencilSize):
                Idx4 = [0,0,0,0]
                FDCoeff = (sp.factorial(matrixcol)*Minv[(i,matrixcol)]) * \
                          (sp.factorial(matrixcol)*Minv[(j,matrixcol)])
                
                if FDCoeff != 0:
                    FDCoeffs.append(FDCoeff)
                    gridpos1 = i + int(MinIdx)
                    gridpos2 = j + int(MinIdx)
                    direction1 = int(dstring[len(dstring)-2])
                    direction2 = int(dstring[len(dstring)-1])
                    Idx4[direction1] = gridpos1
                    Idx4[direction2] = gridpos2
                    FDStencil.append(Idx4)
                    
    return FDCoeffs, FDStencil

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

def disymb(E, dstring, Accuracy = 2):
    DX = idx.zerorank1(DIM=5)
    DX[0], DX[1], DX[2], DX[3] = symbols('dx[0] dx[1] dx[2] dx[3]')
    delta = DX[int(dstring[len(dstring)-1])]
    if 'DD' in dstring:
        delta *= DX[int(dstring[len(dstring)-2])]
    fdcoeffs, fdstencl = FiniteDiff(dstring, Accuracy)
    shiftE = 0
    if E != 0:
        for i in range(len(fdcoeffs)):
            shiftE += fdcoeffs[i]*shift(E,fdstencl[i])
        shiftE *= delta**(-1) 
    return  shiftE


def dsymb(symb, dstring):
    Shift = 0
    direction = dstring[len(dstring) - 1]
    if str(symb)[0] == '-':
        symbout = symbols(str(symb)[1:])
    else:
        symbout = symb
        
    symbname, rankval, componentval= ssp.parse_name_rank_component(str(symbout))
    
    if 'dup' in dstring:
        if 'dupD' in dstring or 'dup1D' in dstring:
            prestring = 'dupD'
        else:
            prestring = 'dup' + dstring[3] + 'D'
    elif 'ddn' in dstring:
        if 'ddnD' in dstring or 'ddn1D' in dstring:
            prestring = 'ddnD'
        else:
            prestring = 'ddn' + dstring[3] + 'D'
    
    elif 'dKOD' in dstring:
        prestring = 'dKOD'
        
    else:
        prestring = 'dD'
    
    if 'DD' in dstring:
        prestring += 'D'
        if direction < dstring[len(dstring) - 2]:
            direction = direction + dstring[len(dstring) - 2]
        else:
            direction = dstring[len(dstring) - 2] + direction 
        rankval += 'LL'
        componentval += direction
        prestring = ssp.unparse_name_rank_component(prestring + symbname, rankval, componentval)
    elif 'dKODFull' in dstring:
        prestring = ssp.unparse_name_rank_component(prestring + symbname, rankval, componentval)
    else:
        rankval += 'L'
        componentval += direction
        prestring = ssp.unparse_name_rank_component(prestring+symbname, rankval, componentval)
    
    if str(symb)[0] == '-':
        return -symbols(prestring)
    elif symb == 0:
        return 0
    else:
        return symbols(prestring)
    
def dsymb_name_rank(symb, dstring):
    dstring += '0'
    if 'DD' in dstring:
        dstring += '0'
    rootname, rank = symb.rootname, symb.rank
    dsymb_parsed = ssp.parse_name_rank_component(str(dsymb(rootname, dstring)))
    return dsymb_parsed[0], rank+dsymb_parsed[1]
