import sympy as sp
import sys

def myEij(i,j):
    return sp.sign(sp.Eijk(i,j))

def myEijk(i,j,k):
    return sp.sign(sp.Eijk(i,j,k))

def myEijkl(i,j,k,l):
    return sp.sign(sp.Eijk(i,j,k,l))

def reducebysymmetry2(i, j, sym):
    istrue = False
    if sym == 'sym01':
        istrue = (i <= j)
    elif sym == 'asym01':
        istrue = (i < j)
    elif sym == 'diag01':
        istrue = (i == j)
    elif sym == 'nosym':
        istrue = True
    return istrue

def reducebysymmetry3(i,j,k,sym):
    istrue = True
    if 'asym01' in sym:
        istrue *= (i < j)
    elif 'sym01' in sym:
        istrue *= (i <= j)
    elif 'diag01' in sym:
        istrue *= (i == j)
    if 'asym12' in sym:
        istrue *= (j < k)
    elif 'sym12' in sym:
        istrue = (j <= k)
    elif 'diag12' in sym:
        istrue *= (j == k)
    if 'asym02' in sym:
        istrue *= (i < k)
    elif 'sym02' in sym:
        istrue *= (i <= k)
    elif 'diag02' in sym:
        istrue *= (i == k)
    elif sym == 'nosym':
        istrue = True
    return istrue

def reducebysymmetry4(i, j, k, l, sym):
    istrue = True
    if 'asym01' in sym:
        istrue *= (i < j)
    elif 'sym01' in sym:
        istrue *= (i <= j)
    elif 'diag01' in sym:
        istrue *= (i == j)
    if 'asym02' in sym:
        istrue *= (i < k)
    elif 'sym02' in sym:
        istrue *= (i <= k)
    elif 'diag02' in sym:
        istrue *= (i == k)
    if 'asym03' in sym:
        istrue *= (i < l)
    elif 'sym03' in sym:
        istrue *= (i <= l)
    elif 'diag03' in sym:
        istrue *= (i == l)
    if 'asym12' in sym:
        istrue *= (j < k)
    elif 'sym12' in sym:
        istrue *= (j <= k)
    elif 'diag12' in sym:
        istrue *= (j == k)
    if 'asym13' in sym:
        istrue *= (j < l)
    elif 'sym13' in sym:
        istrue *= (j <= l)
    elif 'diag13' in sym:
        istrue *= (j == l)
    if 'asym23' in sym:
        istrue *= (k < l)
    elif 'sym23' in sym:
        istrue *= (k <= l)
    elif 'diag23' in sym:
        istrue *= (k == l)
    if sym == 'nosym':
        istrue = True
    return istrue

def zerorank1(DIM=3):
    return [sp.sympify(0) for i in range(DIM)]

def zerorank2(DIM=3):
    return [[sp.sympify(0) for i in range(DIM)] for j in range(DIM)]

def zerorank3(DIM=3):
    return [[[sp.sympify(0) for i in range(DIM)] for j in range(DIM)] for k in range(DIM)]

def zerorank4(DIM=3):
    return [[[[sp.sympify(0) for i in range(DIM)] for j in range(DIM)] for k in range(DIM)] for l in range(DIM)]

def declarerank1(objname, DIM=3):
    return [sp.sympify(objname + "_" + str(i)) for i in range(DIM)]

def declarerank2(objname, symmetry_option = 'nosym', DIM=3):
    IDX_OBJ_TMP = [[sp.sympify(objname + "_" + str(i) + str(j)) for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            if symmetry_option == "sym01":
                if (j < i):
                    # j<i in g_{ij} would indicate, e.g., g_{21}.
                    #  By this convention, we must set
                    #  g_{21} = g_{12}:
                    IDX_OBJ_TMP[i][j] = IDX_OBJ_TMP[j][i]
                    
            elif symmetry_option == "asym01":
                if (j <= i):
                    IDX_OBJ_TMP[i][j] = myEij(i,j)*IDX_OBJ_TMP[j][i]
            elif symmetry_option == "diag01":
                if (i != j):
                    IDX_OBJ_TMP[i][j] = 0
            elif symmetry_option == "nosym":
                pass
            else:
                print("Error: symmetry option " + symmetry_option + " unsupported.")
                sys.exit(1)
    return IDX_OBJ_TMP

def declarerank3(objname, symmetry_option = 'nosym', DIM=3):
    IDX_OBJ_TMP = [[[sp.sympify(objname + "_" + str(i) + str(j) + str(k)) for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                if symmetry_option == 'nosym':
                    pass
                else:
                    if "asym01" in symmetry_option:
                        if j <= i:
                            IDX_OBJ_TMP[i][j][k] = myEij(i,j)*IDX_OBJ_TMP[j][i][k]
                    elif "sym01" in symmetry_option:
                        if j < i:
                            IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[j][i][k]
                    elif "diag01" in symmetry_option:
                        if j != i:
                            IDX_OBJ_TMP[i][j][k] = 0
                    if "asym12" in symmetry_option:
                        if k <= j:
                            IDX_OBJ_TMP[i][j][k] = myEij(j,k)*IDX_OBJ_TMP[i][k][j]
                    elif "sym12" in symmetry_option:
                        if k < j:
                            IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[i][k][j]
                    elif "diag12" in symmetry_option:
                        if k != j:
                            IDX_OBJ_TMP[i][j][k] = 0
                    if "asym02" in symmetry_option:
                        if k <= i:
                            IDX_OBJ_TMP[i][j][k] = myEij(i,j)*IDX_OBJ_TMP[k][j][i]
                    elif "sym02" in symmetry_option:
                        if k < i:
                            IDX_OBJ_TMP[i][j][k] = IDX_OBJ_TMP[k][j][i]  
                    elif "diag02" in symmetry_option:
                        if k != i:
                            IDX_OBJ_TMP[i][j][k] = 0
                        
                    
    return IDX_OBJ_TMP

def declarerank4(objname, symmetry_option= 'nosym', DIM=3):
    IDX_OBJ_TMP = [[[[sp.sympify(objname + "_" + str(i) + str(j) + str(k) + str(l)) for l in range(DIM)] for k in range(DIM)] for j in range(DIM)] for i in range(DIM)]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    if symmetry_option == 'nosym':
                        pass
                    else:
                        if "asym01" in symmetry_option:
                            if(j <= i):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(i,j)*IDX_OBJ_TMP[j][i][k][l]
                        elif "sym01" in symmetry_option:
                            if(j < i):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[j][i][k][l]
                        elif "diag01" in symmetry_option:
                            if(j != i):
                                IDX_OBJ_TMP[i][j][k][l] = 0
                        if "asym02" in symmetry_option:
                            if(k <= i):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(i,k)*IDX_OBJ_TMP[k][j][i][l]
                        elif "sym02" in symmetry_option:
                            if(k < i):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[k][j][i][l]
                        elif "diag02" in symmetry_option:
                            if(k != i):
                                IDX_OBJ_TMP[i][j][k][l] = 0
                        if "asym03" in symmetry_option:
                            if(l <= i):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(i,l)*IDX_OBJ_TMP[l][j][k][i]
                        elif "sym03" in symmetry_option:
                            if(l < i):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[l][j][k][i]
                        elif "diag03" in symmetry_option:
                            if(l != i):
                                IDX_OBJ_TMP[i][j][k][l] = 0
                        if "asym12" in symmetry_option:
                            if(k <= j):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(j,k)*IDX_OBJ_TMP[i][k][j][l]
                        elif "sym12" in symmetry_option:
                            if(k < j):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][k][j][l]
                        elif "diag12" in symmetry_option:
                            if(k != j):
                                IDX_OBJ_TMP[i][j][k][l] = 0
                        if "asym13" in symmetry_option:
                            if(l <= j):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(j,l)*IDX_OBJ_TMP[i][l][k][j]
                        elif "sym13" in symmetry_option:
                            if(l < j):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][l][k][j]
                        elif "diag13" in symmetry_option:
                            if(l != j):
                                IDX_OBJ_TMP[i][j][k][l] = 0
                        if "asym23" in symmetry_option:
                            if(l <= k):
                                IDX_OBJ_TMP[i][j][k][l] = myEij(k,l)*IDX_OBJ_TMP[i][j][l][k]
                        elif "sym23" in symmetry_option:
                            if(l < k):
                                IDX_OBJ_TMP[i][j][k][l] = IDX_OBJ_TMP[i][j][l][k]
                        elif "diag23" in symmetry_option:
                            if(l != k):
                                IDX_OBJ_TMP[i][j][k][l] = 0
    return IDX_OBJ_TMP

def declarerank(objname, symmetry_option = 'nosym', DIM = 3):
    if objname.find('_') == -1:
        return sp.sympify(objname)
    else:
        index_type_string = objname[objname.find('_')+1:]
        
        if len(index_type_string) == 1:
            return declarerank1(objname, DIM)
        elif len(index_type_string) == 2:
            return declarerank2(objname, symmetry_option, DIM)
        elif len(index_type_string) == 3:
            return declarerank3(objname, symmetry_option, DIM)
        elif len(index_type_string) == 4:
            return declarerank4(objname, symmetry_option, DIM)
        else:
            print("Only tensors up to rank 4 are supported")

def symm_matrix_inverter3x3(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using symmetry_optionPy's built-in functions, since the matrix is symmetric.
    outDET = -a[0][2]**2*a[1][1] + 2*a[0][1]*a[0][2]*a[1][2] - \
                a[0][0]*a[1][2]**2 - a[0][1]**2*a[2][2] + \
                a[0][0]*a[1][1]*a[2][2]

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    # First fill in the upper-triangle of the gPhysINV matrix...
    outINV[0][0] = (-a[1][2]**2              + a[1][1]*a[2][2])/outDET
    outINV[0][1] = (+a[0][2]*a[1][2] - a[0][1]*a[2][2])/outDET
    outINV[0][2] = (-a[0][2]*a[1][1] + a[0][1]*a[1][2])/outDET
    outINV[1][1] = (-a[0][2]**2              + a[0][0]*a[2][2])/outDET
    outINV[1][2] = (+a[0][1]*a[0][2] - a[0][0]*a[1][2])/outDET
    outINV[2][2] = (-a[0][1]**2              + a[0][0]*a[1][1])/outDET
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    return outINV

def symm_matrix_det3x3(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using symmetry_optionPy's built-in functions, since the matrix is symmetric.
    outDET = -a[0][2]**2*a[1][1] + 2*a[0][1]*a[0][2]*a[1][2] - \
                a[0][0]*a[1][2]**2 - a[0][1]**2*a[2][2] + \
                a[0][0]*a[1][1]*a[2][2]

    return outDET

def symm_matrix_inverter4x4(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using symmetry_optionPy's built-in functions, since the matrix is symmetric.
    outDET = + a[0][2]*a[0][2]*a[1][3]*a[1][3] + a[0][3]*a[0][3]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][3]*a[2][3] \
             - a[0][0]*a[1][3]*a[1][3]*a[2][2] - a[0][3]*a[0][3]*a[1][1]*a[2][2] - a[0][0]*a[1][1]*a[2][3]*a[2][3] \
             - 2*(+ a[0][1]*a[0][2]*a[1][3]*a[2][3] - a[0][0]*a[1][2]*a[1][3]*a[2][3]                              \
                  - a[0][3]*(- a[0][2]*a[1][2]*a[1][3] + a[0][1]*a[1][3]*a[2][2]                                   \
                             + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3]))                                 \
             - a[3][3] * (+ a[0][2]*a[0][2]*a[1][1] - a[0][1]*a[0][2]*a[1][2] - a[0][1]*a[0][2]*a[1][2]            \
                          + a[0][0]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][2] - a[0][0]*a[1][1]*a[2][2])

    outINV = [[sp.sympify(0) for i in range(4)] for j in range(4)]

    # First fill in the upper-triangle of the gPhysINV matrix...
              
    outINV[0][0] = (-a[1][3]*a[1][3]*a[2][2] + 2*a[1][2]*a[1][3]*a[2][3] - a[1][1]*a[2][3]*a[2][3] - a[1][2]*a[1][2]*a[3][3] + a[1][1]*a[2][2]*a[3][3])/outDET
    outINV[1][1] = (-a[0][3]*a[0][3]*a[2][2] + 2*a[0][2]*a[0][3]*a[2][3] - a[0][0]*a[2][3]*a[2][3] - a[0][2]*a[0][2]*a[3][3] + a[0][0]*a[2][2]*a[3][3])/outDET
    outINV[2][2] = (-a[0][3]*a[0][3]*a[1][1] + 2*a[0][1]*a[0][3]*a[1][3] - a[0][0]*a[1][3]*a[1][3] - a[0][1]*a[0][1]*a[3][3] + a[0][0]*a[1][1]*a[3][3])/outDET
    outINV[3][3] = (-a[0][2]*a[0][2]*a[1][1] + 2*a[0][1]*a[0][2]*a[1][2] - a[0][0]*a[1][2]*a[1][2] - a[0][1]*a[0][1]*a[2][2] + a[0][0]*a[1][1]*a[2][2])/outDET
    outINV[0][1] = (+a[0][3]*a[1][3]*a[2][2] -   a[0][3]*a[1][2]*a[2][3] - a[0][2]*a[1][3]*a[2][3] + a[0][1]*a[2][3]*a[2][3] + a[0][2]*a[1][2]*a[3][3] - a[0][1]*a[2][2]*a[3][3])/outDET
    outINV[0][2] = (-a[0][3]*a[1][2]*a[1][3] +   a[0][2]*a[1][3]*a[1][3] + a[0][3]*a[1][1]*a[2][3] - a[0][1]*a[1][3]*a[2][3] - a[0][2]*a[1][1]*a[3][3] + a[0][1]*a[1][2]*a[3][3])/outDET
    outINV[0][3] = (-a[0][2]*a[1][2]*a[1][3] +   a[0][1]*a[1][3]*a[2][2] + a[0][3]*a[1][2]*a[1][2] - a[0][3]*a[1][1]*a[2][2] + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3])/outDET
    outINV[1][2] = (+a[0][3]*a[0][3]*a[1][2] +   a[0][0]*a[1][3]*a[2][3] - a[0][3]*a[0][2]*a[1][3] - a[0][3]*a[0][1]*a[2][3] + a[0][1]*a[0][2]*a[3][3] - a[0][0]*a[1][2]*a[3][3])/outDET
    outINV[1][3] = (+a[0][2]*a[0][2]*a[1][3] +   a[0][1]*a[0][3]*a[2][2] - a[0][0]*a[1][3]*a[2][2] + a[0][0]*a[1][2]*a[2][3] - a[0][2]*a[0][3]*a[1][2] - a[0][2]*a[0][1]*a[2][3])/outDET
    outINV[2][3] = (+a[0][2]*a[0][3]*a[1][1] -   a[0][1]*a[0][3]*a[1][2] - a[0][1]*a[0][2]*a[1][3] + a[0][0]*a[1][2]*a[1][3] + a[0][1]*a[0][1]*a[2][3] - a[0][0]*a[1][1]*a[2][3])/outDET

    # Then we fill the lower triangle of the symmetric matrix
    outINV[1][0] = outINV[0][1]
    outINV[2][0] = outINV[0][2]
    outINV[2][1] = outINV[1][2]
    outINV[3][0] = outINV[0][3]
    outINV[3][1] = outINV[1][3]
    outINV[3][2] = outINV[2][3]
    
    return outINV

def symm_matrix_det4x4(a):
    # It is far more efficient to write out the matrix determinant and inverse by hand
    #   instead of using symmetry_optionPy's built-in functions, since the matrix is symmetric.
    outDET = + a[0][2]*a[0][2]*a[1][3]*a[1][3] + a[0][3]*a[0][3]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][3]*a[2][3] \
             - a[0][0]*a[1][3]*a[1][3]*a[2][2] - a[0][3]*a[0][3]*a[1][1]*a[2][2] - a[0][0]*a[1][1]*a[2][3]*a[2][3] \
             - 2*(+ a[0][1]*a[0][2]*a[1][3]*a[2][3] - a[0][0]*a[1][2]*a[1][3]*a[2][3]                              \
                  - a[0][3]*(- a[0][2]*a[1][2]*a[1][3] + a[0][1]*a[1][3]*a[2][2]                                   \
                             + a[0][2]*a[1][1]*a[2][3] - a[0][1]*a[1][2]*a[2][3]))                                 \
             - a[3][3] * (+ a[0][2]*a[0][2]*a[1][1] - a[0][1]*a[0][2]*a[1][2] - a[0][1]*a[0][2]*a[1][2]            \
                          + a[0][0]*a[1][2]*a[1][2] + a[0][1]*a[0][1]*a[2][2] - a[0][0]*a[1][1]*a[2][2])
    
    return outDET


# symmetry_optionPy's generic matrix inverter is highly inefficient for 3x3 matrices, so here we have an optimized version.
def generic_matrix_inverter3x3(a):
    outDET = -a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] + \
              a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - \
              a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]

    outINV = [[sp.sympify(0) for i in range(3)] for j in range(3)]

    outINV[0][0] = -a[1][2]*a[2][1] + a[1][1]*a[2][2]
    outINV[0][1] =  a[0][2]*a[2][1] - a[0][1]*a[2][2]
    outINV[0][2] = -a[0][2]*a[1][1] + a[0][1]*a[1][2]
    outINV[1][0] =  a[1][2]*a[2][0] - a[1][0]*a[2][2]
    outINV[1][1] = -a[0][2]*a[2][0] + a[0][0]*a[2][2]
    outINV[1][2] =  a[0][2]*a[1][0] - a[0][0]*a[1][2]
    outINV[2][0] = -a[1][1]*a[2][0] + a[1][0]*a[2][1]
    outINV[2][1] =  a[0][1]*a[2][0] - a[0][0]*a[2][1]
    outINV[2][2] = -a[0][1]*a[1][0] + a[0][0]*a[1][1]

    for i in range(3):
        for j in range(3):
            outINV[i][j] /= outDET

    return outINV

def generic_matrix_det3x3(a):
    outDET = -a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] + \
              a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] - \
              a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]


    return outDET
