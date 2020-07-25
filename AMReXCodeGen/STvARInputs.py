def StoppingCriteria(max_step = 10, max_time = 100.0):
    header = "# Stopping Criteria \n"
    string = "max_step = " + str(max_step) + "\n"
    string += "stop_time = " + str(max_time) + "\n\n"
    return header + string

def Geometry(NumCells = [32, 32, 32], LowerBound = [-10, -10, -10], UpperBound = [10, 10, 10], PeriodicBCs = [0, 0, 0], CoordSys = "Cartesian", dim = 3):
    header = "# Problem Size & Geometry \n"
    ncellstring = "amr.n_cell = "
    geomperiodicstring = "geometry.is_periodic = "
    lostring = "geometry.prob_lo = "
    histring = "geometry.prob_hi = "
    for i in range(dim):
        ncellstring += str(NumCells[i]) + " "
        geomperiodicstring += str(PeriodicBCs[i]) + " "
        lostring += str(LowerBound[i]) + " "
        histring += str(UpperBound[i]) + " "
    ncellstring += "\n\n"
    geomperiodicstring += "\n"
    lostring += "\n"
    histring += "\n"
    return header + ncellstring + geomperiodicstring + lostring + histring + "\n"

def AMRVerbosity(Verbosity = 1):
    return "# Turn on verbosity in amr \n" + "amr.v = " + str(Verbosity) + "\n\n"

def RefinementPar(NumLevels = 1, MaxGridCells = 32, RefRatio = 2, NumBufferCells = 1, BlockingFactor = 8, RegridPeriod = 2):
    header = "# Refinement Parameters \n"
    numlevstring = "amr.max_level = " + str(NumLevels-1) + "\n"
    maxgridcellsstring = "amr.max_grid_size = " + str(MaxGridCells) + "\n"
    refratiostring = "amr.ref_ratio = " + str(RefRatio) + "\n"
    buffstring = "amr.n_error_buf = " + str(NumBufferCells) + "\n"
    blockstring = "amr.blocking_factor = " + str(BlockingFactor) + "\n"
    regridperstring = "amr.regrid_int = " + str(RegridPeriod) + "\n"
    return header + numlevstring + maxgridcellsstring + refratiostring + buffstring + blockstring + regridperstring + "\n"

def Interpolation(Type = 6):
    header = "# AMR Interpolation \n"
    typestring = "amr.interpolation_type = " +str(Type) + "\n"
    return header + typestring + "\n"

def BoundaryConditions(BCLO = [0, 0, 0], BCHI = [0, 0, 0], dim = 3):
    header = "# Problem specific boundary conditions \n"
    lostring = "domain_lo_bc_types = "
    histring = "domain_hi_bc_types = "
    for i in range(dim):
        lostring += str(BCLO[i]) + " "
        histring += str(BCLO[i]) + " "
    lostring += "\n"
    histring += "\n" 
    return header + lostring + histring + "\n"
        
def TaggingThresholds(TaggingArray, ErrorComponent = 0):
    header = "# Problem specific tagging for refinement \n"
    string = "problem.s_error = "
    for i in range(len(TaggingArray)):
        string += str(TaggingArray[i]) + " "
    string += "\n"
    string += "problem.error_comp = " + str(ErrorComponent) + "\n"
    return header + string + "\n"

def ProblemType(cfl = 0.25, is_elliptic = 0):
    header = "# Problem specific inputs \n"
    ellipticstring = "problem.elliptic = " + str(is_elliptic) + "\n"
    cflstring = "problem.cfl = " + str(cfl) + "\n"
    return header + ellipticstring + cflstring + "\n"
    
def PlottingCriteria(PlotPeriod = 1, DiagPeriod = 1, CheckpointPeriod = 10):
    header = "# Plotting steps \n"
    pltstring = "amr.plot_int = " + str(PlotPeriod) + "\n"
    diagstring = "amr.diag_int = " + str(DiagPeriod) + "\n"
    chkstring = "amr.chk_int = " + str(CheckpointPeriod) + "\n"
    return header + pltstring + diagstring + chkstring + "\n"
    
def RestartFromCheckPoint(restartfile, is_initial_data = 0):
    header = "# Restart from checkpoint or use as initial condition? \n"
    restartfilestring = "amr.restart = " + str(restartfile) + "\n"
    isinitialstring = "amr.restart_is_initial_data = " + str(is_initial_data) + "\n"
    return header + restartfilestring + isinitialstring + "\n"
    
def IntegrationType(RKType = 1, rkweights = 1, rknodes = 0, rktableau = 0.0, is_fpe_trap_invalid = 1):
    header1 = "## integration.type can take on the following values: \n## 0 = Forward Euler \n## 1 = Explicit Runge Kutta \n"
    if RKType == 1:
        typestring = "integration.type = 0 \n\n"
    else:
        typestring = "integration.type = 1 \n\n"
    r1 = header1 + typestring
    header2 = "## Explicit Runge-Kuta parameters \n# \n## integration.rk.type can take the following values:\n### 0 = User-specified Butcher Tableau\n### 1 = Forward Euler\n### 2 = Trapezoid Method\n### 3 = SSPRK3 Method\n### 4 = RK4 Method \n"
    rktypestring = "integration.rk.type = " + str(RKType) + "\n\n"
    r2 = header2 + rktypestring
    header3 = "## If using a user-specified Butcher Tableau, then\n## set nodes, weights, and table entries here:\n#\n## The Butcher Tableau is read as a flattened,\n## lower triangular matrix (but including the diagonal)\n## in row major format. \n"
    wntstring = "integration.rk.weights = " + str(rkweights) + "\n"
    wntstring += "integration.rk.nodes = " + str(rknodes) + "\n"
    wntstring += "integration.rk.tableau = " + str(rktableau) + "\n\n"
    wntstring += "amrex.fpe_trap_invalid = " + str(is_fpe_trap_invalid) + "\n"
    r3 = header3 + wntstring
    return r1 + r2 + r3
        
    
    
    
    
    
    