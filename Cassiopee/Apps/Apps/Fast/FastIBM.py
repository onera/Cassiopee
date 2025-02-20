from Generator.IBM import buildOctree, generateIBMMesh, createRefinementBodies

from Connector.IBM import prepareIBMData, dist2wallIBM, blankingIBM, buildFrontIBM, setInterpDataIBM, initializeIBM

from Geom.IBM import setSnear, _setSnear, setDfar, _setDfar, setIBCType, _setIBCType, _setFluidInside, _setFluidOutside, snearFactor, _snearFactor, setIBCType, changeIBCType, _changeIBCType, initOutflow, _initOutflow, initInj, _initInj, setFluidInside, setFluidOutside, flatPlate, bumpInChannel, naca0012

from Generator.IBMmodelHeight import computeModelisationHeight, computeSnearOpt

from Post.IBM import extractPressureGradients, extractPressureHighOrder, extractYplusAtImagePoints, prepareSkinReconstruction, computeSkinVariables, _computeSkinVariables, computeAerodynamicLoads, computeAerodynamicCoefficients

import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Converter.Filter as Filter
import Generator.PyTree as G
import Post.PyTree as P
import Geom.PyTree as D
import Post.Mpi as Pmpi
import Transform.PyTree as T
import Intersector.PyTree as XOR
import Post.Probe as Probe
import Distributor2.PyTree as D2

import numpy
import math
import os

###############
# Pressure outlet control
###############
# For more info : Mouton, S. (2023). Numerical Simulation of the Flow in the ONERA F1 Wind Tunnel. Journal of Aircraft, 60(2), 355-367.

# ∆p_(2n) = G∆M-(1n) + D[∆M_(1n) − ∆M_(n−ncontrol) ]
# G = -0.7 *  ∂p2/∂M1 || G: proportional coefficeint
# 0.7 :: damping coefficient for stability
# ∂p2/∂M1 ≈ −p_(2,is) γ M_(2,is) A_2 Σ′(M_1)/(τ(M_(2,is)) A_1 Σ′(M_(2,is)))
# D = G(n_charac/n_control)
# n_charac = approx. characteristic response time of the CFD of the wind tunnel - a prior run to determine this
# τ :: tau_function
# Σ′:: sigma_function_inv

def tauFunction__(mach, gamma=1.4):
    """Return :math:`1 + \\frac{\\gamma - 1}{2} M^2`

    Parameters
    ----------
    mach : array_like
        Value of Mach number :math:`M`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    """
    return 1. + 0.5*(gamma - 1.)*mach**2

def sigmaFunction__(mach, gamma=1.4):
    """Return :math:`\\Sigma(M)`

    Parameters
    ----------
    mach : array_like
        Value of Mach number :math:`M`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    """
    return ( 2./(gamma + 1.) + (gamma - 1.)/(gamma + 1.) * mach**2 )**(0.5*(gamma + 1.)/(gamma - 1.)) / mach

def sigmaFunctionDerivate__(mach, gamma=1.4):
    """Return :math:`\\Sigma'(M)`, first derivative of :math:`\\Sigma(M)`

    Parameters
    ----------
    mach : array_like
        Value of Mach number :math:`M`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    """
    return ( 2./(gamma + 1.) + (gamma - 1.)/(gamma + 1.) * mach**2 )**(0.5*(gamma + 1.)/(gamma - 1.)) \
        * (-1/mach**2 + 1./(2./(gamma + 1.) + (gamma - 1.)/(gamma + 1.) * mach**2))

def _scalarSigmaFunctionInv__(s, gamma=1.4, range='subsonic'):
    import scipy.optimize
    eps = numpy.finfo(float).eps
    if range == 'subsonic':
        sol = scipy.optimize.root_scalar(lambda M: sigmaFunction__(M, gamma) - s,
                                         x0=0.5, bracket=(2*eps, 1. - 2*eps), method='brentq')
    elif range == 'supersonic':
        sol = scipy.optimize.root_scalar(lambda M: sigmaFunction__(M, gamma) - s,
                                         x0=1.5, bracket=(1. + 2*eps, 1e3), method='brentq') # it is unlikely that user require Mach number above 1000.
    else:
        raise RuntimeError("Unexpected value for `range`: {:s}".format(str(range)))
    return sol.root

def sigmaFunctionInv__(s, gamma=1.4, range='subsonic'):
    # This method vectorizes _scalarSigmaFunctionInv__
    """Return the inverse of the function :math:`\\Sigma(M)`

    Parameters
    ----------
    s : array_like
        Value of :math:`\\Sigma`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    range : ['subsonic', 'supersonic']
        Range over which the inverse is to be looked for
    """
    with numpy.nditer([s, None],
                      op_flags=[['readonly'], ['writeonly', 'allocate', 'no_broadcast']],
                      op_dtypes=['float64', 'float64']) as it:
        for x, y in it:
            y[...] = _scalarSigmaFunctionInv__(s, gamma=gamma)
        return it.operands[1]

def getInfoOutletPressure(tcase, familyName, gamma=1.4):

    familyNameExists=False
    for z in Internal.getZones(tcase):
        FamNode = Internal.getNodeFromType(z, 'FamilyName_t')
        if FamNode is not None:
            FamName = Internal.getValue(FamNode)
            if FamName == familyName:
                familyNameExists = True
                solDef           = Internal.getNodeFromName(z, '.Solver#define')
                controlProbeName = Internal.getNodeFromName(solDef, 'probeName')
                controlProbeName = Internal.getValue(controlProbeName)
                A1               = Internal.getNodeFromName(solDef, 'AtestSection')[1][0]
                A2               = Internal.getNodeFromName(solDef, 'AOutPress')[1][0]
                m1               = Internal.getNodeFromName(solDef, 'machTarget')[1][0]
                pi1              = Internal.getNodeFromName(solDef, 'pStatTarget')[1][0]
                ti1              = Internal.getNodeFromName(solDef, 'tStatTarget')[1][0]
                lbda             = Internal.getNodeFromName(solDef, 'lmbd')[1][0]
                cxSupport        = Internal.getNodeFromName(solDef, 'cxSupport')[1][0]
                sSupport         = Internal.getNodeFromName(solDef, 'sSupport')[1][0]
                itExtrctPrb      = Internal.getNodeFromName(solDef, 'itExtrctPrb')[1][0]
                break

    if not familyNameExists:
        print('Family name %s is not an outlet pressure family name'%(familyName))
        print('Exiting')
        exit()

    ## Calculated Values
    _tau        = tauFunction__(m1)                                   # Eqn. (10) || τ(M)= 1 + (γ-1)/2 M²
    m2is        = sigmaFunctionInv__(A2/A1 * sigmaFunction__(m1))     # Eqn. (11) || M_2,is=Σ⁻¹(A_2 /A_1 Σ(M_1))
    p2is        = pi1 * tauFunction__(m2is)**(-gamma/(gamma - 1.))    # Eqn. (12) || p_2,is=p_i1 τ(M_2,is)^(−γ∕(γ−1))

    # coefficient de perte de charge entre l'entrée et la sortie du domaine,
    # ici uniquement du au support, et calculé à partir d'une estimation de la traînée du support
    if lbda <0: lbda = cxSupport * sSupport/A1

    p0          = pi1 / _tau **(gamma/(gamma - 1.))                  # Pa - static pressure for m1
    q0          = 0.5*p0*gamma*m1**2                                 # Pa - dynamic pressure for m1
    p2          = p2is - q0*lbda                                     # Eqn. (13) || Pa - outlet static pressure || λ = (p_2,is-p_2)/q_1

    values4gain  =[p2,
                   m1,
                   p2is*gamma*m2is,
                   tauFunction__(m2is),
                   A2/A1,
                   sigmaFunctionDerivate__(m1),
                   sigmaFunctionDerivate__(m2is),
                   0,
                   0]

    return values4gain,controlProbeName,itExtrctPrb

def _setUpOutletPressure(values4gain,itValues4gain):
    # G : proportional gain Eqn. (15) || G=-0.7* (∂p_2/∂M_1)
    #     The 0.7 constant is for stability of the algorithm and is independent of the WT
    values4gain[7]=-0.7 * (-values4gain[2] / values4gain[3] * values4gain[4] * values4gain[5] / values4gain[6])

    # D: derivative coefficient || D=G(n_charac/n_control)
    values4gain[8]=values4gain[7] * float(itValues4gain[2]) / itValues4gain[1]
    return None

def _controlOutletPressureMachProbe(tc,dctProbes,controlProbeName,values4gain,it,familyName):
    probe_tmp = dctProbes[controlProbeName]
    proc = probe_tmp._proc

    if Cmpi.rank == proc:
        print("   iteration {:06d}: adjusting back pressure...".format(it), end='')
        time      = numpy.concatenate([Internal.getNodeFromName(z,'time')[1] for z in probe_tmp._probeZones])
        mach      = numpy.concatenate([Internal.getNodeFromName(z,'Mach')[1] for z in probe_tmp._probeZones])

        select    = (time >= 0.)  # select values that have been recorded: array is prefilled with -1
        time      = time[select]
        mach      = mach[select]

        index_current  = -1
        index_previous = -2 #-1 - (itValues4gain[1] // itExtractProbe)
        current_it     = time[index_current]
        current_mach   = mach[index_current]

        previous_it    = time[index_previous]
        previous_mach  = mach[index_previous]
        #print(current_it, previous_it, current_mach, previous_mach,flush=True)
        oldPressure    = values4gain[0]
        values4gain[0] += values4gain[7] * (current_mach - values4gain[1]) + values4gain[8] * (current_mach - previous_mach)
        print("Control of Output Pressure:: target Ma = {:.5f}, current Ma = {:.5f}|| old Pout = {:.1f} Pa, new Pout = {:.1f} Pa, ".format(values4gain[1], current_mach, oldPressure, values4gain[0]))
    Cmpi.bcast(values4gain[0], root=proc)
    _initOutflow(tc, familyName, values4gain[0])
    Cmpi.barrier()
    return None

###############
# Point Probes
###############

def createPointProbes(probe_in, probePointsList):
    listOfZones = []

    for cpt, (x,y,z) in enumerate(probePointsList):
        point = D.point((x,y,z))
        point[0] = "point_{:03d}".format(cpt)
        listOfZones.append(point)

    probe = C.newPyTree(['Base', listOfZones])

    if Cmpi.rank == 0: C.convertPyTree2File(probe, probe_in)

    Cmpi.barrier()

    return None

def initPointProbes(t, probe_in, fields, bufferSize=100, append=False, historyDirectory='.'):
    if isinstance(probe_in, str): probes = C.convertFile2PyTree(probe_in)
    else: probes = Internal.copyTree(probe_in)

    fields = ['centers:'+fname if 'centers' not in fname else fname for fname in fields] #extraction from cell-centered t

    dictOfProbes = {}
    for z in Internal.getZones(probes):
        name = z[0]
        xnode = Internal.getNodeFromName(z,'CoordinateX')
        ynode = Internal.getNodeFromName(z,'CoordinateY')
        znode = Internal.getNodeFromName(z,'CoordinateZ')
        point = [Internal.getValue(xnode),Internal.getValue(ynode),Internal.getValue(znode)]
        dictOfProbes[name] = point

    if 'centers:Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for pname in dictOfProbes.keys():
        point = dictOfProbes[pname]
        filename = "probe_{:s}.cgns".format(pname)
        filename = os.path.join(historyDirectory, filename)

        probe = Probe.Probe(filename, t, X=point, fields=fields, bufferSize=bufferSize, append=append)
        dictOfProbes[pname] = probe

    if 'centers:Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return dictOfProbes

def _updatePointProbes(t, dictOfProbes, it, fields):
    fields = ['centers:'+fname if 'centers' not in fname else fname for fname in fields] #extraction from cell-centered t

    if 'centers:Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for name, probe in dictOfProbes.items():
        probe.extract(t, time=it)

    if 'centers:Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'centers:Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return None

###############
# Surface Probes
###############

def generateIsoXSurface__(tb, x):
    bbox = G.bbox(tb)
    alpha= 0.05
    DY = bbox[4]-bbox[1]; DZ = bbox[5]-bbox[2]
    YMIN = bbox[1]-alpha*DY ; ZMIN = bbox[2]-alpha*DY
    LY = DY + 2*alpha*DY ; LZ = DZ + 2*alpha*DZ
    NJ = 51; NK = 51
    a = G.cart((x, YMIN, ZMIN), (1, LY/(NJ-1), LZ/(NK-1)), (1, NJ, NK))
    a = C.convertArray2Tetra(a)
    zones = Internal.getZones(tb)+[a]
    z = T.join(zones)
    z = XOR.conformUnstr(z, tol=1e-10, itermax=1)
    zones = T.splitManifold(z)
    candidates = []
    eps = 1e-4
    for z in zones:
        bboxz = G.bbox(z)
        if abs(bboxz[0]-x)<eps and abs(bboxz[3]-x)<eps:
            if bboxz[1]>=bbox[1]-eps and bboxz[4]<=bbox[4]+eps and bboxz[2]>=bbox[2]-eps and bboxz[5]<=bbox[5]+eps:
                candidates.append(z)
    for i,z in enumerate(candidates):
        z[0] = "zone_{:02d}".format(i)
    candidates = T.join(candidates)
    surface = C.newPyTree(["X_{:5.3f}".format(x), candidates])
    return surface

def createSurfaceProbes(tb, surface_in, probeSurfaceList):
    surfaces = []

    for x in probeSurfaceList:
        ts = generateIsoXSurface__(tb, x)
        G._getSmoothNormalMap(ts)
        surfaces.append(ts)

    ts = Internal.merge(surfaces)

    if Cmpi.size > 1:
        for b in Internal.getBases(ts):
            T._splitNParts(b, Cmpi.size)
            D2._distribute(b, Cmpi.size)

    if Cmpi.rank == 0: C.convertPyTree2File(ts, surface_in)

    Cmpi.barrier()

    return None

def initSurfaceProbes(t, tc, surface_in, fields, bufferSize=100, historyDirectory='.'):
    if isinstance(surface_in, str):
        if Cmpi.size > 1: probes = Cmpi.convertFile2PyTree(surface_in, proc=Cmpi.rank)
        else: probes = C.convertFile2PyTree(surface_in)
    else: probes = Internal.copyTree(surface_in)

    dictOfProbes = {}

    if 'Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    tcs = Internal.rmNodesFromType(tc, 'ZoneSubRegion_t')
    if Cmpi.size <= 1: Cmpi._setProc(tcs, 0) # Security for Probe functions
    for var in fields+['cellN']: C._cpVars(t, 'centers:'+var, tcs, var)

    for b in Internal.getBases(probes):
        sname = b[0]
        tbs_loc = C.newPyTree([sname, Internal.getZones(b)])

        filename = "surface_{:s}.cgns".format(sname)
        filename = os.path.join(historyDirectory, filename)

        probe = Probe.Probe(filename, tPermeable=tbs_loc, fields=fields, bufferSize=bufferSize)
        tcs_loc = probe.prepare(tcs)

        dictOfProbes[sname] = [probe, tbs_loc, tcs_loc]

    Cmpi.barrier()

    return dictOfProbes

def _updateSurfaceProbes(t, dictOfProbes, fields):
    if 'Mach' in fields: P._computeVariables(t, ['centers:Mach'])
    if 'Pressure' in fields: P._computeVariables(t, ['centers:Pressure'])

    for key in dictOfProbes:
        probe, tbs_loc, tcs_loc = dictOfProbes[key]
        for var in fields:
            C._cpVars(t, 'centers:'+var, tcs_loc, var)
            C._initVars(tbs_loc, var, 1)

        probe.extract(tcs_loc, time=1, onlyTransfer=True)

    if 'Mach' in fields: C._rmVars(t, ['centers:Mach'])
    if 'Pressure' in fields: C._rmVars(t, ['centers:Pressure'])

    Cmpi.barrier()

    return None

def getMassflow(t):
    C._initVars(t, '{massflow}={Density}*({sx}*{VelocityX}+{sy}*{VelocityY}+{sz}*{VelocityZ})')
    massflow = abs(Pmpi.integ(t, 'massflow')[0])
    return massflow

def integrateSurfaceProbes(dictOfProbes):
    massflows = []
    for key in dictOfProbes:
        probe, tbs_loc, tcs_loc = dictOfProbes[key]
        massflow_loc = getMassflow(tbs_loc) # intégration de la masse volumique sur chaque surface probe
        massflows.append(massflow_loc)

    Cmpi.barrier()

    return massflows
