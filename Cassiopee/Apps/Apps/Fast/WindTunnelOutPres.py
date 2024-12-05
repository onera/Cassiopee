import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Post.Probe as Probe
import Post.PyTree as Post
import Geom.IBM as D_IBM
import sys, os
import numpy

R_GAZ       = 287.05       # J/kg/K
GAMMA       = 1.4

def tauFunction__(mach, gamma=GAMMA):
    """Return :math:`1 + \\frac{\\gamma - 1}{2} M^2`

    Parameters
    ----------
    mach : array_like
        Value of Mach number :math:`M`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    """
    return 1. + 0.5*(gamma - 1.)*mach**2


def sigmaFunction__(mach, gamma=GAMMA):
    """Return :math:`\\Sigma(M)`

    Parameters
    ----------
    mach : array_like
        Value of Mach number :math:`M`
    gamma : float
        Specific heat ratio :math:`\\gamma`
    """
    return ( 2./(gamma + 1.) + (gamma - 1.)/(gamma + 1.) * mach**2 )**(0.5*(gamma + 1.)/(gamma - 1.)) / mach


def sigmaFunctionDerivate__(mach, gamma=GAMMA):
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


def _scalarSigmaFunctionInv__(s, gamma=GAMMA, range='subsonic'):
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



## ∆p_(2n) = G∆M-(1n) + D[∆M_(1n) − ∆M_(n−ncontrol) ]
## G = -0.7 *  ∂p2/∂M1 || G: proportional coefficeint
## 0.7 :: damping coefficient for stability
## ∂p2/∂M1 ≈ −p_(2,is) γ M_(2,is) A_2 Σ′(M_1)/(τ(M_(2,is)) A_1 Σ′(M_(2,is)))
## D = G(n_charac/n_control)
## n_charac = approx. characteristic response time of the CFD of the wind tunnel - a prior run to determine this
## τ :: tau_function
## Σ′:: sigma_function_inv

def getInfo(tcase,familyName):
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
    p2is        = pi1 * tauFunction__(m2is)**(-GAMMA/(GAMMA - 1.))    # Eqn. (12) || p_2,is=p_i1 τ(M_2,is)^(−γ∕(γ−1))

    # coefficient de perte de charge entre l'entrée et la sortie du domaine,
    # ici uniquement du au support, et calculé à partir d'une estimation de la traînée du support
    if lbda <0: lbda = cxSupport * sSupport/A1

    p0          = pi1 / _tau **(GAMMA/(GAMMA - 1.))                  # Pa - static pressure for m1
    q0          = 0.5*p0*GAMMA*m1**2                                 # Pa - dynamic pressure for m1
    p2          = p2is - q0*lbda                                     # Eqn. (13) || Pa - outlet static pressure || λ = (p_2,is-p_2)/q_1

    values4gain  =[p2,
                   m1,
                   p2is*GAMMA*m2is,
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


def getPointsFromTree(tree):
    dct_points = {}
    for z in Internal.getZones(tree):
        name = z[0]
        xnode = Internal.getNodeFromName(z,'CoordinateX')
        ynode = Internal.getNodeFromName(z,'CoordinateY')
        znode = Internal.getNodeFromName(z,'CoordinateZ')
        point = [Internal.getValue(xnode),Internal.getValue(ynode),Internal.getValue(znode)]
        dct_points[name] = point
    return dct_points


def setupMachProbe(t,buffer_size,isRestart,DIRECTORY_PROBES):
    Post._computeVariables(t, ['centers:Mach'])

    dct_probe_point       = {}
    dct_points_for_probes = getPointsFromTree(C.convertFile2PyTree(os.path.join(DIRECTORY_PROBES, "probes.cgns")))

    for name, point in dct_points_for_probes.items():
        probe = Probe.Probe(os.path.join(DIRECTORY_PROBES, "probe_{:s}.cgns".format(name)), t, X=point, fields=['centers:Mach'], bufferSize=buffer_size, append=isRestart)
        dct_probe_point[name] = probe

    C._rmVars(t, ['centers:Mach'])
    return dct_points_for_probes,dct_probe_point


def recordDataMach(t,dct_probe_point,it):
    Post._computeVariables(t, ['centers:Mach'])
    for name, probe in dct_probe_point.items(): probe.extract(t, time=it)
    C._rmVars(t, ['centers:Mach'])
    return dct_probe_point


def _controlOutletPressureMachProbe(tc,dctProbes,controlProbeName,DIRECTORY_PROBES,itValues4gain,values4gain,itExtractProbe,it,familyName):
    for name, probe in dctProbes.items():
        probe.flush()
    Cmpi.barrier()

    if Cmpi.rank == 0:
        print("   iteration {:06d}: adjusting back pressure...".format(it), end='')
        probe_tmp = C.convertFile2PyTree(os.path.join(DIRECTORY_PROBES, "probe_{:s}.cgns".format(controlProbeName)))
        time      = numpy.concatenate([node[1] for node in Internal.getNodesFromName(probe_tmp,'time')])
        mach      = numpy.concatenate([node[1] for node in Internal.getNodesFromName(probe_tmp,'Mach')])
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
    Cmpi.bcast(values4gain[0], root=0)
    D_IBM._initOutflow(tc, familyName, values4gain[0])
    Cmpi.barrier()
    return None
