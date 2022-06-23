# - IBM specific post-processing -
from . import PyTree as P
import Connector.OversetData as XOD
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM
import Connector.connector as connector
import Converter
import Converter.Distributed as Distributed
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import copy
import math
import numpy

def _add_gradxi_P(z):
    subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
    for zsr in subRegions:
        nameSubRegion = zsr[0]
        if nameSubRegion[:4] == "IBCD":
            pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
            gradxP   = Internal.getNodeFromName(zsr, 'gradxPressure')
            if gradxP is not None:
                gradyP   = Internal.getNodeFromName(zsr, 'gradyPressure')
                gradzP   = Internal.getNodeFromName(zsr, 'gradzPressure')

            nIBC    = pressure.shape[0]
            Nlength = numpy.zeros((nIBC),numpy.float64)
            if gradxP is  None:                        
                zsr[2].append(['gradxPressure' , Nlength            , [], 'DataArray_t'])
                zsr[2].append(['gradyPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['gradzPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
    return None


def extractIBMInfo(tc_in, filename_out='IBMInfo.cgns'):
    """Extracts the geometrical information required for the IBM (i.e. wall points, target points, and image points).
    Usage: extractIBMInfo (tc_in, filename_out)"""
    if isinstance(tc_in, str): tc = Cmpi.convertFile2PyTree(tc_in)
    else: tc = tc_in

    tibm = TIBM.extractIBMInfo(tc)
    rank = Cmpi.rank
    Distributed._setProc(tibm,rank)
    if isinstance(filename_out, str): Cmpi.convertPyTree2File(tibm, filename_out)
    return tibm


#===========================================================
# IN: ts: geometry tree with solution
# IN: Sref: reference surface area
# IN: alpha: angle of attack (x-y plane)
# IN: beta: angle of attack (x-z plane)
# IN: dimPb: dimension of the problem
# IN: verbose: verbose for pressure & friction forces (& coeffcients) & total cl & total cd I/O to screen
# OUT: geometry tree with additional solution fields
# I/O screen:Cl & Cd
#===========================================================
def loads0(ts, Sref=None, alpha=0., beta=0., dimPb=3, verbose=False):
    if Sref is None:
        C._initVars(ts, '__ONE__',1.)
        Sref = P.integ(ts, '__ONE__')[0];
        C._rmVars(ts, ['__ONE__', 'centers:vol'])

    RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
    PInf     = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
    RoInf    = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
    VxInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
    VyInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
    VzInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
    VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
    VInf     = math.sqrt(VInf2)

    q      = 0.5*RoInf*VInf2
    qinv   = 1./q
    alpha  = math.radians(alpha)
    beta   = math.radians(beta)
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)
    
    #===========================
    # Compute pressure forces
    #===========================
    zw = Internal.getZones(ts)
    zw = T.join(zw)
    Internal._rmNodesFromType(ts,'Zone_t')
    ts[2][1][2].append(zw)
    FSN = Internal.getNodeFromName(ts,Internal.__FlowSolutionNodes__)
    isPresent=False
    if Internal.getNodeFromName(FSN,'Cp') is None:
        C._initVars(ts, '{Cp}=-({Pressure}-%g)*%g'%(PInf,qinv))

    res = P.integNorm(ts, 'Cp')[0]
    res = [i/Sref for i in res]
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)
    if dimPb == 3:
        cd = res[0]*calpha*cbeta + res[1]*salpha*cbeta - res[2]*sbeta
        cl = res[1]*calpha       - res[0]*salpha
    else:
        cd = res[0]*calpha + res[1]*salpha
        cl = res[1]*calpha - res[0]*salpha
    if verbose:
        print("Normalized pressure drag = %g and lift = %g"%(cd, cl))
        print("Vector of pressure loads: (Fx_P,Fy_P,Fz_P)=(",res[0],res[1],res[2],")")
    cd_press = cd ; cl_press = cl
    
    #======================================
    # Compute viscous forces
    #======================================
    C._initVars(ts, '{tau_wall}=0.')
    G._getNormalMap(ts)

    # Calcul de yplus au point image
    C._initVars(ts, '{dist_cw}=sqrt(({CoordinateX_PW}-{CoordinateX_PC})**2 + ({CoordinateY_PW}-{CoordinateY_PC})**2 + ({CoordinateZ_PW}-{CoordinateZ_PC})**2)')
    C._initVars(ts, '{dist_iw}=sqrt(({CoordinateX_PW}-{CoordinateX_PI})**2 + ({CoordinateY_PW}-{CoordinateY_PI})**2 + ({CoordinateZ_PW}-{CoordinateZ_PI})**2)')
    C._initVars(ts, '{yplus_i}=({yplus}/{dist_cw})*{dist_iw}')

    variables = ['Density', 'Cp', 'tau_wall', 'Pressure','VelocityX','VelocityY','VelocityZ']
    FSN = Internal.getNodeFromType(ts,'FlowSolution_t')
    if Internal.getNodeFromName1(FSN,'utau') is not None:
        variables += ['utau','yplus','tau_wall']
        C._initVars(ts, '{tau_wall}={Density}*{utau}*{utau}')

    isGradP = False
    if Internal.getNodeFromName1(FSN,'gradxPressure') is not None:
        variables += ['gradxPressure','gradyPressure','gradzPressure']
        isGradP = True

    if Internal.getNodeFromName1(FSN,'yplus_i') is not None:
        variables += ['yplus_i']

    ts = C.node2Center(ts, variables)
    Internal._rmNodesFromName(ts,'FlowSolution')
    C._normalize(ts, ['centers:sx','centers:sy','centers:sz'])
    # Calc tangent vector
    C._initVars(ts, '{centers:tx}={centers:VelocityX}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sx}')
    C._initVars(ts, '{centers:ty}={centers:VelocityY}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sy}')
    C._initVars(ts, '{centers:tz}={centers:VelocityZ}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sz}')
    C._normalize(ts, ['centers:tx','centers:ty','centers:tz'])

    C._initVars(ts, '{centers:tauxx}=2*{centers:tau_wall}*{centers:tx}*{centers:sx}')
    C._initVars(ts, '{centers:tauyy}=2*{centers:tau_wall}*{centers:ty}*{centers:sy}')
    C._initVars(ts, '{centers:tauzz}=2*{centers:tau_wall}*{centers:tz}*{centers:sz}')
    C._initVars(ts, '{centers:tauxy}={centers:tau_wall}*({centers:tx}*{centers:sy}+{centers:ty}*{centers:sx})')
    C._initVars(ts, '{centers:tauxz}={centers:tau_wall}*({centers:tx}*{centers:sz}+{centers:tz}*{centers:sx})')
    C._initVars(ts, '{centers:tauyz}={centers:tau_wall}*({centers:ty}*{centers:sz}+{centers:tz}*{centers:sy})')

    # Cacl friction forces
    C._initVars(ts, '{centers:Fricx}={centers:tauxx}*{centers:sx}+{centers:tauxy}*{centers:sy}+{centers:tauxz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricy}={centers:tauxy}*{centers:sx}+{centers:tauyy}*{centers:sy}+{centers:tauyz}*{centers:sz}')
    C._initVars(ts, '{centers:Fricz}={centers:tauxz}*{centers:sx}+{centers:tauyz}*{centers:sy}+{centers:tauzz}*{centers:sz}')

    # Calc coefficient of friction
    C._initVars(ts, '{centers:Cf}=(sqrt({centers:Fricx}**2+{centers:Fricy}**2+{centers:Fricz}**2))/%g'%q)

    # Update of pressure gradients 
    if isGradP:
        C._initVars(ts, '{centers:gradnP}={centers:gradxPressure}*{centers:sx}+{centers:gradyPressure}*{centers:sy}+{centers:gradzPressure}*{centers:sz}')
        C._initVars(ts, '{centers:gradtP}={centers:gradxPressure}*{centers:tx}+{centers:gradyPressure}*{centers:ty}+{centers:gradzPressure}*{centers:tz}')

    G._getVolumeMap(ts)
    effortX = P.integ(ts, 'centers:Fricx')[0]
    effortY = P.integ(ts, 'centers:Fricy')[0]
    effortZ = P.integ(ts, 'centers:Fricz')[0]

    QADIMI = 1./(q*Sref)
    if dimPb == 3:
        cd = (effortX*calpha*cbeta + effortY*salpha*cbeta - effortZ*sbeta)*QADIMI
        cl = (effortY*calpha       - effortX*salpha)*QADIMI
    else:
        cd = (effortX*calpha + effortY*salpha)*QADIMI
        cl = (effortY*calpha - effortX*salpha)*QADIMI

    vars = ['centers:sx'   , 'centers:sy'   , 'centers:sz',
            'centers:tx'   , 'centers:ty'   , 'centers:tz',
            'centers:tauxx', 'centers:tauyy', 'centers:tauzz',
            'centers:tauxy', 'centers:tauxz', 'centers:tauyz',
            'centers:Fricx', 'centers:Fricy', 'centers:Fricz']
    C._rmVars(ts, vars)

    if verbose:
        print("Normalized skin friction drag = %g and lift = %g"%(cd, cl))
        print("Vector of skin friction loads: (Fx_f,Fy_f,Fz_f)=(",effortX*QADIMI, effortY*QADIMI, effortZ*QADIMI,")")

    ##################################################
    if verbose:
        print("****************************************")
        print("Total Drag :", cd_press+cd)
        print("Total Lift :", cl_press+cl)
        print("****************************************")

    return ts


#=============================================================================
# Post efforts
# IN: t_case: geometry tree
# IN: tc_in: connectivity tree
# IN: tc2_in: second connectivity tree (when using 2 image points)
# OUT: wall_out or None: file for the output of the forces on the wall at the centers
# IN: alpha: angle for the computation of the forces
# IN: beta: angle for the computation of the forces
# IN: gradP: calculate the pressure gradient
# IN: order: order of the extrapolation of pressure
# IN: Sref: reference area

# NOTE: if tc_in = None, t_case is the geometry tree with the projected solution
#==============================================================================
def loads(t_case, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., gradP=False, order=1, Sref=None, famZones=[]):
    """Computes the viscous and pressure forces on the IB. If tc_in=None, t_case must also contain the projection of the flow field solution onto the IB.
    Usage: loads(t_case, tc_in, tc2_in, wall_out, alpha, beta, gradP, order, Sref, famZones)"""
    if tc_in is not None:
        if isinstance(tc_in, str):
            tc = C.convertFile2PyTree(tc_in)
        else: tc = tc_in
    else: tc = None

    if tc2_in is not None:
        if isinstance(tc2_in, str):
            tc2 = C.convertFile2PyTree(tc2_in)
        else: tc2 = tc2_in
    else: tc2 = None

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    if Sref is None:
        C._initVars(tb, '__ONE__',1.)
        Sref = P.integ(tb, '__ONE__')[0]; print(Sref)
        C._rmVars(tb, ['__ONE__', 'centers:vol'])

    #====================================
    # Wall pressure correction
    #====================================
    if gradP:
        # add gradP fields in tc if necessary
        if tc is not None:
            for z in Internal.getZones(tc):
                _add_gradxi_P(z)

            if order < 2:
                tc = extractPressureHO(tc)
            else:
                tc = extractPressureHO2(tc)

        # add gradP fields in tc2 if necessary
        if tc2 is not None: 
            
            for z in Internal.getZones(tc2):
                _add_gradxi_P(z)
                
            if order < 2:
                tc2 = extractPressureHO(tc2)
            else:
                tc2 = extractPressureHO2(tc2)
            

    #====================================
    # Extraction des grandeurs a la paroi
    #====================================
    if tc is None: 
        zw = Internal.getZones(tb)
        zw = T.join(zw)
    else:
        zw = TIBM.extractIBMWallFields(tc, tb=tb, famZones=famZones)

    #====================================
    # Extract pressure info from tc2 to tc
    #====================================
    if tc2 is not None:
        zw2 = TIBM.extractIBMWallFields(tc2, tb=tb, famZones=famZones, front=1)
        zones_zw  = []
        zones_zw2 = []
        for zone in Internal.getZones(zw): zones_zw.append(zone[0])
        for zone in Internal.getZones(zw2): zones_zw2.append(zone[0])
        nbZones = len(zones_zw)

        for i in range(nbZones): # for multi corps
            szw  = Internal.getNodeFromName(zw, zones_zw[i])
            szw2 = Internal.getNodeFromName(zw2, zones_zw2[i])

            Internal.getNodeFromName(szw, 'Pressure')[1] = Internal.getNodeFromName(szw2, 'Pressure')[1]
            Internal.getNodeFromName(szw, 'Density')[1]  = Internal.getNodeFromName(szw2, 'Density')[1]

            Internal.getNodeFromName(szw, 'gradxPressure')[1] = Internal.getNodeFromName(szw2, 'gradxPressure')[1]
            Internal.getNodeFromName(szw, 'gradyPressure')[1] = Internal.getNodeFromName(szw2, 'gradyPressure')[1]
            Internal.getNodeFromName(szw, 'gradzPressure')[1] = Internal.getNodeFromName(szw2, 'gradzPressure')[1]

    dimPb = Internal.getValue(Internal.getNodeFromName(tb, 'EquationDimension'))

    if dimPb == 2: T._addkplane(zw)

    zw = C.convertArray2Tetra(zw)
    zw = T.reorderAll(zw, 1)

    ts = C.newPyTree(['SKIN']); ts[2][1][2]=zw[2][1][2]

    #==============================
    # Reference state
    #==============================
    RefState = Internal.getNodeFromType(tb,'ReferenceState_t')
    ts[2][1][2].append(RefState)
    ts = loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose=True)

    if dimPb == 2: # reextrait en 2D
        ts = P.isoSurfMC(ts, "CoordinateZ", 0.)
        nodes = Internal.getNodesFromName(ts, 'CoordinateX')
        xmin = numpy.min(nodes[0][1])
        xmax = numpy.max(nodes[0][1])
        dxi = 1./(xmax-xmin)
        C._initVars(ts, 'xsc=({CoordinateX}-%g)*%g'%(xmin, dxi))

    if isinstance(wall_out, str): C.convertPyTree2File(ts, wall_out)
    return ts


#=============================================================================
# Post forces 
# IN: tb: geometry file with solution projected onto it
# IN: Sref: Reference surface area
# OUT: wall_out ou None: fichier pour sortie des efforts sur la paroi aux centres
# IN: alpha: angle pour les efforts
# IN: beta: angle pour les efforts
#==============================================================================
def unsteadyLoads(tb, Sref=None, alpha=0., beta=0.):
    """Computes the viscous and pressure forces on the IB during the computation of the solution. 
    Usage: unsteadyLoads(tb, Sref, alpha, beta)"""
    tp = Internal.copyRef(tb)
    _unsteadyLoads(tp, Sref=Sref, alpha=alpha, beta=beta)
    return tp
 

def _unsteadyLoads(tb, Sref=None, alpha=0., beta=0.):
    """Computes the viscous and pressure forces on the IB during the computation of the solution. 
    Usage: _unsteadyLoads(tb, Sref, alpha, beta)"""
    zones = KCOMM.allgather(Internal.getZones(tb))
    ts    = Distributed.setZonesInTree(tb, zones)
    dimPb = Internal.getValue(Internal.getNodeFromName(ts, 'EquationDimension'))
    return _loads0(ts, Sref=Sref, alpha=alpha, beta=beta, dimPb=dimPb, verbose =False)


#=============================================================================
# Post - General
# IN: t_case: geometry file name or tree
# IN: t_in: result file name or tree
# IN: tc_in: connectivity file name or tree
# OUT: t_out ou None: output file name - values at nodes
# OUT: wall_out ou None: output file name - wall values
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out):
    if isinstance(t_in, str): t = C.convertFile2PyTree(t_in)
    else: t = t_in
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    #=============================
    # Deleting unnecessary fields
    #=============================
    vars =['centers:TurbulentDistance',
           'centers:Density_M1'       , 'centers:Temperature_M1',
           'centers:VelocityX_M1'     , 'centers:VelocityY_M1'  , 'centers:VelocityZ_M1',
           'centers:Density_P1'       , 'centers:Temperature_P1',
           'centers:VelocityX_P1'     , 'centers:VelocityY_P1'  , 'centers:VelocityZ_P1']
    C._rmVars(t, vars)

    #=============================
    # Connectivity tree
    #=============================
    if isinstance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    Internal._rmNodesByName(tc, 'GridCoordinates')

    #==========================================================
    # Compute Cp, Cf, ... on the geometry surface (interpolation)
    #==========================================================
    tb = C.convertArray2Tetra(tb)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
    model = Internal.getValue(model)

    if model == 'Euler': bcType = 0 # slip
    elif model =='NSLaminar': bcType = 1 # noslip
    else: bcType = 3 # Musker

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
     ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,Mus, Cs, Ts, Pr] = C.getState(tb)

    varType = 2 # IBM updated variables (rho,u,t)
    varsIBC = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    vars    = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    

    if model != 'Euler':
        vars += ['ViscosityEddy']
        if model == 'NSTurbulent':
            vars += ['TurbulentSANuTilde']
            varsIBC += ['TurbulentSANuTilde']
            varType = 21

    for z in Internal.getNodesFromType2(t, "Zone_t"):
        zc = Internal.getNodeFromName(tc, z[0])
        for v in varsIBC: C._cpVars(z, 'centers:'+v, zc, v)
        
    X._setInterpTransfers(t, tc, variables=vars,
                          variablesIBC=varsIBC, bcType=bcType,
                          varType=varType, storage=1,
                          Gamma=Gamma, Cv=cvInf, MuS=Mus, 
                          Cs=Cs, Ts=Ts)
    
    zw = TIBM.extractIBMWallFields(tc, tb=tb)
    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':
        C._initVars(zw,'{Cf}=2*%f*{Density}*{utau}**2*%f'%(RoInf,RoUInf2I))

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)

    #================================
    # For 2D, extract a single k plane
    #================================
    if dimPb == 2:
        t = T.subzone(t, (1,1,1), (-1,-1,1))
        C._initVars(t, 'CoordinateZ', 0.) # forced

    #=================================
    # Calc. mu_t/mu in the flow field
    #=================================
    if model != 'Euler':
        betas = Mus*(Ts+Cs)/(Ts**(3./2.))
        C._initVars(t, '{centers:ViscosityMolecular} = %20.16g*sqrt({centers:Temperature})/(1.+%20.16g/{centers:Temperature})'%(betas,Cs))
        C._initVars(t, '{centers:mutsmu}=({centers:ViscosityEddy})/({centers:ViscosityMolecular})-1.')

    #======================================
    # Output of flow solution at cell nodes
    #======================================
    vars = ['centers:Density','centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature','centers:ViscosityEddy',
            'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    #vars = ['centers:Density','centers:VelocityX', 'centers:Temperature','centers:ViscosityEddy',
    #        'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    t = C.center2Node(t, vars)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw





#=============================================================================
# Post-processing pressure extrapolation (at wall) 1st order
# IN: connectivity tree
# OUT: same as input
#=============================================================================
def extractPressureHO(tc):
    """1st order extrapolation of the pressure at the immersed boundary (IB).
    Usage: extractPressureHO(tc)"""
    tp = Internal.copyRef(tc)
    _extractPressureHO(tp)
    return tp


def _extractPressureHO(tc):
    """1st order extrapolation of the pressure at the immersed boundary (IB).
    Usage: _extractPressureHO(tc)"""
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):
                gradxPressure = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                CoordinateX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                Pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density  = Internal.getNodeFromName(zsr, 'Density')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm
                    ny = ny/norm
                    nz = nz/norm
                    
                    nGradP = nx*gradxPressure[i] + ny*gradyPressure[i] + nz*gradzPressure[i]
                    
                    bx   = CoordinateX_PI[i] - CoordinateX[i]
                    by   = CoordinateY_PI[i] - CoordinateY[i]
                    bz   = CoordinateZ_PI[i] - CoordinateZ[i]
                    beta = math.sqrt(bx*bx + by*by + bz*bz)
                    
                    Density[i]  = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta)
                    Pressure[i] = Pressure[i] - nGradP*beta
    return None

#=============================================================================
# Post-processing pressure extrapolation (at wall) 2nd order
# IN: connectivity tree
# OUT: same as input
#=============================================================================    
def extractPressureHO2(tc):
    """2nd order extrapolation of the pressure at the immersed boundary (IB).
    Usage: extractPressureHO2(tc)"""
    tp = Internal.copyRef(tc)
    _extractPressureHO2(tp)
    return tp


def _extractPressureHO2(tc):
    """2nd order extrapolation of the pressure at the immersed boundary (IB).
    Usage: _extractPressureHO2(tc)"""
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):
                CoordinateX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                gradxPressure  = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure  = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure  = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                gradxxPressure = Internal.getNodeFromName(zsr, 'gradxxPressure')[1]
                gradxyPressure = Internal.getNodeFromName(zsr, 'gradxyPressure')[1]
                gradxzPressure = Internal.getNodeFromName(zsr, 'gradxzPressure')[1]

                gradyxPressure = Internal.getNodeFromName(zsr, 'gradyxPressure')[1]
                gradyyPressure = Internal.getNodeFromName(zsr, 'gradyyPressure')[1]
                gradyzPressure = Internal.getNodeFromName(zsr, 'gradyzPressure')[1]

                gradzxPressure = Internal.getNodeFromName(zsr, 'gradzxPressure')[1]
                gradzyPressure = Internal.getNodeFromName(zsr, 'gradzyPressure')[1]
                gradzzPressure = Internal.getNodeFromName(zsr, 'gradzzPressure')[1]

                Pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density  = Internal.getNodeFromName(zsr, 'Density')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm
                    ny = ny/norm
                    nz = nz/norm

                    nGradP   = nx*gradxPressure[i] + ny*gradyPressure[i] + nz*gradzPressure[i]

                    nnxGradP = nx*gradxxPressure[i] + ny*gradxyPressure[i] + nz*gradxzPressure[i]
                    nnyGradP = nx*gradyxPressure[i] + ny*gradyyPressure[i] + nz*gradyzPressure[i]
                    nnzGradP = nx*gradzxPressure[i] + ny*gradzyPressure[i] + nz*gradzzPressure[i]
                    nnGradP  = nx*nnxGradP + ny*nnyGradP + nz*nnzGradP

                    bx = CoordinateX_PI[i] - CoordinateX[i]
                    by = CoordinateY_PI[i] - CoordinateY[i]
                    bz = CoordinateZ_PI[i] - CoordinateZ[i]
                    beta = math.sqrt(bx*bx + by*by + bz*bz)

                    Density[i]  = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2)
                    Pressure[i] = Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2

    return None


#=============================================================================
# Compute convective terms (TBLE FULL)
#=============================================================================
def extractConvectiveTerms(tc):
    """Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc.
    Usage: extractConvectiveTerms(tc)"""
    tp = Internal.copyRef(tc)
    _extractConvectiveTerms(tp)
    return tp


def _extractConvectiveTerms(tc):
    """Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc.
    Usage: _extractConvectiveTerms(tc)"""
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):

                CoordinateX    = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                CoordinateY    = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                CoordinateZ    = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                CoordinateX_PC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                CoordinateY_PC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                CoordinateZ_PC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                CoordinateX_PI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                CoordinateY_PI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                CoordinateZ_PI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                Pressure       = Internal.getNodeFromName(zsr, 'Pressure')[1]
                Density        = Internal.getNodeFromName(zsr, 'Density')[1]

                gradxPressure  = Internal.getNodeFromName(zsr, 'gradxPressure')[1]
                gradyPressure  = Internal.getNodeFromName(zsr, 'gradyPressure')[1]
                gradzPressure  = Internal.getNodeFromName(zsr, 'gradzPressure')[1]

                gradxVelocityX = Internal.getNodeFromName(zsr, 'gradxVelocityX')[1]
                gradyVelocityX = Internal.getNodeFromName(zsr, 'gradyVelocityX')[1]
                gradzVelocityX = Internal.getNodeFromName(zsr, 'gradzVelocityX')[1]

                gradxVelocityY = Internal.getNodeFromName(zsr, 'gradxVelocityY')[1]
                gradyVelocityY = Internal.getNodeFromName(zsr, 'gradyVelocityY')[1]
                gradzVelocityY = Internal.getNodeFromName(zsr, 'gradzVelocityY')[1]

                gradxVelocityZ = Internal.getNodeFromName(zsr, 'gradxVelocityZ')[1]
                gradyVelocityZ = Internal.getNodeFromName(zsr, 'gradyVelocityZ')[1]
                gradzVelocityZ = Internal.getNodeFromName(zsr, 'gradzVelocityZ')[1]

                VelocityX      = Internal.getNodeFromName(zsr, 'VelocityX')[1]
                VelocityY      = Internal.getNodeFromName(zsr, 'VelocityY')[1]
                VelocityZ      = Internal.getNodeFromName(zsr, 'VelocityZ')[1]

                nIBC = numpy.shape(CoordinateX)[0]
                conv1  = numpy.zeros((nIBC),numpy.float64)
                conv2  = numpy.zeros((nIBC),numpy.float64)

                for i in range(nIBC):
                    nx = CoordinateX_PC[i] - CoordinateX[i]
                    ny = CoordinateY_PC[i] - CoordinateY[i]
                    nz = CoordinateZ_PC[i] - CoordinateZ[i]
                    norm_n = math.sqrt(nx*nx + ny*ny + nz*nz)
                    nx = nx/norm_n
                    ny = ny/norm_n
                    nz = nz/norm_n

                    uscaln = (VelocityX[i]*nx + VelocityY[i]*ny + VelocityZ[i]*nz)
                    tx = VelocityX[i] - uscaln*nx
                    ty = VelocityY[i] - uscaln*ny
                    tz = VelocityZ[i] - uscaln*nz
                    norm_t = math.sqrt(tx*tx + ty*ty + tz*tz)
                                        
                    tx = tx/norm_t
                    ty = ty/norm_t
                    tz = tz/norm_t

                    tgradU =  (tx*gradxVelocityX[i] + ty*gradyVelocityX[i] + tz*gradzVelocityX[i])*tx
                    tgradU += (tx*gradxVelocityY[i] + ty*gradyVelocityY[i] + tz*gradzVelocityY[i])*ty
                    tgradU += (tx*gradxVelocityZ[i] + ty*gradyVelocityZ[i] + tz*gradzVelocityZ[i])*tz

                    ngradU =  (nx*gradxVelocityX[i] + ny*gradyVelocityX[i] + nz*gradzVelocityX[i])*tx
                    ngradU += (nx*gradxVelocityY[i] + ny*gradyVelocityY[i] + nz*gradzVelocityY[i])*ty
                    ngradU += (nx*gradxVelocityZ[i] + ny*gradyVelocityZ[i] + nz*gradzVelocityZ[i])*tz

                    norm_n = math.sqrt((uscaln*nx)**2 + (uscaln*ny)**2 + (uscaln*nz)**2)

                    conv1[i] = norm_t*tgradU
                    conv2[i] = norm_n*ngradU

                zsr[2].append(['conv1'  , conv1  , [], 'DataArray_t'])
                zsr[2].append(['conv2'  , conv2  , [], 'DataArray_t'])

    return None


#==========================================================================================
# IN:  ts            : skin (TRI zones) distributed already (partial tree here)
# IN:  tc            : connectivity tree
# OUT: tl            : NODE-type zones of IBM points to be projected locally on ts
# OUT: graphWPOST    : graph of coms between tc and tl
# OUT: interDictWPOST: dictionary of intersection domains
#==========================================================================================
def _prepareSkinReconstruction(ts, tc):
    tBBs=Cmpi.createBBoxTree(ts)
    procDictBBs = Cmpi.getProcDict(tBBs)

    basename=Internal.getName(Internal.getBases(ts)[0])
    tl = C.newPyTree([basename])
    hmin = 0.
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")                  
        GCnode  = Internal.getNodeFromType(zc,"GridCoordinates_t")
        XN      = Internal.getNodeFromName(GCnode,'CoordinateX')

        for IBCD in allIBCD:
            if XN is not None:
                if XN[1].shape[0]>1:
                    hx = C.getValue(zc,'CoordinateX',1)-C.getValue(zc,'CoordinateX',0)
                    hy = C.getValue(zc,'CoordinateY',1)-C.getValue(zc,'CoordinateY',0)
                    hloc = max(abs(hx),abs(hy))
                    hmin = max(hloc,hmin)

            zname = Internal.getValue(IBCD)
            XW    = Internal.getNodeFromName(IBCD,'CoordinateX_PW')[1]
            YW    = Internal.getNodeFromName(IBCD,'CoordinateY_PW')[1]
            ZW    = Internal.getNodeFromName(IBCD,'CoordinateZ_PW')[1]

            zsize = numpy.empty((1,3), numpy.int32, order='F')
            zsize[0,0] = XW.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
            z  = Internal.newZone(name='IBW_Wall_%s_%s'%(zc[0],zname),zsize=zsize,ztype='Unstructured')
            gc = Internal.newGridCoordinates(parent=z)
            coordx = ['CoordinateX',XW,[],'DataArray_t']
            coordy = ['CoordinateY',YW,[],'DataArray_t']
            coordz = ['CoordinateZ',ZW,[],'DataArray_t']
            gc[2] = [coordx,coordy,coordz]
            n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
            Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
            Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
            FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,gridLocation='Vertex', parent=z)
            pressNP = []; utauNP = []; yplusNP = []; densNP = []
            vxNP = []; vyNP = []; vzNP = []

            PW = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
            if PW is not None: pressNP.append(PW[1])
            RHOW = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
            if RHOW is not None: densNP.append(RHOW[1])
            UTAUW = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
            if UTAUW is not None: utauNP.append(UTAUW[1])
            YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
            if YPLUSW is not None: yplusNP.append(YPLUSW[1])

            VXW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
            if VXW is not None: vxNP.append(VXW[1])
            VYW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
            if VYW is not None: vyNP.append(VYW[1])
            VZW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
            if VZW is not None: vzNP.append(VZW[1])

            FSN[2].append([XOD.__PRESSURE__,pressNP[0], [],'DataArray_t'])
            FSN[2].append([XOD.__DENSITY__,densNP[0], [],'DataArray_t'])
            if utauNP != []:
                FSN[2].append([XOD.__UTAU__,utauNP[0], [],'DataArray_t'])
            if yplusNP != []:
                FSN[2].append([XOD.__YPLUS__,yplusNP[0], [],'DataArray_t'])

            if vxNP != []:
                FSN[2].append([XOD.__VELOCITYX__,vxNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYY__,vyNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYZ__,vzNP[0], [],'DataArray_t'])

            Cmpi._setProc(z,Cmpi.rank)          
            tl[2][1][2].append(z)

    tlBB=Cmpi.createBBoxTree(tl, tol=hmin)
    procDictWPOST = Cmpi.getProcDict(tlBB)
    interDictWPOST = X.getIntersectingDomains(tlBB, tBBs)
    graphWPOST = Cmpi.computeGraph(tlBB, type='bbox3',intersectionsDict=interDictWPOST,
                                   procDict=procDictWPOST, procDict2=procDictBBs, t2=tBBs)
 
    RefStateNode = Internal.getNodeFromName(ts,'ReferenceState')
    tl[2][1][2].append(RefStateNode)
    FES =  Internal.getNodeFromName(ts,'FlowEquationSet')
    tl[2][1][2].append(FES)

    C._initVars(ts,XOD.__PRESSURE__,0.)
    C._initVars(ts,XOD.__DENSITY__,0.)
    C._initVars(ts,XOD.__VELOCITYX__,0.)
    C._initVars(ts,XOD.__VELOCITYY__,0.)
    C._initVars(ts,XOD.__VELOCITYZ__,0.)    
    if Internal.getValue(Internal.getNodeFromType1(FES,'GoverningEquations_t'))!= 'Euler':
        C._initVars(ts,XOD.__UTAU__,0.)
        C._initVars(ts,XOD.__YPLUS__,0.)
    
    return tl, graphWPOST, interDictWPOST


# Projected solution onto ts.
# ts: skin (TRI zones) distributed already (partial tree here)
def _computeSkinVariables(ts, tc, tl, graphWPOST, interDictWPOST):
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            PW = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
            RHOW = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
            UTAUW = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
            YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
            VXW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
            VYW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
            VZW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
            
            zname = Internal.getValue(IBCD)
            znamepostw = 'IBW_Wall_%s_%s'%(zc[0],zname)
            zpostw = Internal.getNodeFromName(tl,znamepostw)
            FSP = Internal.getNodeFromType(zpostw,'FlowSolution_t')
            PW2 = Internal.getNodeFromName1(FSP,XOD.__PRESSURE__)
            RHOW2 = Internal.getNodeFromName1(FSP,XOD.__DENSITY__)
            PW2[1]=PW[1]; RHOW2[1]=RHOW[1]

            UTAUW2 = Internal.getNodeFromName1(FSP,XOD.__UTAU__)
            if UTAUW2 is not None:
                YPLUSW2 = Internal.getNodeFromName1(FSP, XOD.__YPLUS__)
                UTAUW2[1]=UTAUW[1]; YPLUSW2[1]=YPLUSW[1]
            VXW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYX__)     
            if VXW2 is not None:
                VYW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYY__)
                VZW2 = Internal.getNodeFromName1(FSP, XOD.__VELOCITYZ__)
                VXW2[1]=VXW[1]
                VYW2[1]=VYW[1]
                VZW2[1]=VZW[1]


    tdl = Cmpi.addXZones(tl, graphWPOST)
    tdl = Cmpi.convert2PartialTree(tdl)
    for nobs in range(len(ts[2])):
        if Internal.getType(ts[2][nobs])=='CGNSBase_t':
            for nozs in range(len(ts[2][nobs][2])):
                zs = ts[2][nobs][2][nozs]
                if Internal.getType(zs)=='Zone_t':
                    cloud = []
                    for zl in Internal.getZones(tdl):
                        if zl != [] and zl is not None and zs[0] in interDictWPOST[zl[0]]:
                            zl = C.convertArray2Node(zl)
                            cloud.append(zl)
    
                    if cloud != []:
                        cloud = T.join(cloud)
                        
                        ts[2][nobs][2][nozs] = P.projectCloudSolution(cloud, zs, dim=3)
                        
    return None


## IMPORTANT NOTE !!
## The functions below will become decrepit after Jan. 1 2023
#====================================================================================
def prepareWallReconstruction(tw, tc):
    # MLS interpolation order
    LSorder = 2

    dimPb = Internal.getNodeFromName(tw, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    rank = Cmpi.rank
    # GLOBAL
    # pour eviter les trous lies au split des zones IBM en nuages de points
    tolBB = 0.
    for snear in Internal.getNodesFromName(tw,"snear"):
        snearval = Internal.getValue(snear)
        tolBB = max(tolBB,2*snearval)

    tcw = TIBM.createIBMWZones(tc,variables=[])

    nzones = len(Internal.getZones(tcw))

    dimPb = Internal.getNodeFromName(tw, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    if dimPb ==2: C._initVars(tcw,"CoordinateZ",0.0)
    # 
    # Create the bbtree 
    tcw_bb = C.newPyTree(["Base"])
    zonesl = []
    for base in Internal.getBases(tcw):
        for z in Internal.getZones(base):
            zbb = G.BB(z)

            # Clean up (zoneSubRegion)
            Internal._rmNodesFromType(zbb, 'ZoneSubRegion_t')
            valxm = C.getValue(zbb,'CoordinateX',0)-tolBB
            valxp = C.getValue(zbb,'CoordinateX',1)-tolBB
            C.setValue(zbb,'CoordinateX',0, valxm)
            C.setValue(zbb,'CoordinateX',1, valxp)
            valxm = C.getValue(zbb,'CoordinateY',0)-tolBB
            valxp = C.getValue(zbb,'CoordinateY',1)-tolBB
            C.setValue(zbb,'CoordinateY',0, valxm)
            C.setValue(zbb,'CoordinateY',1, valxp)
            if dimPb == 3:
                valxm = C.getValue(zbb,'CoordinateZ',0)-tolBB
                valxp = C.getValue(zbb,'CoordinateZ',1)-tolBB
                C.setValue(zbb,'CoordinateZ',0, valxm)
                C.setValue(zbb,'CoordinateZ',1, valxp)
            Cmpi._setProc(zbb,rank)
            zonesl.append(zbb)    

    zonesBB = Cmpi.allgatherZones(zonesl)
    tcw_bb[2][1][2] += zonesBB
    graph={}
    for zw in Internal.getZones(tw):
        zwBB = G.BB(zw)
        for zdBB in Internal.getZones(tcw_bb):
            if G.bboxIntersection(zwBB,zdBB,isBB=True):
                popp = Cmpi.getProc(zdBB)
                zdname = zdBB[0]
                Distributed.updateGraph__(graph, popp, rank, zdname)

    allGraph = Cmpi.KCOMM.allgather(graph)
    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j])) 

    Cmpi._addXZones(tcw, graph, subr=False)
    interDict = X.getIntersectingDomains(tw,tcw_bb)
    procDict = Cmpi.getProcDict(tcw)

    graphX = {}
    for zw in Internal.getZones(tw):
        zwname = Internal.getName(zw)
        for zdname in interDict[zwname]:
            zdbb = Internal.getNodeFromName2(tcw_bb,zdname)
            popp = Cmpi.getProc(zdbb)
            Distributed.updateGraph__(graphX, rank, popp, zdname)

    allGraph = Cmpi.KCOMM.allgather(graphX)
    graphX = {}
    for i in allGraph:
        for k in i:
            if not k in graphX: graphX[k] = {}
            for j in i[k]:
                if not j in graphX[k]: graphX[k][j] = []
                graphX[k][j] += i[k][j]
                graphX[k][j] = list(set(graphX[k][j])) 
    del tcw_bb 

    EXTRAP = numpy.array([],numpy.int32)
    VOL = numpy.array([],numpy.float64)
    ORPHAN = numpy.array([],numpy.float64)
    datas = {}

    for zw in Internal.getZones(tw):
        zwname = Internal.getName(zw)
        dnrZones=[]
        dnrZoneNames=[]
        for zdname in interDict[zwname]:
            zd = Internal.getNodeFromName2(tcw,zdname)     
            dnrZones.append(zd)
            dnrZoneNames.append(zdname)

        coordsD = C.getFields(Internal.__GridCoordinates__, dnrZones, api=1)
        coordsR = C.getFields(Internal.__GridCoordinates__, zw, api=1)[0]        
        if coordsR[1].shape[1]>0:
            coordsD = Converter.convertArray2Node(coordsD)
            ret = connector.setInterpData_IBMWall(coordsD, coordsR, dimPb, LSorder)

            allPLR = ret[0]; allPLD = ret[1]; allCOEFS=ret[3]; allITYPE=ret[2]

            nozd = 0
            for zd in dnrZones:
                zdname = zd[0]
                XOD._createInterpRegion__(zd, zw[0], allPLD[nozd], allPLR[nozd], allCOEFS[nozd], allITYPE[nozd], \
                                          VOL, EXTRAP, ORPHAN, \
                                          tag='Donor', loc='nodes',itype='chimera', prefix='IDW_')
                nozd+=1

                destProc = procDict[zdname]

                IDs = []
                for i in zd[2]:
                    if i[0][0:3] == 'IDW':
                        if Internal.getValue(i)==zwname: IDs.append(i)

                if IDs != []:
                    if destProc == rank:
                        zD = Internal.getNodeFromName2(tcw, zdname)
                        zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []
    for i in range(Cmpi.size):
        if i not in datas: datas[i]=[]
    Cmpi._rmXZones(tcw)
    destDatas = Cmpi.sendRecv(datas, graphX)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tcw, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}
    return tcw
 

# reconstruit les champs parietaux a partir des infos stockees dans tw et les champs de tc
def _computeWallReconstruction(tw, tcw, tc, procDictR=None, procDictD=None, graph=None, 
                               variables=['Pressure','Density','utau','yplus']):
    if procDictR is None:
        procDictR={}
        for z in Internal.getZones(tw): procDictR[z[0]]=0
    if procDictD is None:
        procDictD={}
        for z in Internal.getZones(tcw): procDictD[z[0]]=0

    if graph is None: 
        graph = Cmpi.computeGraph(tcw, type='POST',procDict=procDictD, procDict2=procDictR, t2=tw)

    cellNVariable=''; varType=1; compact=0

    datas={}
    # interpolation from a cloud of points
    for zc_w in Internal.getZones(tcw):
        infos = []
        zcw_name = Internal.getName(zc_w)
        zcw_namel = zcw_name.split('#')
        zc_orig_name = zcw_namel[0]; zsrname = zcw_namel[1]

        FSN = Internal.getNodeFromName2(zc_w,Internal.__FlowSolutionNodes__)
        if FSN is None:
            FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                           gridLocation='Vertex', parent=zc_w)

        # put the fields in corresponding zcw zone 
        zc_orig = Internal.getNodeFromName2(tc,zc_orig_name)
        zsr_orig = Internal.getNodeFromName2(zc_orig,zsrname)
        for varname in variables:
            varnode = Internal.getNodeFromName(zsr_orig,varname)
            FSN[2].append(varnode)

        # search for IDW_zrname node
        for zsr in Internal.getNodesFromType1(zc_w, "ZoneSubRegion_t"):
            sname = zsr[0][0:3]
            if sname=='IDW' and variables is not None:
                dname = Internal.getValue(zsr)
                idn = Internal.getNodeFromName1(zsr, 'InterpolantsDonor')
                if idn is not None: # la subRegion decrit des interpolations
                    zoneRole = Internal.getNodeFromName2(zsr, 'ZoneRole')
                    zoneRole = Internal.getValue(zoneRole)
                    if zoneRole == 'Donor':
                        location = Internal.getNodeFromName1(zsr, 'GridLocation') # localisation des donnees des receveurs

                        if location is not None: 
                            location = Internal.getValue(location)
                            if location == 'CellCenter': loc = 'centers'
                            else: loc = 'nodes'
    
                        else: loc = 'nodes'
    
                        Coefs     = idn[1]
                        DonorType = Internal.getNodeFromName1(zsr,'InterpolantsType')[1]
                        ListDonor = Internal.getNodeFromName1(zsr,'PointList')[1]
                        ListRcv   = Internal.getNodeFromName1(zsr,'PointListDonor')[1]

                        arrayT = connector._setInterpTransfersD(zc_w, variables, ListDonor, DonorType, Coefs, varType, compact,                                                                
                                                                cellNVariable,
                                                                Internal.__GridCoordinates__, 
                                                                Internal.__FlowSolutionNodes__, 
                                                                Internal.__FlowSolutionCenters__)   
                        infos.append([dname,arrayT,ListRcv,loc])
        for n in infos:
            rcvName = n[0]
            proc = procDictR[rcvName]
            if proc == Cmpi.rank:
                fields = n[1]
                if fields != []:
                    listIndices = n[2]
                    zr = Internal.getNodeFromName2(tw,rcvName)     
                    C._updatePartialFields(zr, [fields], [listIndices], loc=n[3])
            else:
                rcvNode = procDictR[rcvName]
                if rcvNode not in datas: datas[rcvNode] = [n]
                else: datas[rcvNode] += [n] 

    #send numpys using graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # put contribution of donor zone to the interpolated fields in the receptor zone
    for i in rcvDatas:
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodeFromName2(tw, rcvName)
                C._updatePartialFields(z,[field], [listIndices], loc=n[3])
    return None


#====================================================================================
