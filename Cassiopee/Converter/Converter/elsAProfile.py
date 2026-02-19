# Module for exporting to elsA
from . import PyTree as C
from . import Internal
import Connector.PyTree as X
from . import Converter
import numpy
import math

__CONSERVATIVE__ = ["Density", "MomentumX", "MomentumY", "MomentumZ", "EnergyStagnationDensity"]
__TURBULENT__ = ["TurbulentEnergyKineticDensity", "TurbulentDissipationDensity", "TurbulentSANuTildeDensity"]
__COMMONS__ = ["Pressure", "Mach", "Temperature"]
__COMMONSNS__ = ["Viscosity_EddyMolecularRatio"]
__WALLDISTANCE__ = ["TurbulentDistance","TurbulentDistanceIndex"]
__XYZ__ = ['CoordinateX', 'CoordinateY', 'CoordinateZ']

__FAMOVERLAPBC__ = 'F_OV_'
__FAMOVERLAPDDBC__ = 'F_OVDD_'
__CHIMGROUPNAME__ = 'ChimGroup_' # reference to donors for suffixed base name
__CHIMGROUPNAMEDD__ = 'ChimGroupDD_' # reference to donors for suffixed base name

##############################################################################

# Traduction : keyword elsA -> keyword CGNS
keyselsA2CGNS = {\
    'config'         :'EquationDimension'             , \
    'fluid'          :'GasModel'                      , \
    ''               :'GasModelType'                  , \
    'pg'             :'Ideal'                         , \
    'gamma'          :'SpecificHeatRatio'             , \
    'cv'             :'SpecificHeatVolume'            , \
    '_'              :'IdealGasConstant'              , \
    '__'             :'ThermalConductivityModel'      , \
    '___'            :'TurbulenceClosure'             , \
    'prandtl'        :'ConstantPrandtl'               , \
    'prandtltb'      :'PrandtlTurbulent'              , \
    'phymod'         :"GoverningEquations"            , \
    'turbmod'        :'TurbulenceModel'               , \
    'euler'          :'Euler'                         , \
    'nstur'          :'NSTurbulent'                   , \
    'nslam'          :'NSLaminar'                     , \
    'spalart'        :'OneEquation_SpalartAllmaras'   , \
    'kepsjl'         :'TwoEquation_JonesLaunder'      , \
    'komega_menter'  :'TwoEquation_MenterSST'         , \
    'komega_wilcox'  :'TwoEquation_Wilcox'            , \
    'komega_kok'     :'UserDefined'                   , \
    'smith'          :'UserDefined'                   , \
    'visclaw'        :'ViscosityModel'                , \
    'sutherland'     :'Sutherland'                    , \
    'suth_const'     :'SutherlandLawConstant'         , \
    'suth_muref'     :'ViscosityMolecularReference'   , \
    'suth_tref'      :'TemperatureReference'          , \
    'walladia'       :'BCWall'                        , \
    'wallslip'       :'BCWallInviscid'                , \
    'cell'           :'CellCenter'                    , \
    'node'           :'Vertex'                        , \
    'x'              :'CoordinateX'                   , \
    'y'              :'CoordinateY'                   , \
    'z'              :'CoordinateZ'                   , \
    'X'              :'CoordinateX'                   , \
    'Y'              :'CoordinateY'                   , \
    'Z'              :'CoordinateZ'                   , \
    'ro'             :'Density'                       , \
    'rou'            :'MomentumX'                     , \
    'rov'            :'MomentumY'                     , \
    'row'            :'MomentumZ'                     , \
    'rovx'           :'MomentumX'                     , \
    'rovy'           :'MomentumY'                     , \
    'rovz'           :'MomentumZ'                     , \
    'roe'            :'EnergyStagnationDensity'       , \
    'roE'            :'EnergyStagnationDensity'       , \
    'rok'            :'TurbulentEnergyKineticDensity' , \
    'roeps'          :'TurbulentDissipationDensity'   , \
    'ronutilde'      :'TurbulentSANuTildeDensity'     , \
    'mach'           :'Mach'                          , \
    'psta'           :'Pressure'                      , \
    'tsta'           :'Temperature'                   , \
    'viscrapp'       :'Viscosity_EddyMolecularRatio'  , \
    'walldistance'   :'TurbulentDistance'             , \
    'wallglobalindex':'TurbulentDistanceIndex'        , \
}

# Traduction : keyword CGNS -> keyword elsA
# keysCGNS2elsA = dict((v, k) for k, v in keyselsA2CGNS.iteritems())
keysCGNS2elsA={
    'EquationDimension'             :'config'       , \
    'GasModel'                      :'fluid'        , \
    'GasModelType'                  :''             , \
    'Ideal'                         :'pg'           , \
    #'CalloricallyPerfect'           :'pg'           , \
    'SpecificHeatRatio'             :'gamma'        , \
    'SpecificHeatVolume'            :'cv'           , \
    'IdealGasConstant'              :'_'            , \
    'ThermalConductivityModel'      :'__'           , \
    'TurbulenceClosure'             :'___'          , \
    'ConstantPrandtl'               :'prandtl'      , \
    'PrandtlTurbulent'              :'prandtltb'    , \
    "GoverningEquations"            :'phymod'       , \
    'TurbulenceModel'               :'turbmod'      , \
    'Euler'                         :'euler'        , \
    'NSTurbulent'                   :'nstur'        , \
    'NSLaminar'                     :'nslam'        , \
    'OneEquation_SpalartAllmaras'   :'spalart'      , \
    'TwoEquation_JonesLaunder'      :'kepsjl'       , \
    'TwoEquation_MenterSST'         :'komega_menter', \
    'TwoEquation_Wilcox'            :'komega_wilcox', \
    'UserDefined'                   :'komega_kok'   , \
    'UserDefined'                   :'smith'        , \
    'ViscosityModel'                :'visclaw'      , \
    'Sutherland'                    :'sutherland'   , \
    'SutherlandLawConstant'         :'suth_const'   , \
    'ViscosityMolecularReference'   :'suth_muref'   , \
    'TemperatureReference'          :'suth_tref'    , \
    'BCWall'                        :'walladia'     , \
    'BCWallInviscid'                :'wallslip'     , \
    'CellCenter'                    :'cell'         , \
    'Vertex'                        :'node'         , \
    'CoordinateX'                   :'x'            , \
    'CoordinateY'                   :'y'            , \
    'CoordinateZ'                   :'z'            , \
    'Density'                       :'ro'           , \
    'MomentumX'                     :'rou'          , \
    'MomentumY'                     :'rov'          , \
    'MomentumZ'                     :'row'          , \
    'EnergyStagnationDensity'       :'roe'          , \
    "TurbulentEnergyKineticDensity" :'rok'          , \
    "TurbulentDissipationDensity"   :'roeps'        , \
    "TurbulentSANuTildeDensity"     :'ronutilde'    , \
    'Mach'                          :'mach'         , \
    'Pressure'                      :'psta'         , \
    'Temperature'                   :'tsta'         , \
    'Viscosity_EddyMolecularRatio'  :'viscrapp'     , \
}

#==============================================================================
# clean=True: delete useless nodes for elsA (Periodic_t)
# the current window to the
# In the CGNS Std (connectMatchPeriodic is compliant) :
# RotationAngle defines the angle from the current interface to the connecting interface
# Translation defines the vector from the current interface to the connecting interface
# In elsA: it is the opposite
#===============================================================================
def adaptPeriodicMatch(t, clean=False):
    """Convert Periodic Match Grid Connectivity (GC) data to be compliant with elsA solver."""
    tp = Internal.copyRef(t)
    _adaptPeriodicMatch(tp, clean=clean)
    return tp

def _adaptPeriodicMatch(t, clean=False):
    jointype = numpy.array([c for c in 'join'], 'c')
    jtopo = numpy.array([c for c in 'periodic'], 'c')
    ptype_rot = numpy.array([c for c in 'rot'], 'c')
    ptype_tra = numpy.array([c for c in 'tra'], 'c')
    jtype = numpy.array([c for c in 'match'], 'c')

    for z in Internal.getZones(t):
        dim = Internal.getZoneDim(z)
        ni = dim[1]; nj = dim[2]; nk = dim[3]; ninj = ni*nj
        connect = Internal.getNodesFromType(z, 'GridConnectivity1to1_t')
        for c in connect:
            donorName = Internal.getValue(c)
            zopp = Internal.getNodeFromName(t, donorName)
            dimopp = Internal.getZoneDim(zopp)
            niopp = dimopp[1]; njopp = dimopp[2]; nkopp = dimopp[3]; ninjopp = niopp*njopp
            periodic = Internal.getNodeFromType2(c, 'Periodic_t')
            if periodic is not None:
                rotAngleNode = Internal.getNodeFromName1(periodic, 'RotationAngle')
                if rotAngleNode is not None:
                    #rotAngle = Internal.getValue(rotAngleNode)
                    rotAngle = numpy.zeros(3,numpy.float64)
                    rotAngleDeg = Internal.getRotationAngleValueInDegrees(rotAngleNode)
                    for i in range(3): rotAngle[i]=rotAngleDeg[i]
                    Internal._rmNodesByNameAndType(rotAngleNode,'DimensionalUnits','DimensionalUnits_t')
                else:
                    rotAngle = numpy.empty(0,numpy.float64)

                rotCenter = Internal.getNodeFromName1(periodic, 'RotationCenter')
                if rotCenter is not None: rotCenter=rotCenter[1]
                else: rotCenter = numpy.zeros(3,numpy.float64)

                translation = Internal.getNodeFromName1(periodic, 'Translation')
                if translation is not None: translation=translation[1]
                else: translation = numpy.empty(0,numpy.float64)

                solverProperty = c[2][len(c[2])-1]
                sp  = [Internal.createNode('type','DataArray_t',value=jointype)]
                sp += [Internal.createNode('jtopo','DataArray_t',value=jtopo)]
                sp += [Internal.createNode('jtype','DataArray_t',value=jtype)]

                if rotAngle.any(): # at least one rotation angle not nul => periodicity by rotation
                    if len(numpy.where(rotAngle)[0]) != 1:
                        print("Warning: adaptPeriodicMatch__: rotation angle must have one non-zero component.")
                        continue
                    axis = numpy.zeros(3,numpy.float64)
                    axis[numpy.where(rotAngle)[0][0]]=1.0
                    angle = -rotAngle[numpy.where(rotAngle)[0][0]] # angle of the periodicity by rotation
                    angle1 = int(round(360./angle))
                    angle2 = 1                                    # number of channels defined by the component which includes the zone. Has to be set by user afterward
                    pointRange = Internal.getNodeFromName1(c,'PointRange')
                    pointRangeDonor = Internal.getNodeFromName1(c,'PointRangeDonor')
                    if pointRange is None or pointRangeDonor is None:
                        print("Warning: adaptPeriodicMatch__: missing PointRange(Donor) for join ",c[0],".")
                        continue

                    [wincur_imin, wincur_imax, wincur_jmin, wincur_jmax, wincur_kmin, wincur_kmax] = Internal.range2Window(pointRange[1])
                    [winopp_imin, winopp_imax, winopp_jmin, winopp_jmax, winopp_kmin, winopp_kmax] = Internal.range2Window(pointRangeDonor[1])
                    index_cur = ninj*(wincur_kmin-1)    + ni*(wincur_jmin-1)    + (wincur_imin-1) # Index of point Pcur
                    index_opp = ninjopp*(winopp_kmin-1) + niopp*(winopp_jmin-1) + (winopp_imin-1) # Index of point Popp
                    # Sign of the periodicity angle wrt the trigonometrical convention
                    if angle > 0.: pangle=1
                    else: pangle=-1

                    sp += [Internal.createNode('ptype','DataArray_t',value=ptype_rot)]
                    sp += [Internal.createNode('pangle','DataArray_t',value=pangle)]
                    # angle and center are managed by node .Solver#Param, child of the zone
                    _addPeriodicDataInSolverParam(z, rotationCenter=rotCenter, rotationAngle=rotAngle, isChimera=False)

                else:
                    if translation.any(): # translation vector not null => periodicity by translation
                        sp += [Internal.createNode('ptype','DataArray_t',value=ptype_tra)]
                        sp += [Internal.createNode('xtran','DataArray_t',value=-translation[0])]
                        sp += [Internal.createNode('ytran','DataArray_t',value=-translation[1])]
                        sp += [Internal.createNode('ztran','DataArray_t',value=-translation[2])]

                Internal.createNode('.Solver#Property','UserDefinedData_t',value=None,children=sp, parent=c)
                if clean:
                    Internal._rmNodesFromType(c, 'Periodic_t')

    return None

#---------------------------------------------------------------------------------------------------------
# add periodic data for grid connectivity in zone in a .Solver#Param node
# direction: for periodic Chimera
#---------------------------------------------------------------------------------------------------------
def _addPeriodicDataInSolverParam(a, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.],
                                  NAzimutalSectors=0, isChimera=False):
    if isChimera: periodicDir=3
    else: periodicDir=0
    if sum(rotationAngle)==0.: return None
    else:
        xc = float(rotationCenter[0])
        yc = float(rotationCenter[1])
        zc = float(rotationCenter[2])

        diraxis = numpy.where(rotationAngle)[0][0]
        axis=[0.,0.,0.]; axis[diraxis] = 1.
        if NAzimutalSectors>0: angle1 = int(NAzimutalSectors)
        else: angle1 = int(round(360./abs(rotationAngle[diraxis]))) # number of angular sectors of the component which includes the zone in radian
        angle2 = 1

        paramNames = ['axis_ang_1','axis_ang_2','axis_pnt_x','axis_pnt_y','axis_pnt_z','axis_vct_x','axis_vct_y','axis_vct_z']
        paramValues = [angle1,angle2,xc,yc,zc,axis[0],axis[1],axis[2]]
        if isChimera: paramNames.append('periodic_dir'); paramValues.append(periodicDir)

        for z in Internal.getZones(a):
            solverParam = Internal.getNodeFromName(z, '.Solver#Param')
            if solverParam is None:
                nodes=[]
                for nop, v in enumerate(paramValues):
                    nodes.append(Internal.createNode(paramNames[nop],'DataArray_t',value=v))
                Internal._createChild(z,".Solver#Param",'UserDefinedData_t',value=None, children=nodes)
            else:
                for nop, v in enumerate(paramValues):
                    node = Internal.getNodeFromName(solverParam, paramNames[nop])
                    if node is None:
                        Internal._createChild(solverParam, paramNames[nop], 'DataArray_t', value=v)
                    else:
                        Internal.setValue(node, v)

    return None

def addPeriodicDataInSolverParam(a, rotationCenter=[0.,0.,0.],
                                 rotationAngle=[1.,0.,0.], NAzimutalSectors=0, isChimera=False):
    """Add periodic data for grid connectivity in zone in a .Solver#Param node."""
    tp = Internal.copyRef(a)
    _addPeriodicDataInSolverParam(tp,rotationCenter=rotationCenter,rotationAngle=rotationAngle,
                                  NAzimutalSectors=NAzimutalSectors, isChimera=isChimera)
    return tp

# -----------------------------------------------------------------------------
def getCGNSkeys(key, verbose=True):
    """Return the CGNS name (if it exists) corresponding to the elsA key.
    Usage: getCGNSkeys(key, verbose=True)"""
    if key in keyselsA2CGNS: return keyselsA2CGNS[key]
    elif key in keyselsA2CGNS.values(): return key
    else:
        if verbose: print('Warning: getCGNSkeys: the given key %s cannot be translated in a CGNS key.'%key)
        return key

# -----------------------------------------------------------------------------
def _addOutput(a, Dict, name='', update=False):
    """Add a node '.Solver#Output'+name to a node a to extract required values from that node a.
    Usage: _addOutput(a, Dict, name='', update=False)"""
    if not isinstance(Dict, dict):
        raise TypeError("_addOutput: second argument is not a Python dictionary.")

    outputName = '.Solver#Output'+name
    outputnode = Internal.getNodeFromName1(a, outputName)

    if outputnode is None:
        nodeso = []
        for each in Dict: nodeso+=[Internal.createNode(each,'DataArray_t',value=Dict[each])]
        a[2].append(Internal.createNode(outputName, 'UserDefinedData_t',children=nodeso))
    else:
        if update:
            Internal._rmNodesByName(a, outputName)
            nodeso = []
            for each in Dict: nodeso+=[Internal.createNode(each,'DataArray_t',value=Dict[each])]
            a[2].append(Internal.createNode(outputName, 'UserDefinedData_t', children=nodeso))
        else:
            for childname in Dict:
                childnode = Internal.getNodeFromName1(outputnode, childname)
                if childnode is None:
                    newchild = Internal.createNode(childname, 'DataArray_t', value=Dict[childname])
                    outputnode[2].append(newchild)
                else:# value exists-> updated
                    Internal.setValue(childnode, Dict[childname])
    return None

def addOutput(a, Dict, name='', update=False):
    """Add a node '.Solver#Output'+name to extract required value from a node.
    Usage: addOutput(a, Dict,name='',update=False)"""
    ap = Internal.copyRef(a)
    _addOutput(ap, Dict, name=name, update=update)
    return ap

# -----------------------------------------------------------------------------
def addOutputForces(node, name="", var=None, loc=4, writingmode=1,
                    period=None, pinf=None, fluxcoef=None,
                    torquecoef=None, xyztorque=None, frame=None,
                    governingEquations="NSTurbulent", xtorque=None,
                    ytorque=None, ztorque=None):
    """Add an output node for forces.
    Usage: addOutputForces(node,name,var,loc,writingmode,period,pinf,fluxcoef, torquecoef,wyztorque,frame,governingEquations,xtorque,ytorque,ztorque)"""
    nodep = Internal.copyRef(node)
    _addOutputForces(nodep, name, var, loc, writingmode,
                     period, pinf, fluxcoef,
                     torquecoef, xyztorque, frame,
                     governingEquations, xtorque, ytorque, ztorque)
    return nodep

def _addOutputForces(node, name="", var=None, loc=4, writingmode=1,
                     period=None, pinf=None, fluxcoef=None,
                     torquecoef=None, xyztorque=None, frame=None,
                     governingEquations="NSTurbulent", xtorque=None,
                     ytorque=None, ztorque=None):
    """Add an output node for forces.
    Usage: _addOutputForces(node,name,var,loc,writingmode,period,pinf,fluxcoef, torquecoef,wyztorque,frame,governingEquations,xtorque,ytorque,ztorque)"""
    Dict = {}
    Dict['var'] = "flux_rou flux_rov flux_row convflux_rou convflux_rov convflux_row diffflux_rou diffflux_rov diffflux_row torque_rou torque_rov torque_row convtorque_rou convtorque_rov convtorque_row difftorque_rou difftorque_rov difftorque_row"
    Dict["writingmode"] = writingmode
    Dict["loc"] = loc
    if var is not None: Dict["var"] = var
    elif governingEquations == "Euler": Dict["var"] = "convflux_rou convflux_rov convflux_row convtorque_rou convtorque_rov convtorque_row"
    if period is not None: Dict["period"] = period
    if pinf is not None: Dict["pinf"] = pinf
    if fluxcoef is not None: Dict["fluxcoeff"] = fluxcoef
    if torquecoef is not None: Dict["torquecoeff"] = torquecoef
    if xtorque is not None:
        Dict["xtorque"] = xtorque
        deprecation("'xtorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if ytorque is not None:
        Dict["ytorque"] = ytorque
        deprecation("'ytorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if ztorque is not None:
        Dict["ztorque"] = ztorque
        deprecation("'ztorque' is obsolete. Please use 'xyztorque' argument instead.", stacklevel=3)
    if xyztorque is not None:
        Dict["xtorque"] = xyztorque[0]
        Dict["ytorque"] = xyztorque[1]
        Dict["ztorque"] = xyztorque[2]
    if writingmode is not None: Dict["writingmode"] = writingmode
    if frame is not None:
        if frame not in ['relative','absolute']:
            raise AttributeError('Frame should be in %s')%(['relative','absolute'])
        Dict["writingframe"] = frame
    _addOutput(node, Dict, ":Forces"+name)
    return None

#==============================================================================
def addOutputFriction(node, name="", var=None, loc=4, writingmode=1,
                      period=None, fluxcoef=None,
                      torquecoef=None, writingframe=None):
    """Add an output node for frictions.
    Usage: addOutputFriction(node,name,var,loc,writingmode,period,fluxcoef,torquecoef,writingframe)"""
    nodep = Internal.copyRef(node)
    _addOutputFriction(node, name, var, loc, writingmode,
                       period, fluxcoef,
                       torquecoef, writingframe)
    return nodep

def _addOutputFriction(node, name="", var=None, loc=4, writingmode=1,
                       period=None, fluxcoef=None,
                       torquecoef=None, writingframe=None):
    """Add an output node for frictions.
    Usage: _addOutputFriction(node,name,var,loc,writingmode,period,fluxcoef,torquecoef,writingframe)"""
    Dict = {}
    Dict['var'] = "SkinFrictionX SkinFrictionY SkinFrictionZ SkinFrictionMagnitude WallCellSize"
    Dict["writingmode"] = writingmode
    Dict["loc"] = loc
    if var is not None: Dict["var"] = var
    if period is not None: Dict["period"] = period
    if fluxcoef is not None: Dict["fluxcoeff"] = fluxcoef
    if writingmode is not None: Dict["writingmode"] = writingmode
    if writingframe is not None: Dict["writingframe"] = writingframe
    _addOutput(node, Dict, ":Friction"+name)
    return None

#=============================================================================
def addGlobalConvergenceHistory(t, normValue=0):
    """Create a node for global convergence history storage for each base.
    Usage: addGlobalConvergenceHistory(t, normValue)"""
    tp = Internal.copyRef(t)
    _addGlobalConvergenceHistory(tp, normValue)
    return tp

def _addGlobalConvergenceHistory(t, normValue=0):
    for base in Internal.getBases(t):
        gcvhist=Internal.getNodeFromType1(base,'ConvergenceHistory_t')
        if gcvhist is None:
            child = Internal.createNode("NormDefinitions", "Descriptor_t", "ConvergenceHistory")
            base[2].append(Internal.createNode('GlobalConvergenceHistory','ConvergenceHistory_t',value=normValue, children=[child]))
        else:
            Internal.setValue(gcvhist,normValue)
    return None

#=============================================================================
def addReferenceState(t, conservative=None, temp=None, turbmod='spalart',
                      name='ReferenceState', comments=None):
    """Add a reference state for each base.
    Usage: addReferenceState(t, conservative, temp, turbmod, name, comments)"""
    tp = Internal.copyRef(t)
    _addReferenceState(tp, conservative, temp, turbmod, name, comments)
    return tp

def _addReferenceState(t, conservative=None, temp=None, turbmod='spalart',
                       name='ReferenceState', comments=None):
    """Add a reference state for each base.
    Usage: _addReferenceState(t, conservative, temp, turbmod, name, comments)"""
    strcomments = str(comments)

    if conservative is None:
        raise ValueError("addReferenceState: conservative arg must be defined as a list of strings.")
    nvars = len(conservative)
    if nvars < 5:
        raise ValueError("addReferenceState: conservative must be a list of length > 4.")

    varnames = ['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity']
    if nvars == 5: pass
    elif nvars == 6:
        if turbmod == 'spalart': varnames+=['TurbulentSANuTildeDensity']
        else: raise ValueError("Inconsistent model for 6-variable equations. Must be 'spalart'.")
    elif nvars == 7:
        if turbmod[0:6] =='komega':
            varnames.append('TurbulentEnergyKineticDensity')
            varnames.append('TurbulentDissipationRateDensity')
        elif turbmod[0:4]=='keps' or turbmod == 'chien' or turbmod == 'asm':
            varnames.append('TurbulentEnergyKineticDensity')
            varnames.append('TurbulentDissipationDensity')
        elif  turbmod == 'smith':
            varnames.append('TurbulentEnergyKineticDensity')
            varnames.append('TurbulentLengthScaleDensity')
        elif turbmod == 'kkl' or turbmod[0:5] == 'earsm':
            varnames.append('TurbulentEnergyKineticDensity')
            varnames.append('TurbulentEnergyKineticPLSDensity')
        else:
            raise ValueError("addReferenceState: %s is not a valid 2-equation turbulence model"%turbmod)
    elif nvars==12:
        if turbmod=='rsm':
            varnames.append("ReynoldsStressXX")
            varnames.append("ReynoldsStressXY")
            varnames.append("ReynoldsStressXZ")
            varnames.append("ReynoldsStressYY")
            varnames.append("ReynoldsStressYZ")
            varnames.append("ReynoldsStressZZ")
            varnames.append('ReynoldsStressDissipationScale')
        else:
            raise ValueError("addReferenceState: %s is not a valid turbulence model. Maybe 'rsm' ?")
    else:
        raise ValueError("addReferenceState: %d is not a valid number of variables. "%nvars)

    for b in Internal.getBases(t):

        childrenRS = []
        novar = 0
        for varname in varnames:
            childrenRS.append(Internal.createNode(varname,'DataArray_t',value=conservative[novar]))
            novar+=1

        if temp is not None:
            childrenRS.append(Internal.createNode('Temperature','DataArray_t',value=temp))

        childrenRS.append(Internal.createNode('ReferenceStateDescription','Descriptor_t',strcomments))
        NoeudReferenceState = Internal.createNode(name, 'ReferenceState_t',children=childrenRS)
        b[2].append(NoeudReferenceState)
    return None

#=============================================================================
def addFlowSolutionEoR(t, name='', variables=None, governingEquations=None,
                       writingFrame='relative', addBCExtract=False,
                       protocol="end"):
    """Add a node to extract the flow solution at the end of the run for each zone of the pyTree.
    Usage: addFlowSolutionEoR(t,name,variables,governingEquations,writingFrame,addBCExtract,protocol)"""
    tp = Internal.copyRef(t)
    _addFlowSolutionEoR(tp, name, variables, governingEquations,
                        writingFrame, addBCExtract, protocol)
    return tp

def _addFlowSolutionEoR(t, name='', variables=None, governingEquations=None,
                        writingFrame='relative', addBCExtract=False,
                        protocol="end"):
    _addFlowSolution(t, name=name+'#EndOfRun', loc='CellCenter', governingEquations=governingEquations, writingFrame=writingFrame, variables=variables, addBCExtract=addBCExtract, protocol=protocol)
    return None

#===============================================================================
def addFlowSolution(t, name='', loc='CellCenter', variables=None,
                    governingEquations=None, writingMode=None,
                    writingFrame='relative', period=None, output=None,
                    addBCExtract=False, protocol="end"):
    """Add a node to extract the flow solution.
    Usage: addFlowSolution(t, name, loc, variables, governingEquations, writingFrame, period, addBCExtract, protocol)"""
    tp = Internal.copyRef(t)
    _addFlowSolution(tp, name, loc, variables,
                     governingEquations, writingMode,
                     writingFrame, period, output,
                     addBCExtract, protocol)
    return tp

def _addFlowSolution(t, name='', loc='CellCenter', variables=None,
                     governingEquations=None, writingMode=None,
                     writingFrame='relative', period=None, output=None,
                     addBCExtract=False, protocol="end"):
    """Add a node to extract the flow solution."""
    if governingEquations is not None: gE0 = governingEquations
    else:
        gE0 = None
        GE = Internal.getNodeFromType2(t,'GoverningEquations_t')
        if GE is not None: gE0 = Internal.getValue(GE)
    if variables is not None:
        for c, v in enumerate(variables): variables[c] = Internal.getCGNSName(v)
    if output is None: outputDict={}
    else: outputDict=output
    if loc == 'cellfict': outputDict["loc"]=2
    elif loc == 'CellCenter' or loc == "Vertex": pass
    else: raise AttributeError("'loc' attribute should be 'CellCenter','Vertex' or 'cellfict' => loc='%s'" %(loc))
    if writingMode is not None: outputDict['writingmode'] = writingMode
    if period is not None: outputDict['period'] = period
    if writingFrame is not None: outputDict['writingframe'] = writingFrame

    for base in Internal.getBases(t):
        if gE0 is None:
            GE = Internal.getNodeFromType2(base, 'GoverningEquations_t')
            if GE is not None:
                gE0 = Internal.getValue(GE)
        for z in Internal.getZones(base):
            if gE0 is None:
                GE = Internal.getNodeFromType2(z, 'GoverningEquations_t')
                if GE is None:
                    if variables is None:
                        raise ValueError("addFlowSolution: no GoverningEquations in tree. You must specify it in your function parameter.")
                    else: gE0 = None
                else: gE0 = Internal.getValue(GE)

            childrenNodes = []
            childrenNodes.append(Internal.createNode('GridLocation','GridLocation_t',value=loc))
            if variables is None:
                if gE0 is None:
                    if variables is None: variables = []

                else:
                    variables = __CONSERVATIVE__+__COMMONS__
                    if gE0 != 'Euler': variables+=__COMMONSNS__
                    if gE0 == 'NSTurbulent': variables += __TURBULENT__+__WALLDISTANCE__
                    #print('addFlowSolution: extracted variables are: ',variables)
            else:
                if isinstance(variables,str): variables = variables.split()
                if "xyz" in variables:
                    variables+= __XYZ__; novar = variables.index("xyz"); variables.pop(novar)
                if "Conservative" in variables:
                    variables += __CONSERVATIVE__; novar = variables.index("Conservative"); variables.pop(novar)
                if 'Turbulent' in variables:
                    variables+=__TURBULENT__; novar = variables.index("Turbulent"); variables.pop(novar)
                if "Commons" in variables:
                    variables+=__COMMONS__; novar = variables.index("Commons"); variables.pop(novar)
                if 'CommonsNS' in variables:
                    variables+=__COMMONSNS__; novar = variables.index("CommonsNS"); variables.pop(novar)
                if 'WallDistance' in variables:
                    variables+=__WALLDISTANCE__; novar = variables.index("WallDistance"); variables.pop(novar)
                variables= list(set(variables))

            for varname in variables:
                childrenNodes.append(Internal.createNode(varname,'DataArray_t',value=None))

            if addBCExtract:
                subChildren1 = []
                for bc in Internal.getNodesFromType2(z, "BC_t"):
                    prval = Internal.getValue(Internal.getNodeFromName1(bc, "PointRange"))
                    subChildren2 = [Internal.createNode("PointRange","DataArray_t", value=prval)]
                    subChildren2+= [Internal.createNode("Protocol","Descriptor_t",value=protocol)]
                    strval = '\n'.join(['SurfaceSolution/NeumannData/'+variable for variable in variables])
                    subChildren2+= [Internal.createNode("Path","DataArray_t",value=strval)]
                    subChildren1+= [Internal.createNode(Internal.getName(bc), "UserDefinedData_t",children=subChildren2)]

                child = Internal.createNode(".Solver#SolutionSubSetInBC", "UserDefinedData_t",children=subChildren1)
                childrenNodes += [child]

            FS = Internal.createNode("FlowSolution"+name,"FlowSolution_t",children=childrenNodes)
            if outputDict != {}: _addOutput(FS, outputDict)
            # add the flow solution to zone node
            z[2].append(FS) # SP : a mettre en unique or not ?

    return None

#==============================================================================
def _buildMaskFiles(t, keepOversetHoles=True, fileDir='.', prefixBase=False):
    dataFormat='%14.7e'
    c = 0
    for base in Internal.getBases(t):
        basename = Internal.getName(base)
        for z in Internal.getZones(base):
            ho = Internal.getNodeFromType2(z, 'OversetHoles_t')
            if ho is not None:
                pl = Internal.getNodeFromName1(ho,'PointList')
                h = Internal.getValue(pl)
                dim = Internal.getZoneDim(z)
                h = Internal.convertIJKArray21DArray(h, dim[1]-1,dim[2]-1,dim[3]-1)
                hf = numpy.empty(h.shape, dtype=numpy.float64)
                hf[:] = h[:]
                array = ['cell_index', hf, hf.size, 1, 1]
                fileName = fileDir+'/hole_'
                if prefixBase: fileName += basename+'_'
                Converter.convertArrays2File([array], fileName+z[0]+'.v3d',
                                             "bin_v3d", dataFormat=dataFormat)

        if not keepOversetHoles: Internal._rmNodesByName(t, 'OversetHoles')
    return None

#==============================================================================
def buildMaskFiles(t, keepOversetHoles=True, fileDir='.', prefixBase=False):
    """Write the mask files for elsA solver."""
    tp = Internal.copyRef(t)
    _buildMaskFiles(tp, keepOversetHoles, fileDir, prefixBase)
    return tp

#==============================================================================
# set BCOverlap - version 1: a la Thomas
#==============================================================================
def overlapGC2BC(t):
    """ Convert the Overlap boundary conditions from GC (Grid Connectivity) to BC (Boundary Condition) for elsA solver.
    Usage: overlapGC2BC(t)"""
    tp = Internal.copyRef(t)
    _overlapGC2BC(tp)
    return tp

def _overlapGC2BC(t):
    donorDict={}
    for nob in range(len(t[2])):
        if Internal.getType(t[2][nob])=='CGNSBase_t':
            base = t[2][nob]
            baseName = base[0]
            isDD = False; isClassical=False
            oversetBase=False
            for noz in range(len(base[2])):
                zone = base[2][noz]
                if Internal.getType(zone)=='Zone_t':
                    zoneName = Internal.getName(zone)
                    overlapgcs = []
                    overlapbcs = []
                    zgc = Internal.getNodesFromType1(zone,"ZoneGridConnectivity_t")
                    zbc = Internal.getNodeFromType1(zone,'ZoneBC_t')
                    if zgc is not None:
                        gcs = Internal.getNodesFromType1(zgc,"GridConnectivity_t")
                        for gc in gcs:
                            gctype = Internal.getNodeFromType1(gc,"GridConnectivityType_t")
                            gctype = Internal.getValue(gctype)
                            if gctype=='Overset': overlapgcs.append(gc); oversetBase=True

                        for ov in overlapgcs:
                            PR = Internal.getNodeFromName1(ov,"PointRange")
                            PRBC = Internal.createNode(Internal.getName(PR),Internal.getType(PR),value=Internal.getValue(PR))
                            bcval = 'UserDefined'
                            foundDD = False
                            ovsons = [PRBC]
                            for childov in Internal.getChildren(ov):
                                if Internal.getType(childov)=='UserDefinedData_t':
                                    if Internal.getNodeFromName1(childov,"doubly_defined") is not None:
                                        isDD = True
                                        foundDD = True
                                        bcval = 'UserDefinedDD' # pour ne pas les grouper  avec les classiques
                                        dnrList=""
                                        # donors: list of zones or family of zones
                                        for dnr in Internal.getValue(ov).split(','):
                                            if dnrList=="": dnrList+=dnr
                                            else: dnrList+=" "+dnr

                                        # family dd
                                        famOvlpName = __FAMOVERLAPDDBC__+ov[0]
                                        DDnode=Internal.createNode('doubly_defined','DataArray_t',value='active')
                                        NLnode=Internal.createNode('NeighbourList','DataArray_t',value=dnrList)
                                        FBC = Internal.createNode("FamilyBC","FamilyBC_t","UserDefined")
                                        SOV = Internal.createNode(".Solver#Overlap","UserDefinedData_t",children=[NLnode,DDnode])
                                        famovsons = [FBC,SOV]
                                        Internal.createNode(famOvlpName,'Family_t',children=famovsons,parent=t[2][nob])
                                        ovsons.append(Internal.createNode("FamilyName","FamilyName_t",value=famOvlpName))

                            ovbc = Internal.createNode(ov[0],"BC_t",value=bcval,children=ovsons)
                            overlapbcs.append(ovbc)
                            if not foundDD: isClassical = True

                        if zbc is None:
                            newZoneBC = Internal.createNode("ZoneBC","ZoneBC_t",value=None,children=overlapbcs)
                            Internal._addChild(t[2][nob][2][noz],newZoneBC)
                        else:
                            for overlapbc in overlapbcs:
                                Internal._addChild(zbc,overlapbc)

            if oversetBase:
                if isClassical:
                    famOvlpName = __FAMOVERLAPBC__+baseName
                    Internal._groupBCByBCType(base, btype='UserDefined', name=famOvlpName)
                    FamOvlp = Internal.getNodeFromName1(base, famOvlpName)
                    if FamOvlp:
                        NLnode = Internal.createNode('NeighbourList', 'DataArray_t')
                        Internal._createChild(FamOvlp, '.Solver#Overlap', 'UserDefinedData_t',
                                              value=None, children=[NLnode])

            # renommage en FamilySpecified des doubly defined UserDefinedDD
            for bc in Internal.getNodesFromType3(t[2][nob],'BC_t'):
                if Internal.getValue(bc)=='UserDefinedDD':
                    Internal.setValue(bc,'FamilySpecified')

    return None

#==============================================================================
def _rmGCOverlap(t):
    for z in Internal.getZones(t):
        allzgcs = Internal.getNodesFromType1(z,'ZoneGridConnectivity_t')
        for zgc in allzgcs:
            gc = Internal.getNodesFromType1(zgc,'GridConnectivity_t')
            for node in gc:
                gct = Internal.getNodeFromType1(node,'GridConnectivityType_t')
                if Internal.getValue(gct)=='Overset': Internal._rmNode(zgc,node)
    return t

#==============================================================================
def rmGCOverlap(t):
    """ Remove the Overlap boundary conditions described as Grid Connectivities (GC).
    Usage: rmGCOverlap(t)"""
    tp = Internal.copyRef(t)
    _rmGCOverlap(tp)
    return tp

#==============================================================================
# fill NeigbourList by families of zones according to the dict of intersecting
# zones for all zones
# if sameBase=1 allows for donors in same base
# classical overlaps are gathered in a same family F_OV_*
# donor zones are tagged with ChimGroup_*
# doubly defined overlaps: one family per doubly defined overlap BC
# in that case, donors zones must be specified in the NeighbourList already
# as a list of zones and/or families of zones
#==============================================================================
def fillNeighbourList(t, sameBase=0):
    """Fill neighbour list by families of zones according to intersection."""
    tp = Internal.copyRef(t)
    _fillNeighbourList(tp, sameBase=sameBase)
    return tp

def _fillNeighbourList(t, sameBase=0):
    import Generator.PyTree as G
    dictOfNobOfZone={}
    dictOfNozOfZone={}
    for nobR in range(len(t[2])):
        baseR = t[2][nobR]
        if Internal.getType(baseR)=='CGNSBase_t':
            for nozR in range(len(baseR[2])):
                zr = baseR[2][nozR]
                if Internal.getType(zr)=='Zone_t':
                    dictOfNobOfZone[zr[0]]=nobR
                    dictOfNozOfZone[zr[0]]=nozR

    nbOfDD = 1
    nbOfG = 1
    tBBD = G.BB(t)
    tBBR = C.newPyTree()
    for base in Internal.getBases(tBBD):
        [xmin,ymin,zmin,xmax,ymax,zmax] = G.bbox(base)
        a = G.cart((xmin,ymin,zmin),(xmax-xmin,ymax-ymin,zmax-zmin),(2,2,2))
        C._addBase2PyTree(tBBR,base[0])
        a[0] = base[0]
        tBBR[2][-1][2].append(a)

    if sameBase == 1: intersectionDict = X.getIntersectingDomains(tBBR,t2=tBBD,method='AABB',taabb=tBBR,taabb2=tBBD)
    else:
        intersectionDict={}
        for baseR in Internal.getBases(tBBR):
            tBBDL = C.newPyTree()
            for baseD in Internal.getBases(tBBD):
                if baseD[0] != baseR[0]: tBBDL[2].append(baseD)
            intersectionDictR = X.getIntersectingDomains(baseR,t2=tBBDL,method='AABB',taabb=tBBR,taabb2=tBBDL)
            intersectionDict[baseR[0]]=intersectionDictR[baseR[0]]
    for nobR in range(len(t[2])):
        baseR = t[2][nobR]
        basenameR = baseR[0]
        NList=""
        chimgroupname = __CHIMGROUPNAME__+str(nbOfG)
        if Internal.getType(baseR)=='CGNSBase_t':
            listOfDnrBasesNob=[]
            dnrZoneNames0 = intersectionDict[basenameR]
            for zdname in dnrZoneNames0:
                nozd = dictOfNozOfZone[zdname]
                nobd = dictOfNobOfZone[zdname]
                if sameBase == 1 or (sameBase==0 and nobd != nobR):
                    Internal.createNode("FamilyName", "FamilyName_t", value=chimgroupname, parent=t[2][nobd][2][nozd])

                # create ChimGroup for donor base
                if nobd not in listOfDnrBasesNob:
                    listOfDnrBasesNob.append(nobd)
                    chimgroupsons = [Internal.createNode("FamilyBC","FamilyBC_t",value='UserDefined')]
                    chmg=Internal.createNode(chimgroupname,"Family_t",children=chimgroupsons)
                    t[2][nobd][2].append(chmg)
                    basenameD = t[2][nobd][0]
                    if NList=="": NList=basenameD+'/'+chimgroupname
                    else: NList+=" "+basenameD+'/'+chimgroupname
                    #print 'NeighbourList = ', NList, ' for chimgroup ', chimgroupname
                    nbOfG+=1

            if NList != "":
                NLnode = Internal.getNodeFromName1(baseR,__FAMOVERLAPBC__+basenameR)
                if NLnode is not None:
                    NLnode = Internal.getNodeFromName2(NLnode,'NeighbourList')
                    Internal.setValue(NLnode,NList)

            # doubly defined : passer les NL definies sous forme de noms de zones en familles de zones
            allFAMDD = Internal.getNodesFromName(baseR,__FAMOVERLAPDDBC__+'*')
            for famdd in allFAMDD:
                famddname = famdd[0]
                chimgroupnamedd=__CHIMGROUPNAMEDD__+str(nbOfDD)
                NList = Internal.getNodesFromName2(famdd,"NeighbourList")
                if NList != []:
                    donorList= Internal.getValue(NList[0])
                    dnrListNL = ""
                    listOfDnrBasesNob = []
                    for dnrname in donorList.split(' '):
                        if dnrname not in dictOfNozOfZone:
                            # c est une famille de zones, la prefixer par sa base
                            dnrname2 = dnrname.split('/')
                            if len(dnrname2) == 1:
                                for baseLoc in Internal.getBases(t):
                                    myFamLoc = Internal.getNodesFromType1(baseLoc,"Family_t")
                                    if myFamLoc != []:
                                        myFamLoc = Internal.getNodesFromName1(myFamLoc,dnrname2[0])
                                        if myFamLoc != []:
                                            myFamLoc = myFamLoc[0]
                                            if Internal.getNodeFromType1(myFamLoc,"FamilyBC_t") is None:
                                                Internal.createNode("FamilyBC","FamilyBC_t",parent=myFamLoc)
                                            baseName = baseLoc[0]
                                            dnrname = baseName+'/'+dnrname
                                            if dnrListNL =="": dnrListNL=dnrname
                                            else: dnrListNL += " "+dnrname
                            else:
                                if dnrListNL =="": dnrListNL=dnrname
                                else: dnrListNL += " "+dnrname
                        else:
                            # create a chim group
                            nozd = dictOfNozOfZone[dnrname]
                            nobd = dictOfNobOfZone[dnrname]
                            Internal.createNode("FamilyName", "FamilyName_t", value=chimgroupnamedd, parent=t[2][nobd][2][nozd])
                            # create ChimGroup for donor base
                            if nobd not in listOfDnrBasesNob:
                                listOfDnrBasesNob.append(nobd)
                                chimgroupsons = [Internal.createNode("FamilyBC","FamilyBC_t",value='UserDefined')]
                                chmg=Internal.createNode(chimgroupnamedd,"Family_t",children=chimgroupsons)
                                t[2][nobd][2].append(chmg)
                                basenameD = t[2][nobd][0]
                                if dnrListNL == "": dnrListNL=basenameD+'/'+chimgroupnamedd
                                else: dnrListNL+=" "+basenameD+'/'+chimgroupnamedd

                    #print 'NeighbourList = ', NList[0], " : ", dnrListNL, ' for dd chimgroup ', chimgroupnamedd
                    Internal.setValue(NList[0],dnrListNL)
                    nbOfDD+=1

    _cleanDonorFamilyNodes(t, name=__CHIMGROUPNAME__+'*')
    return None

def _cleanDonorFamilyNodes(t, name=__CHIMGROUPNAME__+'*'):
    listOfFamNames = []
    for base in Internal.getBases(t):
        famNodes = Internal.getNodesFromType1(base,'Family_t')
        famNodes = Internal.getNodesFromName(famNodes,name)
        for fam in famNodes:
            listOfFamNames.append(Internal.getName(fam))
    listOfFamNames = list(set(listOfFamNames))

    listOfReferencedFamNames=[]
    for base in Internal.getBases(t):
        for SOV in Internal.getNodesFromName2(base,".Solver#Overlap"):
            for NL in Internal.getNodesFromName1(SOV,"NeighbourList"):
                NLVal = Internal.getValue(NL)
                NList = NLVal.split(' ')
                for NL0 in NList:
                    listOfReferencedFamNames.append(NL0.split('/')[1])

    listOfReferencedFamNames = list(set(listOfReferencedFamNames))
    for famName in listOfFamNames:
        if famName not in listOfReferencedFamNames:
            Internal._rmNodesByNameAndType(t, famName, 'Family_t')
            Internal._rmNodesByNameAndType(t, famName, 'FamilyName_t')

    return None

#==============================================================================
def _addTurbulentDistanceIndex(t):
    """Add the TurbulentDistance index node for elsA solver (in place version)."""
    for z in Internal.getZones(t):
        sol = Internal.getNodeFromName1(z,Internal.__FlowSolutionCenters__)
        if sol is not None:
            distindex = Internal.getNodeFromName1(sol,'TurbulentDistanceIndex')
            if distindex is None:
                dist = Internal.getNodeFromName1(sol,'TurbulentDistance')
                if dist is not None:
                    #indexa = numpy.full(dist[1].shape,-1, dtype=Internal.__E_NPY_INT, order='F')
                    #indexa = numpy.empty(dist[1].shape, dtype=Internal.__E_NPY_INT, order='F')
                    indexa = numpy.empty(dist[1].shape, dtype=numpy.float64, order='F')
                    indexa[:] = -1
                    distindex = Internal.createNode("TurbulentDistanceIndex",'DataArray_t',value=indexa,parent=sol)
    return None

#==============================================================================
def addTurbulentDistanceIndex(t):
    """Add the TurbulentDistance index node for elsA solver.
    Usage: addTurbulentDistanceIndex(t)"""
    tp = Internal.copyRef(t)
    _addTurbulentDistanceIndex(tp)
    return tp

#==============================================================================
def adaptNearMatch(t):
    """Convert Nearmatch Grid Connectivity (GC) to be compliant with elsA solver.
    Usage: adaptNearmatch(t)"""
    tp  = Internal.copyRef(t)
    _adaptNearMatch(tp)
    return tp

def _adaptNearMatch(t):
    for z in Internal.getZones(t):
        zgc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
        if zgc is not None:
            allgcs = Internal.getNodesFromType1(zgc, 'GridConnectivity_t')
            for gc in allgcs:
                udd = Internal.getNodeFromType1(gc, 'UserDefinedData_t')
                if udd is not None:
                    nm = Internal.getNodeFromName1(udd, 'NMRatio')
                    if nm is not None:
                        Internal.setType(gc,'GridConnectivity1to1_t')
                        transfo = Internal.getNodeFromName1(udd,'Transform')
                        PRD = Internal.getNodeFromName1(udd,'PointRangeDonor')
                        if transfo is None or PRD is None:
                            raise ValueError("adaptNearMatch: near match gc %s is not valid."%(gc[0]))

                        PRD2 = Internal.createNode("PointRangeDonor","IndexRange_t",value=Internal.getValue(PRD),parent=gc)
                        transfo2 = Internal.createNode('Transform','\"int[IndexDimension]\"',value=Internal.getValue(transfo),parent=gc)
                        Internal._rmNodesByName(gc,udd[0])
                        Internal._rmNodesByName(gc,"PointListDonor")
                        Internal._rmNodesByName(gc,'GridConnectivityType')
                        # fine or coarse
                        fine = 0; iratio=1; jratio=1; kratio=1
                        nmr = Internal.getValue(nm)
                        if nmr.shape[0] != 3: raise ValueError("adaptNearMatch: NearMatchRatio must be of shape (3,)")
                        nmratioFact=nmr[0]*nmr[1]*nmr[2]
                        matchside = ''
                        if nmratioFact<1: # fine
                            matchside='fine'
                            iratio = int(math.ceil(1./nmr[0]))
                            jratio = int(math.ceil(1./nmr[1]))
                            kratio = int(math.ceil(1./nmr[2]))
                        else:
                            matchside='coarse'
                            iratio = int(nmr[0])
                            jratio = int(nmr[1])
                            kratio = int(nmr[2])

                        #Create .Solver#Property node
                        spsons = []
                        spsons += [Internal.createNode("jtype","DataArray_t",value='nearmatch')]
                        spsons += [Internal.createNode("type","DataArray_t",value='join')]
                        spsons += [Internal.createNode("matchside","DataArray_t",value=matchside)]
                        spsons += [Internal.createNode("i_ratio","DataArray_t",value=iratio)]
                        spsons += [Internal.createNode("j_ratio","DataArray_t",value=jratio)]
                        spsons += [Internal.createNode("k_ratio","DataArray_t",value=kratio)]
                        solverProperty = Internal.createNode(".Solver#Property","UserDefinedData_t",children=spsons,parent=gc)

    return None

#=========================================================================================
def createElsaHybrid(t, method=0, axe2D=0, methodPE=0):
    """Create elsAHybrid node necessary for NGON zones."""
    tp = Internal.copyRef(t)
    _createElsaHybrid(tp, method, axe2D, methodPE)
    return tp

def _createElsaHybrid(t, method=0, axe2D=0, methodPE=0):
    from . import converter
    zones = Internal.getZones(t)
    for z in zones:
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        if found:
            node = GEl[NGON]
            CE = Internal.getNodeFromName1(node, 'ElementConnectivity')
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                Internal._adaptNFace2PE(z, remove=False, methodPE=methodPE)
                PE = Internal.getNodeFromName1(node, 'ParentElements')
            er = Internal.getNodeFromName1(node, 'ElementRange')
            nfaces = er[1][1]-er[1][0]+1
            child = Internal.createUniqueChild(z, ':elsA#Hybrid', 'UserDefinedData_t')
            ESO = Internal.getNodeFromName1(node, 'ElementStartOffset')
            if ESO is not None: ESO = ESO[1]

            # to be removed (only used by elsA for nothing)
            sct = numpy.arange((nfaces), dtype=Internal.E_NpyInt)
            Internal.newDataArray('SortedCrossTable', value=sct, parent=child)
            inct = numpy.empty((nfaces), dtype=Internal.E_NpyInt)
            if ESO is None: Internal.newDataArray('IndexNGONCrossTable', value=inct, parent=child)
            # OK
            ict = -1*numpy.ones((nfaces), dtype=Internal.E_NpyInt)
            bcct = -1*numpy.ones((nfaces), dtype=Internal.E_NpyInt)
            Internal.newDataArray('InversedCrossTable', value=ict, parent=child)
            Internal.newDataArray('BCCrossTable', value=bcct, parent=child)
            if axe2D > 0:
                x = Internal.getNodeFromName2(z, 'CoordinateX')[1]
                y = Internal.getNodeFromName2(z, 'CoordinateY')[1]
                z = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
            else: x = None; y = None; z = None
            (iTRI, iQUADS, eTRI, eQUADS) = converter.createElsaHybrid(
                CE[1], PE[1], ict, bcct, inct, method,
                axe2D, x, y, z, ESO)
            if method == 0:
                Internal.newDataArray('InternalTris', iTRI, parent=child)
                Internal.newDataArray('InternalQuads', iQUADS, parent=child)
                Internal.newDataArray('ExternalTris', eTRI, parent=child)
                Internal.newDataArray('ExternalQuads', eQUADS, parent=child)
            else:
                Internal.newDataArray('InternalElts', iTRI, parent=child)
                Internal.newDataArray('ExternalElts', eTRI, parent=child)
        else:
            pass
            #print('Warning: createElsaHybrid: no NGON node found for zone %s. No :elsAHybrid node created.'%z[0])
    return None

#==============================================================================
# Prefix ID_* zone subregion donors
#==============================================================================
def prefixDnrInSubRegions(t):
    """Prefix zone subregion ID_* donor names with base name."""
    tp = Internal.copyRef(t)
    _prefixDnrInSubRegions(tp)
    return tp

# in place version
def _prefixDnrInSubRegions(t):
    if Internal.getNodeFromType3(t,"ZoneSubRegion_t") is None: return None

    baseDict={}
    for base in Internal.getBases(t):
        basename = Internal.getName(base)
        for z in Internal.getZones(base):
            zname = Internal.getName(z)
            baseDict[zname] = basename

    for z in Internal.getZones(t):
        subRegions= Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        for s in subRegions:
            sname = Internal.getName(s)
            if sname.split('_')[0] == 'ID':
                dnrname = Internal.getValue(s)
                s2 = dnrname.split('/')
                if len(s2)==1:
                    baseName = baseDict[dnrname]
                    newname = baseName+'/'+dnrname
                    Internal.setValue(s,newname)

    return None

# remove OrphanPointList nodes and  subregions with empty PointList nodes (created for orphan points)
def _cleanIDSubregions(t):
    for z in Internal.getZones(t):
        removedNames=[]
        for zsr in Internal.getNodesFromType1(z,'ZoneSubRegion_t'):
            PL = Internal.getNodeFromName1(zsr,'PointList')
            if PL[1].shape[0]==0: removedNames.append(zsr[0])
            Internal._rmNodesFromName(zsr,"OrphanPointList")
        for srname in removedNames:
            Internal._rmNodesFromName1(z,srname)
    return None
#==============================================================================
# Conversion d'un arbre CGNS a la Cassiopee en un arbre de profil elsAxdt
#==============================================================================
def convert2elsAxdt(t, sameBase=0, fileDir='.'):
    """Perform all necessary transformations to obtain a computable tree for elsA."""
    tp = Internal.copyRef(t)
    _convert2elsAxdt(tp, fileDir)
    return tp

def _convert2elsAxdt(t, sameBase=0, fileDir='.'):
    print('1. addTurbulentDistance index')
    _addTurbulentDistanceIndex(t)
    print('2. buildMaskFiles')
    _buildMaskFiles(t, fileDir=fileDir)
    print('3. adaptNearMatch')
    _adaptNearMatch(t)
    print('4. adaptPeriodicMatch')
    _adaptPeriodicMatch(t, clean=True)
    print('5. overlapGC2BC')
    _overlapGC2BC(t)
    print('6. rmGCOverlap')
    _rmGCOverlap(t)
    print('7. fillNeighbourList')
    _fillNeighbourList(t)
    print('8. prefixDnrInSubRegions')
    _prefixDnrInSubRegions(t)
    print('9. clean subregions ID')
    _cleanIDSubregions(t)
    return None

#===============================================================================================================================
#===============================================================================================================================
# PARTIE UTILISEE DANS L ANCIEN CONVERT2ELSAXDT ET NON DOCUMENTEE POUR L INSTANT
#===============================================================================================================================
#===============================================================================================================================
#==============================================================================
def adaptNearmatch__(t):
    print('Warning: elsAProfile: adaptNearmatch__ is obsolete. New function name is adaptNearMatch.')
    return adaptNearMatch(t)
def adaptNearmatch(t):
    print('Warning: elsAProfile: adaptNearmatch is obsolete. New function name is adaptNearMatch.')
    return adaptNearMatch(t)

#==============================================================================
# Obsolete
#==============================================================================
#==============================================================================
# OBSOLETE - SERA SUPPRIME DANS UNE PROCHAINE RELEASE
def addNeighbours__(t, sameBase=0):
    """ Fill the NeighbourList nodes with bounding-box domains intersection.
    Usage: addNeighbours(t,sameBase)"""
    tp = Internal.copyRef(t)
    _addNeighbours__(tp,sameBase=sameBase)
    return tp

#==============================================================================
def _addNeighbours__(t, sameBase=0):
    """ Fill the NeighbourList nodes with bounding-box domains intersection.
    """
    bases = Internal.getBases(t)
    for i in range(len(bases)):
        fams = []
        doms = X.getCEBBIntersectingDomains(bases[i] , bases, sameBase)
        if doms != [[]]:
            for j in range(len(doms)):
                famsByZone=[]
                if doms[j] !=[]:
                    for donorZoneName in doms[j]:
                        z = Internal.getNodesFromName2(t,donorZoneName)
                        if z != []:
                            (donorBase, num) = Internal.getParentOfNode(t,z[0])
                            baseName=donorBase[0]
                            donorFamilyName='F_'+ donorZoneName
                            donorF = Internal.getNodesFromName2(t,donorFamilyName)[0]
                            famsByZone.append(baseName+'/'+donorF[0])
                fams.append(famsByZone)
        zones = Internal.getNodesFromType1(bases[i],'Zone_t')
        for zi in range(len(zones)):
            familyName='F_'+ zones[zi][0]
            F = Internal.getNodesFromName1(bases[i],familyName)[0]
            lOvlp = Internal.getNodesFromName(F,'.Solver#Overlap')
            if lOvlp != []:
                Ovlp = lOvlp[0]
                N = Internal.getNodesFromName(F,'NeighbourList')[0]
                (parentBase, numFamily) = Internal.getParentOfNode(bases[i],F)
                (parentFamily, numOvlp) = Internal.getParentOfNode(F,Ovlp)
                (parentOvlp, numN) = Internal.getParentOfNode(Ovlp,N)
                listFamily = ''
                for l in range(len(fams[zi])):
                    listFamily = listFamily+fams[zi][l]
                    if l != len(fams[zi])-1:
                        listFamily = listFamily + ' '
                v = numpy.array([c for c in listFamily], 'c')

                N[1] = v
                Ovlp[2][numN] = N
                F[2][numOvlp] = Ovlp
                bases[i][2][numFamily] = F
            # famille doubly-defined
            gc = Internal.getNodesFromType2(zones[zi],'GridConnectivity_t')
            for k in range(len(gc)):
                familyNameDD='FDD_'+ zones[zi][0]+'_'+gc[k][0]
                FDD = (Internal.getNodesFromName1(bases[i],familyNameDD))
                if FDD != []:
                    listVal="";domsDD=""
                    if isinstance(gc[k][1], numpy.ndarray):
                        val = gc[k][1].tobytes().decode()
                        listVal = val.split(",")
                    for vi in range(len(listVal)):
                        nvi = (Internal.getNodesFromName(t,listVal[vi]))[0]
                        (pvi, nvi) = Internal.getParentOfNode(t,nvi)
                        namevi = pvi[0]+'/F_'+listVal[vi]
                        if vi != len(listVal)-1:
                            namevi = namevi + ' '
                        domsDD = domsDD + namevi
                    OvlpDD  = Internal.getNodesFromName(FDD[0],'.Solver#Overlap')[0]
                    NDD  = Internal.getNodesFromName(FDD[0],'NeighbourList')[0]
                    (parentBaseDD, numFamilyDD) = Internal.getParentOfNode(bases[i],FDD[0])
                    (parentFamilyDD, numOvlpDD) = Internal.getParentOfNode(FDD[0],OvlpDD)
                    (parentOvlpDD, numNDD) = Internal.getParentOfNode(OvlpDD,NDD)
                    Internal._setValue(NDD, domsDD)
                    OvlpDD[2][numNDD] = NDD
                    FDD[0][2][numOvlpDD] = OvlpDD
                    bases[i][2][numFamilyDD] = FDD[0]

    cgnsv = Internal.getNodeFromType1(t, 'CGNSLibraryVersion_t')
    (parent, pos) = Internal.getParentOfNode(t,cgnsv)
    bases = [t[2][pos]]+bases
    t[2] = bases
    return t
#==============================================================================
def addFamilyBCNode__(t):
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    for i in range(len(bases)):
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for f in families:
            (parentBase, numFamily) = Internal.getParentOfNode(bases[i], f)
            fbc = Internal.getNodesFromType1(f, 'FamilyBC_t')
            if fbc == []:
                Internal._createChild(f, 'FamilyBC', 'FamilyBC_t', value='UserDefined')
            else:
                fbc[0][0] = 'FamilyBC'
            bases[i][2][numFamily] = f

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp, cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases
    return tp

def buildBCOverlap(t):
    """OBSOLETE: use ovelrapGC2BC"""
    print('WARNING: elsAProfile.buildBCOverlap is obsolete. Will be removed in next release. Please replace by overlapGC2BC in your script.')
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)

    c = 0 # compteur pour le nommage des conditions BCOverlap
    for i in range(len(bases)):
        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in range(len(zones)):
            (parentBase, numZone) = Internal.getParentOfNode(tp,zones[j])
            # Creation d'une famille par zone"
            familyName='F_'+zones[j][0]
            hasOverset = 0
            zonegct = Internal.getNodesFromType(zones[j],'GridConnectivityType_t')
            for zgct in zonegct:
                valz = Internal.getValue(zgct)
                if valz == 'Overset': hasOverset = 1
            if hasOverset == 0: # pas de connectivite Overset dans la zone => famille "vide"
                # Creation de la famille
                bases[i] = C.addFamily2Base(bases[i], familyName, 'UserDefined')
                # Ajout pour la zone d'un noeud fils FamilyName, pour faire reference a la famille creee
                famNameArray = numpy.array([c for c in familyName], 'c')
                familyNameList = Internal.getNodesFromType1(zones[j], 'FamilyName_t')
                if familyNameList == []:
                    zones[j][2].append(['FamilyName', famNameArray, [], 'FamilyName_t'])
                else:
                    addFamilyNameList = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')
                    if addFamilyNameList == []:
                        zones[j][2].append(['AddFamilyName', famNameArray, [], 'AdditionalFamilyName_t'])
            else:
                bases[i] = C.addFamily2Base(bases[i], familyName, 'BCOverlap')
                # Creation d'un noeud NeighbourList dans le noeud .Solver#Overlap de la famille creee
                F = Internal.getNodesFromName1(bases[i], familyName)[0]
                (parentBase, numFamily) = Internal.getParentOfNode(bases[i],F)
                Ovlp = Internal.getNodesFromName(F,'.Solver#Overlap')[0]
                (parentFamily, numOvlp) = Internal.getParentOfNode(F,Ovlp)
                F[2][numOvlp][2].append(['NeighbourList', None, [], 'DataArray_t'])
                bases[i][2][numFamily] = F
                # Creation des noeuds ZoneBC_t de type BCOverlap
                gc = Internal.getNodesFromType2(zones[j],'GridConnectivity_t')
                for k in range(len(gc)):
                    # search in parent GridConnectivity if a PointRange is present
                    # if no, connectivity is linked to blanking and not to a BCOverlap
                    prange = Internal.getNodesFromName1(gc[k],'PointRange')
                    if prange != []: # corresponds to a BCOverlap
                        gct = Internal.getNodesFromType3(gc[k],'GridConnectivityType_t')
                        for o in gct:
                            val = Internal.getValue(o)

                            if val == 'Overset':
                                # Recuperation des donnees necessaires a la creation d'un noeud ZoneBC_t
                                #   range
                                lPointRange = Internal.getNodesFromName1(gc[k], 'PointRange')
                                r = Internal.range2Window(lPointRange[0][1])
                                i1=int(r[0]); j1=int(r[2]); k1=int(r[4])
                                i2=int(r[1]); j2=int(r[3]); k2=int(r[5])
                                wrange = [i1,i2,j1,j2,k1,k2]
                                #   nom de la ZoneBC_t
                                overlapName='overlapBC'+str(c); c=c+1
                                #   doubly_defined
                                userDef = Internal.getNodesFromName(gc[k], 'UserDefinedData')
                                if userDef != []:
                                    if len(userDef[0]) == 4:
                                        info = userDef[0][2][0]
                                        if info[0] == 'doubly_defined':
                                            # Creation d'une famille doubly_defined pour la zone"
                                            familyNameDD='FDD_'+zones[j][0]+'_'+gc[k][0]
                                            ListFDD = (Internal.getNodesFromName1(bases[i],familyNameDD))
                                            if ListFDD == []:
                                                bases[i] = C.addFamily2Base(bases[i], familyNameDD, 'BCOverlap')
                                            zones[j] = C.addBC2Zone(zones[j], overlapName, 'FamilySpecified:'+familyNameDD, wrange, rangeDonor='doubly_defined')
                                            zones[j] = C.tagWithFamily(zones[j],familyNameDD)
                                            FDD = Internal.getNodesFromName1(bases[i],familyNameDD)[0]
                                            (parentBaseDD, numFamilyDD) = Internal.getParentOfNode(bases[i],FDD)
                                            #    Creation d'un noeud doubly_defined
                                            OvlpDD  = Internal.getNodesFromName(FDD,'.Solver#Overlap')[0]
                                            (parentFamilyDD, numOvlpDD) = Internal.getParentOfNode(FDD,OvlpDD)
                                            FDD[2][numOvlpDD][2].append(['NeighbourList', None, [], 'DataArray_t'])
                                            Internal._createChild(FDD[2][numOvlpDD], 'doubly_defined', 'DataArray_t', value='active')
                                            bases[i][2][numFamilyDD] = FDD
                                else:
                                    zones[j] = C.addBC2Zone(zones[j], overlapName, 'FamilySpecified:'+familyName, wrange)
                                    familyNameList = Internal.getNodesFromType1(zones[j],'FamilyName_t')
                                    if familyNameList == []: zones[j] = C.tagWithFamily(zones[j],familyName)

            bases[i][2][numZone] = zones[j]

    cgnsv = Internal.getNodesFromType1(tp, 'CGNSLibraryVersion_t')
    if cgnsv != []:
        (parent, pos) = Internal.getParentOfNode(tp,cgnsv[0])
        bases = [tp[2][pos]]+bases
    tp[2] = bases
    return tp

def addBaseToDonorZone__(t):
    tp = Internal.copyRef(t)

    zones = Internal.getZones(tp)
    for z in zones:
        subRegions= Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        interpSubRegions=[]
        for s in subRegions:
            sname = s[0]
            if sname.split('_')[0] == 'ID': interpSubRegions.append(s)
        for s in interpSubRegions:
            donorname = s[1].tobytes().decode()
            donorzone = Internal.getNodeFromName(tp, donorname)
            base,pos = Internal.getParentOfNode(tp, donorzone)
            Internal._setValue(s, base[0]+"/"+donorname)
    return tp

# Remove useless families
def removeExtraFamily__(t, listOfNodes):
    tp = Internal.copyRef(t)
    bases = Internal.getBases(tp)
    toDelete = []
    for i in range(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for j in range(len(families)):
            familyName = families[j][0]
            path="%s/%s"%(baseName,familyName)
            # Checking if the family is not needed  if it is a default family with just a familyBC_t node of type UserDefined
            if path not in listOfNodes and len(families[j][2]) == 1 and families[j][2][0][2] == [] and families[j][2][0][1].tobytes().decode() == 'UserDefined':
                Internal._rmNode(tp,families[j])
                toDelete.append(path)

        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in range(len(zones)):
            familyNameNodes = Internal.getNodesFromType1(zones[j],'FamilyName_t')
            additionalFamilyNameNodes = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')

            for k in range(len(familyNameNodes)):
                familyName = familyNameNodes[k][1]
                if type(familyName) == 'numpy.ndarray':
                    familyName = familyName.tobytes().decode()
                familyPath = "%s/%s"%(baseName,familyName)
                if familyPath in toDelete:
                    Internal._rmNode(tp, familyNameNodes[k])

            for k in range(len(additionalFamilyNameNodes)):
                familyName = additionalFamilyNameNodes[k][1].tobytes().decode()
                familyPath = "%s/%s"%(baseName,familyName)
                if familyPath in toDelete:
                    Internal._rmNode(tp, additionalFamilyNameNodes[k])
    return tp

# Add new merged family and update family names
# - Add new family names
# - remove the old ones and get the .Solver#Overlap nodes
# - Replace value of familyName_t nodes with the new merged families for Zone_t and BC_t
# - Update NeighbourList node
def addMergedFamily__(t, equivalenceNodes):
    tp = Internal.copyRef(t)
    # Looking for overlap BC FamilySpecified!
    bases = Internal.getBases(tp)
    for i in range(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')

        # Adding new family
        for e in sorted(equivalenceNodes.keys()):
            if e.split("/")[0] == baseName:
                Fnew = [e.split("/")[1], None, [], 'Family_t']
                Fnew[2].append(['FamilyBC', 'UserDefined', [], 'FamilyBC_t'])
                ovlp = None
                for familyName in equivalenceNodes[e]:
                    fname = familyName.split("/")[1]
                    F = Internal.getNodeFromName(bases[i],fname)
                    if F is not None:
                        # Getting the .Solver#Overlap node
                        if len(Internal.getNodesFromName(F, '.Solver#Overlap')):
                            ovlp = Internal.getNodesFromName(F, '.Solver#Overlap')[0]
                    Internal._rmNode(tp, F)
                if ovlp:
                    Fnew[2].append(ovlp)
                bases[i][2].append(Fnew)

        # Updating NeighbourList
        for j in range(len(families)):
            familyName = families[j][0]
            path = "%s/%s"%(baseName,familyName)

            Ovlp  = Internal.getNodesFromName(families[j],'.Solver#Overlap')
            for k in range(len(Ovlp)):
                NeighbourList = Internal.getNodesFromName(Ovlp[k], 'NeighbourList')
                nl =  NeighbourList[0][1]
                if nl is not None:# can be None for doubly defined BCs
                    nl = nl.tobytes().decode().split()
                    for e in sorted(equivalenceNodes.keys()):
                        if equivalenceNodes[e] <= set(nl):
                            familyNodes = Internal.getNodesFromName(bases[i],e)
                            snl = set(nl) - equivalenceNodes[e]
                            snl.add(e)
                            Internal._setValue(NeighbourList[0], " ".join(list(snl)))

        zones = Internal.getNodesFromType1(bases[i], 'Zone_t')
        for j in range(len(zones)):
            zoneName = zones[j][0]
            path = "%s/%s"%(baseName,zoneName)
            familyNameNodes = Internal.getNodesFromType1(zones[j],'FamilyName_t')
            additionalFamilyNameNodes = Internal.getNodesFromType1(zones[j],'AdditionalFamilyName_t')

            # Updating FamilyName_t for Zone_t
            for k in range(len(familyNameNodes)):
                familyName = familyNameNodes[k][1]
                if type(familyName) == 'numpy.ndarray':
                    familyName = familyName.tobytes().decode()
                familyPath = "%s/%s"%(baseName,familyName)
                for e in sorted(equivalenceNodes.keys()):
                    if familyPath in equivalenceNodes[e]:
                        Internal._setValue(familyNameNodes[k], e.split("/")[1])

            # Updating AdditionalFamilyName_t for Zone_t
            for k in range(len(additionalFamilyNameNodes)):
                familyName = additionalFamilyNameNodes[k][1].tobytes().decode()
                familyPath = "%s/%s"%(baseName,familyName)
                for e in sorted(equivalenceNodes.keys()):
                    if familyPath in equivalenceNodes[e]:
                        Internal._setValue(additionalFamilyNameNodes[k], e.split("/")[1])


            # Updating familyName_t for familySpecified BC_t
            zonebcs = Internal.getNodesFromType1(zones[j], 'ZoneBC_t')
            for k in range(len(zonebcs)):
                bcs = Internal.getNodesFromType1(zonebcs[k], 'BC_t')
                for l in range(len(bcs)):
                    if bcs[l][1].tobytes().decode() == 'FamilySpecified':
                        familyNameNodes = Internal.getNodesFromType1(bcs[l], 'FamilyName_t')
                        additionalFamilyNameNodes = Internal.getNodesFromType1(bcs[l], 'AdditionalFamilyName_t')
                        familyName = familyNameNodes[0][1].tobytes().decode()
                        familyPath = "%s/%s"%(baseName,familyName)
                        for e in sorted(equivalenceNodes.keys()):
                            if familyPath in equivalenceNodes[e]:
                                Internal._setValue(familyNameNodes[0], e.split("/")[1])

                        # Updating AdditionalFamilyName_t for Zone_t
                        for k in range(len(additionalFamilyNameNodes)):
                            familyName = additionalFamilyNameNodes[k][1].tobytes().decode()
                            familyPath = "%s/%s"%(baseName,familyName)
                            for e in sorted(equivalenceNodes.keys()):
                                if familyPath in equivalenceNodes[e]:
                                    Internal._setValue(additionalFamilyNameNodes[k], e.split("/")[1])


    return tp


# Building graph of neighbours zones
def buildGraph__(t):
    bases = Internal.getBases(t)
    g = {}
    for i in range(len(bases)):
        baseName = bases[i][0]
        families = Internal.getNodesFromType1(bases[i], 'Family_t')
        for j in range(len(families)):
            familyName = families[j][0]
            path = "%s/%s"%(baseName, familyName)
            Ovlp = Internal.getNodesFromName(families[j],'.Solver#Overlap')
            for k in range(len(Ovlp)):
                if path not in g:
                    g[path] = set()
                NeighbourList = Internal.getNodesFromName(Ovlp[k], 'NeighbourList')
                nl = NeighbourList[0][1]
                if nl is not None: # can be None for doubly defined BCs
                    nl = nl.tobytes().decode().split()
                    g[path].update(nl)
                    for elt in nl:
                        if elt not in g: g[elt] = set()
                        g[elt].add(path)
    return g

# Compute equivalent families that can be merged.
# Families are equivalent and can be merged if they have the same neighbours
def buildPart__(g):
    equivalenceValues = {}
    equivalenceNodes = {}
    i = 0
    for k1 in sorted(g.keys()):
        for k2 in sorted(g.keys()):
            if k1 < k2:
                if g[k1] == g[k2]:
                    if g[k1] not in equivalenceValues.values():
                        i = i+1
                        equivalenceValues[i] = g[k1]
                    if i not in equivalenceNodes:
                        equivalenceNodes[i] = set([])
                    equivalenceNodes[i].update([k1,k2])

    equivalenceNames = {}
    i = 0
    for e in equivalenceNodes:
        i += 1
        baseName = list(equivalenceNodes[e])[0].split("/")[0]
        newFamilyName = "F_Merged.%d"%(i)
        newPath = "%s/%s"%(baseName,newFamilyName)
        equivalenceNames[newPath] = equivalenceNodes[e]

    return equivalenceNames

# Nettoie l'arbre
def cleanTree__(t):
    tp = Internal.copyRef(t)
    g = buildGraph__(tp)
    equivalenceNodes = buildPart__(g)
    neededFamily = g.keys()
    # Remove useless families
    tp = removeExtraFamily__(tp, g.keys())
    # Useless call!
    # Add new merged family and update family names
    tp = addMergedFamily__(tp, equivalenceNodes)
    return tp


# Split point list of hybrid matching BC on structured grids
# one BC match per structured face
def _splitHybridBCMatch(t):
    bases = Internal.getBases(t)
    for base in bases:
        zones = Internal.getZones(base)
        for z in zones:
            dims = Internal.getZoneDim(z)
            if dims[0] == 'Structured':
                zgc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
                if zgc is not None:
                    gc = Internal.getNodesFromType1(zgc, 'GridConnectivity1to1_t')
                    for g in gc:
                        PL = Internal.getNodeFromName1(g, 'PointList')
                        PLD = Internal.getNodeFromName1(g, 'PointListDonor')
                        if PL is not None and PLD is not None:
                            SPL = Converter.converter.pointList2SPL(PL[1],PLD[1],dims[1],dims[2],dims[3])
                            name = g[0]
                            Internal._rmNodesFromName(zgc, name)
                            for i in range(6):
                                if SPL[i] is not None:
                                    n = Internal.createChild(zgc, C.getBCName(name), 'GridConnectivity1to1_t', value=Internal.getValue(g))
                                    Internal.createChild(n, 'GridLocation', 'GridLocation_t', value='FaceCenter')
                                    Internal.createChild(n, 'PointList', 'IndexArray_t', value=SPL[i].reshape(1,SPL[i].size))
                                    Internal.createChild(n, 'PointListDonor', 'IndexArray_t', value=SPL[i+6].reshape(1,SPL[i+6].size))
                                    Internal.createChild(n, 'GridConnectivityType', 'GridConnectivityType_t', value='Abutting1to1')

            if dims[0] == 'Unstructured':
                zgc = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
                if zgc is not None:
                    gc = Internal.getNodesFromType1(zgc, 'GridConnectivity_t')
                    for g in gc:
                        zdname = Internal.getValue(g)
                        zd = Internal.getNodeFromName1(base,zdname)
                        dimsd = Internal.getZoneDim(zd)
                        if dimsd[0] == 'Structured':
                            PL = Internal.getNodeFromName1(g, 'PointList')
                            PLD = Internal.getNodeFromName1(g, 'PointListDonor')
                            if PL is not None and PLD is not None:
                                SPL = Converter.converter.pointList2SPL(PLD[1],PL[1],dimsd[1],dimsd[2],dimsd[3])
                                name = g[0]
                                Internal._rmNodesFromName(zgc, name)
                                for i in range(6):
                                    if SPL[i] is not None:
                                        n = Internal.createChild(zgc, C.getBCName(name), 'GridConnectivity_t', value=zdname)
                                        Internal.createChild(n, 'GridLocation', 'GridLocation_t', value='FaceCenter')
                                        Internal.createChild(n, 'PointList', 'IndexArray_t', value=SPL[i+6].reshape(1,SPL[i+6].size))
                                        Internal.createChild(n, 'PointListDonor', 'IndexArray_t', value=SPL[i].reshape(1,SPL[i].size))
                                        Internal.createChild(n, 'GridConnectivityType', 'GridConnectivityType_t', value='Abutting1to1')
    return None


# convert NearMatch elsA to NoMatch elsA
def _convertNearMatch2NoMatch(t):
    for zone in Internal.getZones(t):
        for match in Internal.getNodesFromType2(zone,'GridConnectivity1to1_t'):
            zdname = Internal.getValue(match)
            jtnode = Internal.getNodeFromName2(match,'jtype')
            if jtnode is not None:
                Internal.setType(match,'GridConnectivity_t')
                Internal.createChild(match,'GridConnectivityType','GridConnectivityType_t',value='Abutting')
                Internal._rmNodesByName(match,'PointRangeDonor')
                Internal._rmNodesByName(match,'Transform')
                Internal.setValue(jtnode,'nomatch')
                prop = Internal.getNodeFromName1(match,'.Solver#Property')
                Internal.createChild(prop,'globborder','DataArray_t',value='Glob%s_%s'%(zone[0],zdname))
                Internal.createChild(prop,'globborderdonor','DataArray_t',value='Glob%s_%s'%(zdname,zone[0]))
                Internal.createChild(prop,'nomatch_special','DataArray_t',value='none')
                Internal._rmNodesByName(prop,'matchside')
                Internal._rmNodesByName(prop,'i_ratio')
                Internal._rmNodesByName(prop,'j_ratio')
                Internal._rmNodesByName(prop,'k_ratio')

    return None



#==============================================================================
# VERSION OBSOLETE
#==============================================================================
def convert2elsAxdt__(t, sameBase=0):

    tp = Internal.copyRef(t)

    _addTurbulentDistanceIndex(tp)
    _buildMaskFiles(tp)
    _adaptNearMatch(tp)
    _adaptPeriodicMatch(tp)

    tp = buildBCOverlap__(tp)
    tp = addNeighbours__(tp, sameBase)
    tp = addFamilyBCNode__(tp)
    tp = cleanTree__(tp)
    tp = addBaseToDonorZone__(tp)
    return tp
