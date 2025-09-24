"""Immersed boundary geometry definition module.
"""
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy
import copy
import math
from . import PyTree as D

varsDeleteIBM = ['utau','StagnationEnthalpy','StagnationPressure',
                 'dirx'          ,'diry'          ,'dirz',
                 'gradxPressure' ,'gradyPressure' ,'gradzPressure' ,
                 'gradxVelocityX','gradyVelocityX','gradzVelocityX',
                 'gradxVelocityY','gradyVelocityY','gradzVelocityY',
                 'gradxVelocityZ','gradyVelocityZ','gradzVelocityZ',
                 'KCurv'         ,'yplus'         ,
                 't11_model'     ,'t12_model'     ,'t22_model',
                 't13_model'     ,'t23_model'     ,'t33_model']
varsDeleteIBMRotTmp=['CoordinateX_PC#Init','CoordinateX_PC#Init','CoordinateX_PC#Init',
                     'CoordinateX_PW#Init','CoordinateX_PW#Init','CoordinateX_PW#Init',
                     'CoordinateX_PI#Init','CoordinateX_PI#Init','CoordinateX_PI#Init',
                     'MotionType','omega',
                     'transl_speedX','transl_speedY','transl_speedZ',
                     'axis_pntX'    ,'axis_pntY'    ,'axis_pntZ'    ,
                     'axis_vctX'    ,'axis_vctY'    ,'axis_vctZ'    ]

EPSCART = 1.e-6

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## COMPUTE INFO FOR F42 (e.g. Yplus & modelisation height etc.)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#=============================================================================
# Compute the skin friction coefficient for a given emperical law
#=============================================================================
def compute_Cf(Re, Cf_law='ANSYS'):
    if Cf_law == 'ANSYS':
        return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        return 1.328*Re**(-0.5)


#=============================================================================
# Compute the corresponding yplus of a given modeling height
#=============================================================================
def computeYplus(Re, Cf_law='ANSYS', height=0.1, L=1.):
    h0 = (L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))
    return height/h0


#=============================================================================
# Compute the modeling height
#=============================================================================
def computeModelisationHeight(Re, Cf_law='ANSYS', yplus=100., L=1.):
    return (yplus*L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))


#=============================================================================
# Compute the best modeling height for a given snear
#=============================================================================
def computeBestModelisationHeight(Re, h, Cf_law='ANSYS', L=1., q=1.2):
    h0 = (L*numpy.sqrt(2))/(Re*numpy.sqrt(compute_Cf(Re,Cf_law)))
    hmod = (h0-q*h)/(1.-q)
    return hmod, hmod/h0


def computeYplusOpt(Re=None,tb=None,Lref=1.,q=1.2,snear=None,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else:
                Re = Internal.getValue(Re)
        else: fail = 1
    if fail:
        raise ValueError("computeYplusOpt: requires Reynolds number as a float or in tb.")
    fail = 0
    if snear is None:
        snear = Internal.getNodeFromName(tb,"snear")
        if snear is None: fail=1
        else: snear = Internal.getValue(snear)
    if fail:
        raise ValueError("computeYlusOpt: requires snear as a float or in tb.")

    print("Warning: estimation of the optimum y+ at Reynolds number ", Re, " and snear target at image point ", snear)
    h0 = (1.*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+1
    h_opti = (h0-q*snear)/(1.-q) #Hauteur de modelisation opti
    yplus_opti = h_opti/h0 #yplus opti
    # print('\nInformation for the body-fitted mesh :')
    # print('h_opti     = %1.2e'%(h_opti))
    # print('h0         = %1.2e\n'%(h0))
    # print('Information for the Cartesian mesh :')
    # print('yplus_opti = %d\n'%(int(math.ceil(yplus_opti))))
    return yplus_opti


# compute the near wall spacing in agreement with the yplus target at image points - front42
def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    fail=0
    if Re is None:
        if tb is not None:
            Re = Internal.getNodeFromName(tb,"Reynolds")
            if Re is None: fail=1
            else: Re = Internal.getValue(Re)
        else: fail = 1
    if fail:
        raise ValueError("computeSnearOpt: requires Reynolds number as a float or in tb.")


    print("Estimation of the optimum near-wall spacing at Reynolds number ", Re, " and yplus target at image point ", yplus)
    h_mod = (yplus*Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law)))
    h0    = (Lref*math.sqrt(2.))/(Re*math.sqrt(compute_Cf(Re,Cf_law))) #Taille de maille pour y+=1
    n     = int(math.ceil(math.log(1-yplus*(1-q))/math.log(q))) # number of cells in the BF mesh for the height h
    snear_opti = q**(n-1)*h0 # best snear for the target yplus
    # print('\nInformation for the body-fitted mesh :')
    # print('h           = %1.2e'%(h_mod))
    # print('h0          = %1.2e\n'%(h0))
    # print('Information for the Cartesian mesh :')
    # print('snear_opti  = %1.3e\n'%(snear_opti))
    return snear_opti


def getMinimumCartesianSpacing(t):
    baseC = Internal.getNodeFromName1(t, 'CARTESIAN')
    if baseC is None: return -1.

    zonesC = Internal.getZones(baseC)
    dxmin = 1.e6
    for z in zonesC:
        dx = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
        if dx < dxmin: dxmin = dx

    print('Minimum spacing on Cartesian grids = %f.'%dxmin, flush=True)
    return dxmin

#==============================================================================
# Creation of a case with a symmetry plane
#==============================================================================
def _symetrizePb(t, bodySymName, snear_sym, dir_sym=2):
    import Converter.Mpi as Cmpi
    if dir_sym not in [1,2,3]: raise ValueError('The symmetry direction %d is not supported. Must be 1(x), 2(y), or 3(z). Exiting...'%dir_sym)
    base   = Internal.getNodeFromName(t, bodySymName)
    minval = C.getMinValue(base, ['CoordinateX', 'CoordinateY','CoordinateZ'])
    minval = minval[dir_sym-1]
    if dir_sym==1: symPlane=(minval,0,0)
    elif dir_sym==2: symPlane=(0,minval,0)
    else: symPlane=(0,0,minval)
    _symetrizeBody(base, dir_sym=dir_sym, symPlane=symPlane)
    _addSymPlane(t, snear_sym, dir_sym=dir_sym, midPlane=minval)
    return None

def _symetrizeBody(base, dir_sym=2, symPlane=(0.,0.,0.)):
    import Transform.PyTree as T
    zones = Internal.getZones(base)
    C._initVars(zones,'centers:cellN',1.)
    if dir_sym == 2:
        v1 = (1,0,0); v2=(0,0,1)
    elif dir_sym == 1:
        v1 = (0,1,0); v2=(0,0,1)
    elif dir_sym == 3:
        v1 = (1,0,0); v2=(0,1,0)

    zones_dup = T.symetrize(zones,symPlane,v1, v2)
    for z in zones_dup: z[0] += '_sym'
    T._reorder(zones_dup,(-1,2,3))
    C._initVars(zones_dup,'centers:cellN',0.)
    base[2]+= zones_dup
    return None

def _addSymPlane(tb, snear_sym, dir_sym=2, midPlane=0):
    import Generator.PyTree as G
    snearList=[]; dfarList=[]
    snearFactor=50
    bodies = Internal.getZones(tb)
    for nos, s in enumerate(bodies):
        sdd = Internal.getNodeFromName1(s, ".Solver#define")
        snearl = Internal.getNodeFromName1(sdd, "snear")
        if snearl is not None:
            snearl = Internal.getValue(snearl)
            snearList.append(snearl*snearFactor)

        dfarl = Internal.getNodeFromName1(sdd,"dfar")
        if dfarl is not None:
            dfarl = Internal.getValue(dfarl)
            dfarList.append(dfarl)

    o = G.octree(tb, snearList=snearList, dfarList=dfarList)

    [xmin,ymin,zmin,xmax,ymax,zmax] = G.bbox(o)
    L = 0.5*(xmax+xmin); eps = 0.2*L
    xmin = xmin-eps; ymin = ymin-eps; zmin = zmin-eps
    xmax = xmax+eps; ymax = ymax+eps; zmax = zmax+eps
    if dir_sym==1: coordsym = 'CoordinateX'; xmax=midPlane
    elif dir_sym==2: coordsym = 'CoordinateY'; ymax=midPlane
    elif dir_sym==3: coordsym = 'CoordinateZ'; zmax=midPlane
    a = D.box((xmin,ymin,zmin),(xmax,ymax,zmax))
    C._initVars(a,'{centers:cellN}=({centers:%s}>(%g-1e-8))'%(coordsym,midPlane))
    C._addBase2PyTree(tb,"SYM")
    base = Internal.getNodeFromName(tb,"SYM"); base[2]+=a
    _setSnear(base,snear_sym)
    _setDfar(base,-1.)
    _setIBCType(base, "slip")
    return None

#==============================================================================
# Set snear in zones
#==============================================================================
def setSnear(t, value):
    """Set the value of snear in a geometry tree.
    Usage: setSnear(t, value=X)"""
    tp = Internal.copyRef(t)
    _setSnear(tp, value)
    return tp

def _setSnear(t, value):
    """Set the value of snear in a geometry tree.
    Usage: _setSnear(t,value=X)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', float(value))
    return None

#==============================================================================
# Set dfar in zones
#==============================================================================
def setDfar(t, value):
    """Set the value of dfar in a geometry tree.
    Usage: setDfar(t, value=X)"""
    tp = Internal.copyRef(t)
    _setDfar(tp, value)
    return tp


def _setDfar(t, value):
    """Set the value of dfar in a geometry tree.
        Usage: _setDfar(t, value=X)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', float(value))
    return None

#==============================================================================
#
#==============================================================================
def getDfarOpt(tb, vmin, snear_opt, factor=10, nlevel=-1):
    """Computes the optimal dfar to get the exact snear.
    Usage: getDfarOpt(tb, vmin, snear_opt, factor=10, nlevel=-1)"""
    import Generator.PyTree as G
    bb = G.bbox(tb)
    sizemax = max(bb[3]-bb[0],bb[4]-bb[1],bb[5]-bb[2])
    #print('sizemax',sizemax)
    if factor>0 and nlevel <1:
        dfar= factor*sizemax
        nlevel = int(numpy.ceil((numpy.log(2*dfar+sizemax)-numpy.log((vmin-1)*snear_opt))/numpy.log(2)))
        #print('nb levels=',nlevel)
    elif nlevel > 0:
        pass
    dfaropt = 0.5*(snear_opt*(2**nlevel*(vmin-1)) - sizemax)
    dfaropt = numpy.trunc(dfaropt*10**nlevel)/10**nlevel
    #print('dfar opt',dfaropt)
    #print('2*dfar opt + sizemax =',2*dfar + sizemax)
    #print('2**n*vmin*snear opt =',(vmin-1)*snear_opt*2**nlevel)
    return dfaropt

#==============================================================================
# Multiply the snear by factors XX in zones
#==============================================================================
def snearFactor(t, sfactor):
    """Multiply the value of snear in a geometry tree by a sfactor.
    Usage: snearFactor(t, sfactor)"""
    tp = Internal.copyRef(t)
    _snearFactor(tp, sfactor)
    return tp


def _snearFactor(t, sfactor):
    """Multiply the value of snear in a geometry tree by a sfactor.
    Usage: _snearFactor(t, sfactor)"""
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, sfactor*Internal.getValue(n))
    return None

#==============================================================================
# Set the IBC type in zones
#==============================================================================
def setIBCType(t, value):
    """Set the IBC type in a geometry tree.
    Usage: setIBCType(t, value=X)"""
    tp = Internal.copyRef(t)
    _setIBCType(tp, value)
    return tp


def _setIBCType(t, value):
    """Set the IBC type in a geometry tree.
    Usage: _setIBCType(t, value=X)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None

#==============================================================================
# Set the fluid inside the geometry
#==============================================================================
def setFluidInside(t):
    """Set fluid inside a geometry tree.
    Usage: setFluidInside(t)"""
    tp = Internal.copyRef(t)
    _setFluidInside(tp)
    return tp

def _setFluidInside(t):
    """Set fluid inside a geometry tree.
    Usage: _setFluidInside(t)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'inv', 'DataArray_t', value=1)
    return None

#==============================================================================
# Set the fluid outside the geometry
#==============================================================================
def setFluidOutside(t):
    """Set fluid inside a geometry tree.
    Usage: setFluidInside(t)"""
    tp = Internal.copyRef(t)
    _setFluidInside(tp)
    return tp

def _setFluidOutside(t):
    """Set fluid inside a geometry tree.
    Usage: _setFluidInside(t)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'inv', 'DataArray_t', value=0)
    return None


#==============================================================================
# Set outpress control parameters in zones
#==============================================================================
def setOutPressControlParam(t, probeName='pointOutPress', AtestSection=1, AOutPress=1,
                            machTarget=0.1, pStatTarget=1e05, tStatTarget=298.15,lmbd=0.1,
                            cxSupport=0.6, sSupport=0.1, itExtrctPrb=10):
    """Set the user input parameters for the outpress control algorithm.
    Usage: setOutPressControlParam(t, probeName='X', AtestSection=Y, AOutPress=Z, machTarget=XX,pStatTarget=YY,tStatTarget=ZZ,lmbd=XXX,cxSupport=YYY,sSupport=ZZZ,itExtrctPrb=XXXX)"""
    tp = Internal.copyRef(t)
    _setOutPressControlParam(tp, probeName=probeName, AtestSection=AtestSection, AOutPress=AOutPress,
                             machTarget=machTarget, pStatTarget=pStatTarget, tStatTarget=tStatTarget,lmbd=lmbd,
                             cxSupport=cxSupport, sSupport=sSupport, itExtrctPrb=itExtrctPrb)
    return tp


def _setOutPressControlParam(t, probeName='pointOutPress', AtestSection=1, AOutPress=1,
                             machTarget=0.1, pStatTarget=1e05, tStatTarget=298.15,lmbd=0.1,
                             cxSupport=0.6, sSupport=0.1, itExtrctPrb=10):
    """Set the user input parameters for the outpress control algorithm.
    Usage: _setOutPressControlParam(t, probeName='X', AtestSection=Y, AOutPress=Z, machTarget=XX,pStatTarget=YY,tStatTarget=ZZ,lmbd=XXX,cxSupport=YYY,sSupport=ZZZ,itExtrctPrb=XXXX)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'probeName'   , 'DataArray_t', probeName)
        Internal._createUniqueChild(n, 'AtestSection', 'DataArray_t', AtestSection)
        Internal._createUniqueChild(n, 'AOutPress'   , 'DataArray_t', AOutPress)
        Internal._createUniqueChild(n, 'machTarget'  , 'DataArray_t', machTarget)
        Internal._createUniqueChild(n, 'pStatTarget' , 'DataArray_t', pStatTarget)
        Internal._createUniqueChild(n, 'tStatTarget' , 'DataArray_t', tStatTarget)
        Internal._createUniqueChild(n, 'lmbd'        , 'DataArray_t', lmbd)
        Internal._createUniqueChild(n, 'cxSupport'   , 'DataArray_t', cxSupport)
        Internal._createUniqueChild(n, 'sSupport'    , 'DataArray_t', sSupport)
        Internal._createUniqueChild(n, 'itExtrctPrb' , 'DataArray_t', itExtrctPrb)

    return None
#==============================================================================
# Set the IBC type outpress for zones in familyName
#==============================================================================
def initOutflow(tc, familyName, PStatic, InterpolPlane=None, PressureVar=0, isDensityConstant=True):
    """Set the value of static pressure PStatic for the outflow pressure IBC with family name familyName. 
    A plane InterpolPlane may also be provided with only static pressure variable or various variables with static pressure as the PressureVar (e.g. 2nd) variable)"""
    tc2 = Internal.copyRef(tc)
    _initOutflow(tc2, familyName, PStatic, InterpolPlane=InterpolPlane, PressureVar=PressureVar, isDensityConstant=isDensityConstant)
    return tc2

def _initOutflow(tc, familyName, PStatic, InterpolPlane=None, PressureVar=0, isDensityConstant=True):
    """Set the value of the pressure PStatic for the outflow pressure IBC with family name familyName.
    A plane InterpolPlane may also be provided with various variables with static pressure as the PressureVar (e.g. 2nd) variable)"""
    import Post.PyTree as P
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc, 'IBCD_4_*'):
            FamNode = Internal.getNodeFromType1(zsr, 'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName == familyName:
                    stagPNode = Internal.getNodeFromName(zsr, 'Pressure')
                    sizeIBC = numpy.shape(stagPNode[1])
                    if InterpolPlane:
                        print("Zone: %s | ZoneSubRegion: %s"%(zc[0],zsr[0]))
                        x_wall = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                        y_wall = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                        z_wall = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]
                        list_pnts = []
                        for i in range(sizeIBC[0]): list_pnts.append((x_wall[i],y_wall[i],z_wall[i]))
                        val = P.extractPoint(InterpolPlane, list_pnts, 2)
                        val_flat = []
                        for i in range(len(val)): val_flat.append(val[i][PressureVar])
                        stagPNode[1][:] = val_flat[:]
                    else:
                        stagPNode[1][:] = PStatic
                    if not isDensityConstant:
                        dens = Internal.getNodeFromName(zsr, 'Density')
                        dens[1][:] = -dens[1][:]
    return None

#==============================================================================
#
#==============================================================================
def initIsoThermal(tc, familyName, TStatic):
    """Set the value of static temperature TStatic for the wall no slip IBC with family name familyName.
    Usage: initIsoThermal(tc, familyName, TStatic)"""
    tc2 = Internal.copyRef(tc)
    _initIsoThermal(tc2, familyName, TStatic)
    return tc2

def _initIsoThermal(tc, familyName, TStatic):
    """Set the value of static temperature TStatic for the wall no slip IBC with family name familyName.
    Usage: _initIsoThermal(tc, familyName, TStatic)"""
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc, 'IBCD_12_*'):
            FamNode = Internal.getNodeFromType1(zsr, 'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName == familyName:
                    stagPNode = Internal.getNodeFromName(zsr, 'TemperatureWall')
                    sizeIBC = numpy.shape(stagPNode[1])
                    stagPNode[1][:] = TStatic
                    #Internal.setValue(stagPNode,TStatic*numpy.ones(sizeIBC))

                    stagPNode = Internal.getNodeFromName(zsr, 'Temperature')
                    Internal.setValue(stagPNode, TStatic*numpy.ones(sizeIBC))
    return None

#==============================================================================
#
#==============================================================================
def initHeatFlux(tc, familyName, QWall):
    """Set the value of heat flux QWall for the wall no slip IBC with family name familyName.
    Usage: initHeatFlux(tc,familyName, QWall)"""
    tc2 = Internal.copyRef(tc)
    _initHeatFlux(tc2, familyName, QWall)
    return tc2


def _initHeatFlux(tc, familyName, QWall):
    """Set the value of heat flux QWall for the wall no slip IBC with family name familyName.
    Usage: _initHeatFlux(tc,familyName, QWall)"""
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_13_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName == familyName:
                    stagPNode = Internal.getNodeFromName(zsr,'WallHeatFlux')
                    sizeIBC   = numpy.shape(stagPNode[1])
                    stagPNode[1][:]  = QWall
                    #Internal.setValue(stagPNode,QWall*numpy.ones(sizeIBC))

                    stagPNode =  Internal.getNodeFromName(zsr,'Temperature')
                    Internal.setValue(stagPNode,QWall*numpy.ones(sizeIBC))
    return None

#==============================================================================
# Set the IBC type inj for zones in familyName
#==============================================================================
def initInj(tc, familyName, PTot, HTot, injDir=[1.,0.,0.], InterpolPlane=None, PressureVar=0, EnthalpyVar=0):
    """Set the total pressure PTot, total enthalpy HTot, and direction of the flow injDir for the injection IBC with family name familyName.
    A plane InterpolPlane may also be provided with at least the stagnation pressure and stagnation enthalpy variables with the former and latter as the PressureVar (e.g. 2nd) and EnthalpyVar variables, respectively.)
    Usage: initInj(tc, familyName, PTot, HTot, injDir, InterpolPlane, PressureVar, EnthalpyVar)"""
    tc2 = Internal.copyRef(tc)
    _initInj(tc2, familyName, PTot, HTot, injDir, InterpolPlane=InterpolPlane, PressureVar=PressureVar, EnthalpyVar=EnthalpyVar)
    return tc2


def _initInj(tc, familyName, PTot, HTot, injDir=[1.,0.,0.], InterpolPlane=None, PressureVar=0, EnthalpyVar=0):
    """Set the total pressure PTot, total enthalpy HTot, and direction of the flow injDir for the injection IBC with family name familyName.
    A plane InterpolPlane may also be provided with at least the stagnation pressure and stagnation enthalpy variables with the former and latter as the PressureVar (e.g. 2nd) and EnthalpyVar variables, respectively.)
    Usage: _initInj(tc, familyName, PTot, HTot, injDir, InterpolPlane, PressureVar, EnthalpyVar)"""
    import Post.PyTree as P
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc, 'IBCD_5_*'):
            FamNode = Internal.getNodeFromType1(zsr, 'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName == familyName:
                    node_temp = Internal.getNodeFromName(zsr, 'utau')
                    if node_temp is not None: Internal._rmNode(zsr, node_temp)

                    node_temp = Internal.getNodeFromName(zsr, 'yplus')
                    if node_temp is not None: Internal._rmNode(zsr, node_temp)

                    stagPNode = Internal.getNodeFromName(zsr, 'StagnationPressure')
                    stagHNode = Internal.getNodeFromName(zsr, 'StagnationEnthalpy')
                    dirxNode  = Internal.getNodeFromName(zsr, 'dirx')
                    diryNode  = Internal.getNodeFromName(zsr, 'diry')
                    dirzNode  = Internal.getNodeFromName(zsr, 'dirz')
                    sizeIBC   = numpy.shape(stagHNode[1])
                    if InterpolPlane:
                        print("Zone: %s | ZoneSubRegion: %s"%(zc[0],zsr[0]))
                        x_wall = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                        y_wall = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                        z_wall = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]
                        list_pnts=[]
                        for i in range(sizeIBC[0]): list_pnts.append((x_wall[i],y_wall[i],z_wall[i]))
                        val = P.extractPoint(InterpolPlane, list_pnts, 2)
                        val_flatPtot = []
                        val_flatHtot = []
                        for i in range(len(val)):
                            val_flatPtot.append(val[i][PressureVar])
                            val_flatHtot.append(val[i][EnthalpyVar])
                        stagPNode[1][:] = val_flatPtot[:]
                        stagHNode[1][:] = val_flatHtot[:]
                    else:
                        stagPNode[1][:] = PTot
                        stagHNode[1][:] = HTot
                    dirxNode[1][:] = injDir[0]
                    diryNode[1][:] = injDir[1]
                    dirzNode[1][:] = injDir[2]

    return None

#==============================================================================
# Add variables to the IBC
#==============================================================================
def _addVariablesTcIbc(zsr, ibctype, nIBC, nsModel='NSLaminar'):
    Nlength = numpy.zeros((nIBC),numpy.float64)
    if ibctype in [2, 3, 6, 10, 11, 31, 32, 331, 332]:
        zsr[2].append(['utau' , copy.copy(Nlength), [], 'DataArray_t'])
        zsr[2].append(['yplus', copy.copy(Nlength), [], 'DataArray_t'])
        if ibctype in [331, 332] and nsModel=='NSLaminar':
            zsr[2].append(['t11_model' , copy.copy(Nlength), [], 'DataArray_t'])
            zsr[2].append(['t12_model' , copy.copy(Nlength), [], 'DataArray_t'])
            zsr[2].append(['t22_model' , copy.copy(Nlength), [], 'DataArray_t'])
            zsr[2].append(['t13_model' , copy.copy(Nlength), [], 'DataArray_t'])
            zsr[2].append(['t23_model' , copy.copy(Nlength), [], 'DataArray_t'])
            zsr[2].append(['t33_model' , copy.copy(Nlength), [], 'DataArray_t'])

    if ibctype == 5:
        Internal._createChild(zsr, 'StagnationEnthalpy', 'DataArray_t', value=copy.copy(Nlength))
        Internal._createChild(zsr, 'StagnationPressure', 'DataArray_t', value=copy.copy(Nlength))
        Internal._createChild(zsr, 'dirx'              , 'DataArray_t', value=copy.copy(Nlength))
        Internal._createChild(zsr, 'diry'              , 'DataArray_t', value=copy.copy(Nlength))
        Internal._createChild(zsr, 'dirz'              , 'DataArray_t', value=copy.copy(Nlength))

    elif ibctype == 10 or ibctype == 11:
        zsr[2].append(['gradxPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
        zsr[2].append(['gradyPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
        zsr[2].append(['gradzPressure' , copy.copy(Nlength) , [], 'DataArray_t'])

        if ibctype == 11:
            zsr[2].append(['gradxVelocityX' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradyVelocityX' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradzVelocityX' , copy.copy(Nlength) , [], 'DataArray_t'])

            zsr[2].append(['gradxVelocityY' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradyVelocityY' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradzVelocityY' , copy.copy(Nlength) , [], 'DataArray_t'])

            zsr[2].append(['gradxVelocityZ' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradyVelocityZ' , copy.copy(Nlength) , [], 'DataArray_t'])
            zsr[2].append(['gradzVelocityZ' , copy.copy(Nlength) , [], 'DataArray_t'])

    elif ibctype == 100:
        zsr[2].append(["KCurv" , copy.copy(Nlength) , [], 'DataArray_t'])

    return None

#==============================================================================
#
#==============================================================================
def changeIBCType(tc, oldIBCType, newIBCType):
    """Change the IBC type in a connectivity tree from oldIBCType to newIBCType.
    Usage: changeIBCType(tc, oldIBCType, newIBCType)"""
    tcp = Internal.copyRef(tc)
    _changeIBCType(tcp, oldIBCType, newIBCType)
    return tcp


def _changeIBCType(tc, oldIBCType, newIBCType):
    """Change the IBC type in a connectivity tree from oldIBCType to newIBCType.
    Usage: changeIBCType(tc, oldIBCType, newIBCType)"""
    for z in Internal.getZones(tc):
        govEqn     = Internal.getNodeFromName(z, 'GoverningEquations')
        nsModel    = 'NSLaminar'
        if govEqn: nsModel = Internal.getValue(govEqn)
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == "IBCD":
                ibcType = int(nameSubRegion.split("_")[1])
                if ibcType == oldIBCType:
                    zsr[0] = "IBCD_{}_".format(newIBCType)+"_".join(nameSubRegion.split("_")[2:])
                    pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                    nIBC = pressure.shape[0]

                    for var_local in varsDeleteIBM:
                        Internal._rmNodesByName(zsr, var_local)

                    _addVariablesTcIbc(zsr, newIBCType, nIBC, nsModel)

    return None

#==============================================================================
#
#==============================================================================
def transformTc2(tc2):
    """Change the name of the IBM nodes for the second image point.
    Usage: transformTc2(tc2)"""
    tcp = Internal.copyRef(tc2)
    _transformTc2(tcp)
    return tcp


def _transformTc2(tc2):
    for z in Internal.getZones(tc2):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:6] == "2_IBCD":
                ibctype = int(nameSubRegion.split("_")[2])
                zsr[0] = "IBCD_{}_".format(ibctype)+"_".join(nameSubRegion.split("_")[3:])

                pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                nIBC = pressure.shape[0]

                vars_delete = ['Density','VelocityX','VelocityY','VelocityZ']+varsDeleteIBM
                for var_local in vars_delete:
                    Internal._rmNodesByName(zsr, var_local)

                Nlength = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['Density'   , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityX' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityY' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityZ' , copy.copy(Nlength) , [], 'DataArray_t'])
                _addVariablesTcIbc(zsr,ibctype,nIBC)

    return None

#==============================================================================
#
#==============================================================================
def closeContour(contour, N=2):
    """Closes a 1D contour with a line or a set of lines."""
    import Transform.PyTree as T

    contour2 = C.convertArray2Tetra(contour)
    contour2 = T.splitConnexity(contour2)

    points = []
    for z in contour2:
        connectivity = Internal.getNodeFromName(z, 'ElementConnectivity')[1]
        values, counts = numpy.unique(connectivity, return_counts=True)
        ends = values[numpy.argwhere(counts==1).flatten()]
        if len(ends) == 2:
            i_first = ends[0] - 1
            i_last  = ends[1] - 1
            x = Internal.getNodeFromName(z, "CoordinateX")[1]
            y = Internal.getNodeFromName(z, "CoordinateY")[1]
            z = Internal.getNodeFromName(z, "CoordinateZ")[1]
            points.append((x[i_first], y[i_first], z[i_first]))
            points.append((x[i_last], y[i_last], z[i_last]))

    if len(points) == 0: return contour

    while points != []:
        P1 = points[0]
        minDist = 1.e6
        minP = (0,0,0)
        for P2 in points[1:]:
            x1,y1,z1 = P1
            x2,y2,z2 = P2
            locDist = numpy.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
            if locDist < minDist: minP = P2
        line = D.line(P1,minP,N=N)
        line = C.convertArray2Tetra(line)
        contour2 = T.join(contour2, line)
        points.remove(P1)
        points.remove(minP)

    contour2 = C.convertBAR2Struct(contour2)
    return contour2

#==============================================================================
#
#==============================================================================
def closeSurface(surface):
    """Closes a 2D surface with a patch or a set of patches."""

    import Post.PyTree as P
    import Transform.PyTree as T
    import Generator.PyTree as G

    cont = P.exteriorFaces(surface)
    cont = T.splitConnexity(cont)

    if len(cont) == 0: return surface

    print("Info: closeSurface: closing surface {} with {} patch(es)".format(surface[0], len(cont)))

    listOfZones = [surface]

    for c in cont:
        c = T.reorder(c, (1,2,3))
        try:
            p = G.fittingPlaster(c)
            p = G.gapfixer(c, p)
            listOfZones.append(p)
        except:
            print("Info: closeSurface: one gapfixer failed, moving on to the next")

    surface2 = C.newPyTree(['Base', listOfZones])

    return surface2

#================================================================================
# Selection/determination of tb for closed solid & filament
#================================================================================
def determineClosedSolidFilament__(tb):
    ##OUT - filamentBases        :list of names of open geometries
    ##OUT - isFilamentOnly       :boolean if there is only a filament in tb
    ##OUT - isOrthoProjectFirst  :bolean to do orthonormal projection first
    ##OUT - tb                   :tb of solid geometries only
    ##OUT - tbFilament           :tb of filament geometries only

    ## General case where only a closeSolid is in tb
    ## or tb only has a filament
    filamentBases = []
    isFilamentOnly= False

    for b in Internal.getBases(tb):
        if "IBCFil" in b[0]:filamentBases.append(b[0])

    if len(filamentBases) == len(Internal.getBases(tb)):isFilamentOnly=True
    isOrthoProjectFirst = isFilamentOnly

    ## assume isFilamentOnly=True
    tbFilament = Internal.copyTree(tb)

    ## if tb has 1) no filaments or 2) filaments & closed bodies
    if not isFilamentOnly:
        if len(filamentBases)==0:
            tbFilament=None
        else:
            ##tb        : only closed bodies
            ##tbFilament: filament bodies
            tbFilament = []
            for b in filamentBases:
                node_local = Internal.getNodeFromNameAndType(tb, b, 'CGNSBase_t')
                tbFilament.append(node_local)
                Internal._rmNode(tb,node_local)
                isOrthoProjectFirst = True
            tbFilament = C.newPyTree(tbFilament);

    return tb, tbFilament


def localWMMFlags__(tb,tbFilament):
    isFilamentOnly=False
    isWireModel   =False

    if tbFilament:
        if len(Internal.getBases(tbFilament))==len(Internal.getBases(tb)):
            bb1 = Internal.getBases(tb)
            bb2 = Internal.getBases(tbFilament)
            if bb1[0][0]==bb2[0][0]:
                isFilamentOnly=True
        for z in Internal.getZones(tbFilament):
            soldef  = Internal.getNodeFromName(z,'.Solver#define')
            if soldef is not None:
                ibctype = Internal.getNodeFromName(soldef,'ibctype')
                if ibctype is not None:
                    if Internal.getValue(ibctype) == "wiremodel":
                        isWireModel=True
                        break
    return isFilamentOnly,isWireModel


###############
# Special test-cases
###############
def flatPlate(snear=0.001, ibctype='Musker'):
    """Generate an IBM case for the canonical flat plate test-case."""
    Gamma = 1.4
    UInf  = 0.2
    RoInf = 1.
    TInf  = 1.
    PInf  = 1./Gamma
    RoEInf = PInf/(Gamma-1.) + 0.5*RoInf*UInf*UInf
    cvInf  = (RoEInf/RoInf - 0.5*UInf*UInf)/TInf

    # plate
    z_wall = D.line((0.,0.,0.), (2.,0.,0.), int(2./snear)); z_wall[0] = 'wall'
    _setSnear(z_wall, snear)
    _setIBCType(z_wall, ibctype)
    _setDfar(z_wall, 0)

    # sym
    z_sym = D.line((-0.33,0.,0.), (0.,0.,0.), 1000); z_sym[0] = 'sym'
    _setSnear(z_sym, snear)
    _setIBCType(z_sym, 'slip')
    _setDfar(z_sym, 0)

    # farfield on top (treated as sym and push back at y=5m)
    z_far = D.line((-0.33,5.,0.), (2.,5.,0.), 1000); z_far[0] = 'farfield'
    _setSnear(z_far, 16*snear)
    _setIBCType(z_far, 'slip')
    _setDfar(z_far, 0)

    # inflow
    z_inj = D.line((-0.33,0.,0.), (-0.33,5.,0.), 1000); z_inj[0] = 'inflow'
    _setSnear(z_inj, 16*snear)
    _setIBCType(z_inj, 'inj')
    _setDfar(z_inj, 0)
    n = Internal.getNodeFromName(z_inj, '.Solver#define')
    Internal.createUniqueChild(n, 'StagnationPressure', 'DataArray_t', value=PInf*1.02828)
    Internal.createUniqueChild(n, 'StagnationEnthalpy', 'DataArray_t', value=TInf*1.008*Gamma*cvInf)
    Internal.createUniqueChild(n, 'dirx', 'DataArray_t', value=1.)
    Internal.createUniqueChild(n, 'diry', 'DataArray_t', value=0.)
    Internal.createUniqueChild(n, 'dirz', 'DataArray_t', value=0.)
    C._tagWithFamily(z_inj, 'INLET')

    # outflow
    z_out  = D.line((2.,0.,0.), (2.,5.,0.), 1000); z_out[0] = 'outflow'
    _setSnear(z_out, 16*snear)
    _setIBCType(z_out, 'outpress')
    _setDfar(z_out, 0)
    n = Internal.getNodeFromName(z_out, '.Solver#define')
    Internal.createUniqueChild(n, 'pStatic',           'DataArray_t', value=PInf)
    Internal.createUniqueChild(n, 'isDensityConstant', 'DataArray_t', value=0.)
    C._tagWithFamily(z_out, 'OUTLET')

    t = C.newPyTree([z_wall, z_sym, z_far, z_inj, z_out])

    _setFluidInside(t)

    C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=5.e6,\
                MutSMuInf=0.2, TurbLevelInf=0.0001, EquationDimension=2, GoverningEquations='NSTurbulent')

    return t

def bumpInChannel(snear=0.001, ibctype='Musker'):
    """Generate an IBM case for the canonical bump-in-channel test-case."""
    import math

    Gamma = 1.4
    UInf  = 0.2
    RoInf = 1.
    TInf  = 1.
    PInf  = 1./Gamma
    RoEInf  = PInf/(Gamma-1.) + 0.5*RoInf*UInf*UInf
    cvInf = (RoEInf/RoInf - 0.5*UInf*UInf)/TInf

    # bump
    z_wall = D.line((0.,0.,0.), (1.5,0.,0.), int(1.5/snear)); z_wall[0] = 'wall'
    x = Internal.getNodeFromName(z_wall, 'CoordinateX')[1]
    y = Internal.getNodeFromName(z_wall, 'CoordinateY')[1]
    for i, x_loc in enumerate(x):
        if x_loc >= 0.3 and x_loc <= 1.2:
            y_loc = 0.05*(math.sin(math.pi*x_loc/0.9-(math.pi/3.)))**4
            y[i] = y_loc
    _setSnear(z_wall, snear)
    _setIBCType(z_wall, ibctype)
    _setDfar(z_wall, 0)

    # sym
    z_symL = D.line((-25.,0.,0.), (0.,0.,0.), 1000); z_symL[0] = 'sym_left'
    z_symR = D.line((1.5,0.,0.), (26.5,0.,0.), 1000); z_symR[0] = 'sym_right'
    z_symT = D.line((-25.,5.,0.), (26.5,5.,0.), 1000); z_symT[0] = 'sym_top'
    for z in [z_symL, z_symR, z_symT]:
        _setSnear(z, 16*snear)
        _setIBCType(z, 'slip')
        _setDfar(z, 0)

    # inflow
    z_inj  = D.line((-25.,0.,0.), (-25.,5.,0.), 1000); z_inj[0] = 'inflow'
    _setSnear(z_inj, 16*snear)
    _setIBCType(z_inj, 'inj')
    _setDfar(z_inj, 0)
    n = Internal.getNodeFromName(z_inj, '.Solver#define')
    Internal.createUniqueChild(n, 'StagnationPressure', 'DataArray_t', value=PInf*1.02828)
    Internal.createUniqueChild(n, 'StagnationEnthalpy', 'DataArray_t', value=TInf*1.008*Gamma*cvInf)
    Internal.createUniqueChild(n, 'dirx', 'DataArray_t', value=1.)
    Internal.createUniqueChild(n, 'diry', 'DataArray_t', value=0.)
    Internal.createUniqueChild(n, 'dirz', 'DataArray_t', value=0.)
    C._tagWithFamily(z_inj, 'INLET')

    # outflow
    z_out  = D.line((26.5,0.,0.), (26.5,5.,0.), 1000); z_out[0] = 'outflow'
    _setSnear(z_out, 16*snear)
    _setIBCType(z_out, 'outpress')
    _setDfar(z_out, 0)
    n = Internal.getNodeFromName(z_out, '.Solver#define')
    Internal.createUniqueChild(n, 'pStatic',           'DataArray_t', value=PInf)
    Internal.createUniqueChild(n, 'isDensityConstant', 'DataArray_t', value=0.)
    C._tagWithFamily(z_out, 'OUTLET')

    t = C.newPyTree([z_wall, z_symL, z_symR, z_symT, z_inj, z_out])

    _setFluidInside(t)

    C._addState(t, adim='adim1', MInf=0.2, alphaZ=0., alphaY=0., ReInf=3.e6,\
                MutSMuInf=0.2, TurbLevelInf=0.0001, EquationDimension=2, GoverningEquations='NSTurbulent')

    return t

def naca0012(snear=0.001, ibctype='Musker', alpha=0.):
    """Generate an IBM case for the canonical subsonic NACA0012 test-case."""

    z = D.naca('0012', N=int(1./snear))
    t = C.newPyTree(['Base', z])

    _setSnear(t, snear)
    _setIBCType(t, ibctype)
    _setDfar(t, 100)

    C._addState(t, adim='adim1', MInf=0.15, alphaZ=alpha, alphaY=0., ReInf=6.e6,\
                MutSMuInf=0.2, TurbLevelInf=0.0001, EquationDimension=2, GoverningEquations='NSTurbulent')

    return t

#====================================================================================
#Add .Solver#Define with dirx, diry, dirz, & granularity to the base of the tboneover. tboneover is the
#PyTree that defines the region in space wherein a one over n coarsening will be pursued
#during the automatic cartesian grid generator of FastIBC.
#IN: t: PyTree
#IN: oneOver: list of list of dirx,diry,dirz,granularity for each base in tboneover. E.g. oneOver=[[1,1,2,0],[1,2,1,0],[2,1,1,1]]
#             for a tboneover with 3 bases where the 1st base has dirx=1, diry=1, dirz=2, & granularity=0 (coarse)
#                                                    2nd base has dirx=1, diry=2, dirz=1, & granularity=0 (coarse)
#                                                    3rd base has dirx=2, diry=1, dirz=1, & granularity=0 (fine)
#OUT: Nothing. Rewrite tboneover with the same FileName as that original used
##NOTE # 1: To be run SEQUENTIALLY ONLY. This is ok as we are dealing with a surface geometry which tend to be
##          relatively small.
##NOTE # 2: Generation of tboneover is similar to that used for tbox.
def _addOneOverLocally(t,oneOver):
    count   = 0
    nodes   = Internal.getNodesFromNameAndType(t, '*OneOver*', 'CGNSBase_t')
    for b in nodes:
        Internal._createUniqueChild(b, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(b, '.Solver#define')
        Internal._createUniqueChild(n, 'dirx'       , 'DataArray_t', value=oneOver[count][0])
        Internal._createUniqueChild(n, 'diry'       , 'DataArray_t', value=oneOver[count][1])
        Internal._createUniqueChild(n, 'dirz'       , 'DataArray_t', value=oneOver[count][2])
        Internal._createUniqueChild(n, 'granularity', 'DataArray_t', value=oneOver[count][3])
        count+=1
    return None
