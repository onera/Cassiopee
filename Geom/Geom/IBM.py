"""Immersed boundary geometry definition module.
"""
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter
import numpy
import copy
import Generator.IBMmodelHeight as G_IBM_Height

varsDeleteIBM=['utau','StagnationEnthalpy','StagnationPressure',
               'dirx'          ,'diry'          ,'dirz',
               'gradxPressure' ,'gradyPressure' ,'gradzPressure' ,
               'gradxVelocityX','gradyVelocityX','gradzVelocityX',
               'gradxVelocityY','gradyVelocityY','gradzVelocityY',
               'gradxVelocityZ','gradyVelocityZ','gradzVelocityZ',
            'KCurv','yplus']        

# compute the near wall spacing in agreement with the yplus target at image points - front42
def computeSnearOpt(Re=None,tb=None,Lref=1.,q=1.2,yplus=300.,Cf_law='ANSYS'):
    return G_IBM_Height.computeSnearOpt(Re=Re, tb=tb, Lref=Lref, q=q, yplus=yplus, Cf_law=Cf_law)

# Set snear in zones
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
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', value)
    return None

# Set dfar in zones
def setDfar(t, value):
    """Set the value of dfar in a geometry tree.
    Usage: setDfar(t,value=X)"""
    tp = Internal.copyRef(t)
    _setDfar(tp, value)
    return tp


def _setDfar(t, value):
    """Set the value of dfar in a geometry tree.
        Usage: _setDfar(t,value=X)"""
    zones = Internal.getZones(t)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', value)
    return None

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

# Multiply the snear by factors XX in zones
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


# Set the IBC type in zones
def setIBCType(t, value):
    """Set the IBC type in a geometry tree.
    Usage: setIBCType(t, value=X)"""
    tp = Internal.copyRef(t)
    _setIBCType(tp, value)
    return tp

def _setIBCType(t, value):
    """Set the IBC type in a geometry tree.
    Usage: _setIBCType(t,value=X)"""
    zones = Internal.getZones(t)
    for z in zones:         
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None

# Set the fluid inside the geometry
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

# Set the IBC type outpress for zones in familyName
def initOutflow(tc, familyName, PStatic):
    """Set the value of static pressure PStatic for the outflow pressure IBC with family name familyName.
    Usage: initOutflow(tc,familyName, PStatic)"""
    tc2 = Internal.copyRef(tc)
    _initOutflow(tc2, familyName, PStatic)
    return tc2

def _initOutflow(tc, familyName, PStatic):
    """Set the value of the pressure PStatic for the outflow pressure IBC with family name familyName.
    Usave: _initOutflow(tc,familyName, PStatic)"""    
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_4_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyName:
                    stagPNode =  Internal.getNodeFromName(zsr,'Pressure')    
                    sizeIBC = numpy.shape(stagPNode[1])
                    Internal.setValue(stagPNode,PStatic*numpy.ones(sizeIBC))
    return None

def initIsoThermal(tc, familyName, TStatic):
    """Set the value of static temperature TStatic for the wall no slip IBC with family name familyName.
    Usage: initIsoThermal(tc,familyName, TStatic)"""
    tc2 = Internal.copyRef(tc)
    _initIsoThermal(tc2, familyName, TStatic)
    return tc2


def _initIsoThermal(tc, familyName, TStatic):
    """Set the value of static temperature TStatic for the wall no slip IBC with family name familyName.
    Usage: _initIsoThermal(tc,familyName, TStatic)"""
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_12_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyName:
                    stagPNode =  Internal.getNodeFromName(zsr,'TemperatureWall')    
                    sizeIBC = numpy.shape(stagPNode[1])
                    Internal.setValue(stagPNode,TStatic*numpy.ones(sizeIBC))

                    stagPNode =  Internal.getNodeFromName(zsr,'Temperature')    
                    Internal.setValue(stagPNode,TStatic*numpy.ones(sizeIBC))
    return None


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
                if FamName==familyName:
                    stagPNode =  Internal.getNodeFromName(zsr,'WallHeatFlux')    
                    sizeIBC = numpy.shape(stagPNode[1])
                    Internal.setValue(stagPNode,QWall*numpy.ones(sizeIBC))

                    stagPNode =  Internal.getNodeFromName(zsr,'Temperature')    
                    Internal.setValue(stagPNode,QWall*numpy.ones(sizeIBC))
    return None




# Set the IBC type inj for zones in familyName
def initInj(tc, familyName, PTot, HTot, injDir=[1.,0.,0.]):
    """Set the total pressure PTot, total enthalpy HTot, and direction of the flow injDir for the injection IBC with family name familyName.
    Usave: initInj(tc, familyName, PTot, HTot, injDir)"""
    tc2 = Internal.copyRef(tc)
    _initInj(tc2, familyName, PTot, HTot, injDir)
    return tc2
                 

def _initInj(tc, familyName, PTot, HTot, injDir=[1.,0.,0.]):
    """Set the total pressure PTot, total enthalpy HTot, and direction of the flow injDir for the injection IBC with family name familyName.
    Usage: _initInj(tc, familyName, PTot, HTot, injDir)"""
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_5_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyName:
                    node_temp = Internal.getNodeFromName(zsr,'utau')
                    if node_temp is not None:Internal._rmNode(zsr,node_temp)

                    node_temp = Internal.getNodeFromName(zsr,'yplus')
                    if node_temp is not None:Internal._rmNode(zsr,node_temp)
                    
                    stagPNode =  Internal.getNodeFromName(zsr,'StagnationPressure')
                    stagHNode =  Internal.getNodeFromName(zsr,'StagnationEnthalpy')
                    dirxNode = Internal.getNodeFromName(zsr,'dirx')
                    diryNode = Internal.getNodeFromName(zsr,'diry')
                    dirzNode = Internal.getNodeFromName(zsr,'dirz')
                    sizeIBC = numpy.shape(stagHNode[1])
                    Internal.setValue(stagHNode,HTot*numpy.ones(sizeIBC))
                    Internal.setValue(stagPNode,PTot*numpy.ones(sizeIBC))

                    Internal.setValue(dirxNode, injDir[0]*numpy.ones(sizeIBC))
                    Internal.setValue(diryNode, injDir[1]*numpy.ones(sizeIBC))
                    Internal.setValue(dirzNode, injDir[2]*numpy.ones(sizeIBC))
                    
    return None


# Add variables to the IBC
def _addVariablesTcIbc(zsr, ibctype, nIBC):
    Nlength = numpy.zeros((nIBC),numpy.float64)
    if ibctype in [2, 3, 6, 10, 11]:
        zsr[2].append(['utau' , copy.copy(Nlength), [], 'DataArray_t'])
        zsr[2].append(['yplus', copy.copy(Nlength), [], 'DataArray_t'])

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
                    
                    _addVariablesTcIbc(zsr, newIBCType, nIBC)

    return None

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
                    Internal._rmNodesByName(zsr,vars_delete)

                Nlength = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['Density'   , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityX' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityY' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VelocityZ' , copy.copy(Nlength) , [], 'DataArray_t'])
                _addVariablesTcIbc(zsr,ibctype,nIBC)
                
    return None


