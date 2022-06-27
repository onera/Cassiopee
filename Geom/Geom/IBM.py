"""Immersed boundary geometry definition module.
"""
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter
import numpy
import copy

vars_delete_ibm=['utau','StagnationEnthalpy','StagnationPressure',
                 'dirx'          ,'diry'          ,'dirz',
                 'gradxPressure' ,'gradyPressure' ,'gradzPressure' ,
                 'gradxVelocityX','gradyVelocityX','gradzVelocityX',
                 'gradxVelocityY','gradyVelocityY','gradzVelocityY',
                 'gradxVelocityZ','gradyVelocityZ','gradzVelocityZ',
                 'KCurv','yplus']


def check_input_familyName(t,familyName=[]):
    names = C.getFamilyZoneNames(t)
    for x in (x for x in familyName if x not in names):
        print('ERROR::FamilyName not in tree. FamilyNames in tree are: %s. Currently selected familyNames: %s'%(names,familyName))
        exit()
    return None
        

# Set snear in zones
def setSnear(t, value, familyName=[]):
    """Set the value of snear in a geometry tree.
    Usage: setSnear(t,value=X,familyName=[])"""
    tp = Internal.copyRef(t)
    _setSnear(tp, value,familyName=familyName)
    return tp


def _setSnear(t, value, familyName=[]):
    """Set the value of snear in a geometry tree.
    Usage: _setSnear(t,value=X,familyName=[])"""
    check_input_familyName(t,familyName)
    zones = Internal.getZones(t)
    if familyName:
        zones    = C.getFamilyZones(t, familyName)
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


# Multiply the snear by factors XX in zones
def snearFactor(t, sfactor,familyName=[]):
    """Multiply the value of snear in a geometry tree by a sfactor.
    Usage: snearFactor(t,sfactor,familyName=[])"""
    tp = Internal.copyRef(t)
    _snearFactor(tp, sfactor,familyName=familyName)
    return tp


def _snearFactor(t, sfactor, familyName=[]):
    """Multiply the value of snear in a geometry tree by a sfactor.
    Usage: _snearFactor(t,sfactor,familyName=[])"""
    check_input_familyName(t,familyName)
    zones = Internal.getZones(t)
    if familyName:
        zones    = C.getFamilyZones(t, familyName)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, sfactor*Internal.getValue(n))
    return None


# Set the IBC type in zones
def setIBCType(t, value,familyName=[]):
    """Set the IBC type in a geometry tree.
    Usage: setIBCType(t,value=X,familyName=[])"""
    tp = Internal.copyRef(t)
    _setIBCType(tp, value,familyName=familyName)
    return tp


def _setIBCType(t, value,familyName=[]):
    """Set the IBC type in a geometry tree.
    Usage: _setIBCType(t,value=X,familyName=[])"""
    check_input_familyName(t,familyName)
    zones = Internal.getZones(t)
    if familyName:
        zones    = C.getFamilyZones(t, familyName)
    for z in zones:         
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None


# Set the IBC type outpress for zones in familyName
def initOutflow(tc, familyName, P_static):
    """Set the value of static pressure P_static for the outflow pressure IBC with family name familyName.
    Usage: initOutflow(tc,familyName, P_static)"""
    tc2 = Internal.copyRef(tc)
    _initOutflow(tc2, familyName, P_static)
    return tc2
                 

def _initOutflow(tc, familyName, P_static):
    """Set the value of the pressure P_static for the outflow pressure IBC with family name familyName.
    Usave: _initOutflow(tc,familyName, P_static)"""    
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_4_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyName:
                    stagPNode =  Internal.getNodeFromName(zsr,'Pressure')    
                    sizeIBC = numpy.shape(stagPNode[1])
                    Internal.setValue(stagPNode,P_static*numpy.ones(sizeIBC))
    return None


# Set the IBC type inj for zones in familyName
def initInj(tc, familyName, P_tot, H_tot, injDir=[1.,0.,0.]):
    """Set the total pressure P_tot, total enthalpy H_tot, and direction of the flow injDir for the injection IBC with family name familyName.
    Usave: initInj(tc, familyName, P_tot, H_tot, injDir)"""
    tc2 = Internal.copyRef(tc)
    _initInj(tc2, familyName, P_tot, H_tot, injDir)
    return tc2
                 

def _initInj(tc, familyName, P_tot, H_tot, injDir=[1.,0.,0.]):
    """Set the total pressure P_tot, total enthalpy H_tot, and direction of the flow injDir for the injection IBC with family name familyName.
    Usage: _initInj(tc, familyName, P_tot, H_tot, injDir)"""
    for zc in Internal.getZones(tc):
        for zsr in Internal.getNodesFromName(zc,'IBCD_5_*'):
            FamNode = Internal.getNodeFromType1(zsr,'FamilyName_t')
            if FamNode is not None:
                FamName = Internal.getValue(FamNode)
                if FamName==familyName:
                    stagPNode =  Internal.getNodeFromName(zsr,'StagnationPressure')
                    stagHNode =  Internal.getNodeFromName(zsr,'StagnationEnthalpy')
                    dirxNode = Internal.getNodeFromName(zsr,'dirx')
                    diryNode = Internal.getNodeFromName(zsr,'diry')
                    dirzNode = Internal.getNodeFromName(zsr,'dirz')
                    sizeIBC = numpy.shape(stagHNode[1])
                    Internal.setValue(stagHNode,H_tot*numpy.ones(sizeIBC))
                    Internal.setValue(stagPNode,P_tot*numpy.ones(sizeIBC))

                    Internal.setValue(dirxNode, injDir[0]*numpy.ones(sizeIBC))
                    Internal.setValue(diryNode, injDir[1]*numpy.ones(sizeIBC))
                    Internal.setValue(dirzNode, injDir[2]*numpy.ones(sizeIBC))
                    
    return None


# Change IBC Types
def add_variables_tc_ibc(zsr,ibctype,nIBC):
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

    if ibctype == 10 or ibctype == 11:
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
        
    if ibctype == 100:
        zsr[2].append(["KCurv" , copy.copy(Nlength) , [], 'DataArray_t'])
        
    return zsr


def changeIBCType(tc, oldIBCType, newIBCType):
    """Change the IBC type in a connectivity tree from oldIBCType to newIBCType.
    Usave: changeIBCType(tc, oldIBCType, newIBCType)"""
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

                    for var_local in vars_delete_ibm:
                        Internal._rmNodesByName(zsr,var_local)
                    

                    zsr=add_variables_tc_ibc(zsr,newIBCType,nIBC)

    return tc


def transformTc2(tc2):
    """Change the name of the IBM nodes for the second image point.
    Usage: transformTc2(tc2)"""
    for z in Internal.getZones(tc2):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:6] == "2_IBCD":
                ibctype = int(nameSubRegion.split("_")[2])
                zsr[0] = "IBCD_{}_".format(ibctype)+"_".join(nameSubRegion.split("_")[3:])

                pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                nIBC = pressure.shape[0]
                
                vars_delete = ['Density','VelocityX','VelocityY','VelocityZ']+vars_delete_ibm
                for var_local in vars_delete:
                    Internal._rmNodesByName(zsr,vars_delete)

                Nlength = numpy.zeros((nIBC),numpy.float64)
                zsr[2].append(['Density'    , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VeloicityX' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VeloicityY' , copy.copy(Nlength) , [], 'DataArray_t'])
                zsr[2].append(['VeloicityZ' , copy.copy(Nlength) , [], 'DataArray_t'])

                zsr=add_variables_tc_ibc(zsr,ibctype,nIBC)
                
    return tc2


