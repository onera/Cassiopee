# - IBM specific post-processing -
from tkinter.ttk import Progressbar
from . import PyTree as P
import Connector.OversetData as XOD
import Connector.PyTree as X
import Connector.IBM as X_IBM
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

# IMPORTANT NOTE !!
## The functions of this submodule will become decrepit after Jan. 1 2023
#=============================================================================
from . import IBM_OLDIES
#=============================================================================

#=============================================================================
# Extract the flow field at the IBM target points onto the surface.
# IN: tc (tree): connectivity tree containing IBM information
# IN: tb (tree): geometry tree (IBM bodies) 
#       If tb is None, return the cloud of IBM points
#       If tb is not None, the solution is projected onto the bodies, at the vertices
# IN : famZones (list): if famZones is not empty, only extract some subregion families (['FAM1','FAM2,...])
# IN : extractIBMInfo (boolean): if True, store IB coordinates (PC, PI and PW)
#=============================================================================
def extractIBMWallFields(tc, tb=None, coordRef='wall', famZones=[], IBCNames="IBCD_*", extractIBMInfo=False):
    """Extracts the flow field stored at IBM points onto the surface."""
    xwNP = []; ywNP = []; zwNP = []
    xiNP = []; yiNP = []; ziNP = []
    xcNP = []; ycNP = []; zcNP = []
    xNP = []; yNP = []; zNP = []
    pressNP = []; utauNP = []; yplusNP = []; yplusINP = []; densNP = []
    vxNP = []; vyNP = []; vzNP = []
    KCurvNP = []
    gradxPressureNP = []; gradyPressureNP = []; gradzPressureNP = []
    conv1NP = []; conv2NP = []
    temperatureNP=[];
    if coordRef =='cible':coordRef='target'
    dictOfFamilies={}
    if famZones != []:
        out = []
        allIBCD = Internal.getNodesFromName(tc,IBCNames)
        for zsr in allIBCD:
            fam = Internal.getNodeFromType(zsr,'FamilyName_t')
            if fam is not None:
                famName=Internal.getValue(fam)
                if famName in famZones:
                    if famName not in dictOfFamilies: dictOfFamilies[famName]=[zsr]
                    else: dictOfFamilies[famName]+=[zsr]

        # Creation of a single zone
        zsize = numpy.empty((1,3), dtype=Internal.E_NpyInt, order='F')
        zsize[0,0] = 1; zsize[0,1] = 0; zsize[0,2] = 0
        dictOfZoneFamilies={}
        for z in Internal.getZones(tb):
            famName = Internal.getNodeFromType(z, 'FamilyName_t')
            if famName is not None:
                famName = Internal.getValue(famName)
                if famName in famZones:
                    if famName not in dictOfZoneFamilies: dictOfZoneFamilies[famName]=[z]
                    else: dictOfZoneFamilies[famName]+=[z]
        for famName in dictOfFamilies:
            zd = Internal.newZone(name='ZIBC_%s'%famName, zsize=zsize, ztype='Unstructured')
            zd[2] += dictOfFamilies[famName]
            tb2 = None
            if tb is not None:
                zones = dictOfZoneFamilies[famName]
                if zones!=[]: tb2=C.newPyTree(['Base']); tb2[2][1][2]=zones

            zd = extractIBMWallFields(zd, tb=tb2, coordRef=coordRef, famZones=[])
            out += Internal.getZones(zd)
        return out

    else:
        allZSR = Internal.getNodesFromType(tc,'ZoneSubRegion_t')
        if allZSR != []:
            allIBCD = Internal.getNodesFromName(allZSR,IBCNames)
            for IBCD in allIBCD:
                if coordRef == 'target':
                    xPC = Internal.getNodeFromName1(IBCD,"CoordinateX_PC")[1]
                    yPC = Internal.getNodeFromName1(IBCD,"CoordinateY_PC")[1]
                    zPC = Internal.getNodeFromName1(IBCD,"CoordinateZ_PC")[1]
                    xNP.append(xPC); yNP.append(yPC); zNP.append(zPC)
                elif coordRef == 'image':
                    xPI = Internal.getNodeFromName1(IBCD,"CoordinateX_PI")[1]
                    yPI = Internal.getNodeFromName1(IBCD,"CoordinateY_PI")[1]
                    zPI = Internal.getNodeFromName1(IBCD,"CoordinateZ_PI")[1]
                    xNP.append(xPI); yNP.append(yPI); zNP.append(zPI)                 
                else:
                    xPW = Internal.getNodeFromName1(IBCD,"CoordinateX_PW")[1]
                    yPW = Internal.getNodeFromName1(IBCD,"CoordinateY_PW")[1]
                    zPW = Internal.getNodeFromName1(IBCD,"CoordinateZ_PW")[1]
                    xNP.append(xPW); yNP.append(yPW); zNP.append(zPW)

                if extractIBMInfo:
                    if coordRef != 'target':
                        xPC = Internal.getNodeFromName1(IBCD,"CoordinateX_PC")[1]
                        yPC = Internal.getNodeFromName1(IBCD,"CoordinateY_PC")[1]
                        zPC = Internal.getNodeFromName1(IBCD,"CoordinateZ_PC")[1]
                    xcNP.append(xPC); ycNP.append(yPC); zcNP.append(zPC)

                    if coordRef != 'image':
                        xPI = Internal.getNodeFromName1(IBCD,"CoordinateX_PI")[1]
                        yPI = Internal.getNodeFromName1(IBCD,"CoordinateY_PI")[1]
                        zPI = Internal.getNodeFromName1(IBCD,"CoordinateZ_PI")[1]
                    xiNP.append(xPI); yiNP.append(yPI); ziNP.append(zPI)

                    if coordRef != 'wall':
                        xPW = Internal.getNodeFromName1(IBCD,"CoordinateX_PW")[1]
                        yPW = Internal.getNodeFromName1(IBCD,"CoordinateY_PW")[1]
                        zPW = Internal.getNodeFromName1(IBCD,"CoordinateZ_PW")[1]
                    xwNP.append(xPW); ywNP.append(yPW); zwNP.append(zPW)
                    
                PW = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
                if PW is not None: pressNP.append(PW[1])
                RHOW = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
                if RHOW is not None:
                    maxVal = numpy.max(RHOW[1])
                    minVal = numpy.min(RHOW[1])
                    if maxVal<0 and minVal<0: RHOW[1][:]=-RHOW[1][:]
                    densNP.append(RHOW[1])
                UTAUW = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
                if UTAUW is not None: utauNP.append(UTAUW[1])
                YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
                if YPLUSW is not None: yplusNP.append(YPLUSW[1])
                YPLUSIW = Internal.getNodeFromName1(IBCD, "yplusIP")
                if YPLUSIW is not None: yplusINP.append(YPLUSIW[1])

                VXW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
                if VXW is not None: vxNP.append(VXW[1])
                VYW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
                if VYW is not None: vyNP.append(VYW[1])
                VZW = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
                if VZW is not None: vzNP.append(VZW[1])

                KCURVW = Internal.getNodeFromName1(IBCD, XOD.__KCURV__)
                if KCURVW is not None: KCurvNP.append(KCURVW[1])

                TEMP = Internal.getNodeFromName1(IBCD, XOD.__TEMPERATURE__)
                if TEMP is not None:temperatureNP.append(TEMP[1])

                GRADXPW = Internal.getNodeFromName1(IBCD, XOD.__GRADXPRESSURE__)
                if GRADXPW is not None: gradxPressureNP.append(GRADXPW[1])
                GRADYPW = Internal.getNodeFromName1(IBCD, XOD.__GRADYPRESSURE__)
                if GRADYPW is not None: gradyPressureNP.append(GRADYPW[1])
                GRADZPW = Internal.getNodeFromName1(IBCD, XOD.__GRADZPRESSURE__)
                if GRADZPW is not None: gradzPressureNP.append(GRADZPW[1])

                CONV1PW = Internal.getNodeFromName1(IBCD, XOD.__CONV1__)
                if CONV1PW is not None: conv1NP.append(CONV1PW[1])
                CONV2PW = Internal.getNodeFromName1(IBCD, XOD.__CONV2__)
                if CONV2PW is not None: conv2NP.append(CONV2PW[1])
                
    if pressNP == []: return None
    else:
        pressNP = numpy.concatenate(pressNP)
        if densNP != []: densNP = numpy.concatenate(densNP)
        if utauNP != []: utauNP = numpy.concatenate(utauNP)
        if yplusNP != []: yplusNP = numpy.concatenate(yplusNP)
        if yplusINP != []: yplusINP = numpy.concatenate(yplusINP)
        if vxNP != []: vxNP = numpy.concatenate(vxNP)
        if vyNP != []: vyNP = numpy.concatenate(vyNP)
        if vzNP != []: vzNP = numpy.concatenate(vzNP)

        if KCurvNP != []: KCurvNP = numpy.concatenate(KCurvNP)

        if temperatureNP != []: temperatureNP = numpy.concatenate(temperatureNP)
        
        if gradxPressureNP != []: gradxPressureNP = numpy.concatenate(gradxPressureNP)
        if gradyPressureNP != []: gradyPressureNP = numpy.concatenate(gradyPressureNP)
        if gradzPressureNP != []: gradzPressureNP = numpy.concatenate(gradzPressureNP)

        if conv1NP != []: conv1NP = numpy.concatenate(conv1NP)
        if conv2NP != []: conv2NP = numpy.concatenate(conv2NP)

        xNP = numpy.concatenate(xNP)
        yNP = numpy.concatenate(yNP)
        zNP = numpy.concatenate(zNP)

        if extractIBMInfo:
            xcNP = numpy.concatenate(xcNP)
            ycNP = numpy.concatenate(ycNP)
            zcNP = numpy.concatenate(zcNP)

            xiNP = numpy.concatenate(xiNP)
            yiNP = numpy.concatenate(yiNP)
            ziNP = numpy.concatenate(ziNP)

            xwNP = numpy.concatenate(xwNP)
            ywNP = numpy.concatenate(ywNP)
            zwNP = numpy.concatenate(zwNP)

    # Creation d une seule zone
    zsize = numpy.empty((1,3), dtype=Internal.E_NpyInt, order='F')
    zsize[0,0] = xNP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=z)
    coordx = ['CoordinateX',xNP,[],'DataArray_t']
    coordy = ['CoordinateY',yNP,[],'DataArray_t']
    coordz = ['CoordinateZ',zNP,[],'DataArray_t']
    gc[2] = [coordx,coordy,coordz]
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                   gridLocation='Vertex', parent=z)
    FSN[2].append([XOD.__PRESSURE__,pressNP, [],'DataArray_t'])
    FSN[2].append([XOD.__DENSITY__,densNP, [],'DataArray_t'])
  
    if extractIBMInfo:
        FSN[2].append(["CoordinateX_PW",xwNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateY_PW",ywNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateZ_PW",zwNP, [],'DataArray_t'])
        
        FSN[2].append(["CoordinateX_PI",xiNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateY_PI",yiNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateZ_PI",ziNP, [],'DataArray_t'])

        FSN[2].append(["CoordinateX_PC",xcNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateY_PC",ycNP, [],'DataArray_t'])
        FSN[2].append(["CoordinateZ_PC",zcNP, [],'DataArray_t'])
    
    utauPresent = 0; yplusPresent = 0; yplusIPresent = 0
    if utauNP != []:
        utauPresent = 1
        FSN[2].append([XOD.__UTAU__,utauNP, [],'DataArray_t'])
    if yplusNP != []:
        yplusPresent = 1
        FSN[2].append([XOD.__YPLUS__,yplusNP, [],'DataArray_t'])
    if yplusINP != []:
        yplusIPresent = 1
        FSN[2].append(["yplusIP",yplusINP, [],'DataArray_t'])
    vxPresent = 0
    if vxNP != []:
        vxPresent = 1
        FSN[2].append([XOD.__VELOCITYX__,vxNP, [],'DataArray_t'])
        FSN[2].append([XOD.__VELOCITYY__,vyNP, [],'DataArray_t'])
        FSN[2].append([XOD.__VELOCITYZ__,vzNP, [],'DataArray_t'])

    kcurvPresent = 0
    if KCurvNP != []:
        kcurvPresent = 1
        FSN[2].append([XOD.__KCURV__,KCurvNP, [],'DataArray_t'])

    temperaturePresent = 0
    if temperatureNP != []:
        temperaturePresent = 1
        FSN[2].append([XOD.__TEMPERATURE__,temperatureNP, [],'DataArray_t'])

    gradxPressurePresent = 0
    if gradxPressureNP != []:
        gradxPressurePresent = 1
        FSN[2].append([XOD.__GRADXPRESSURE__,gradxPressureNP, [],'DataArray_t'])
        FSN[2].append([XOD.__GRADYPRESSURE__,gradyPressureNP, [],'DataArray_t'])
        FSN[2].append([XOD.__GRADZPRESSURE__,gradzPressureNP, [],'DataArray_t'])

    conv1Present = 0
    if conv1NP != []:
        conv1Present = 1
        FSN[2].append([XOD.__CONV1__,conv1NP, [],'DataArray_t'])
        FSN[2].append([XOD.__CONV2__,conv2NP, [],'DataArray_t'])

    if tb is None: return z
    else:
        dimPb = Internal.getNodeFromName(tb,'EquationDimension')
        if dimPb is None:
            print('Warning: extractIBMWallFields: pb dimension is set to 3.')
            dimPb = 3
        else:
            dimPb = Internal.getValue(dimPb)
        td = Internal.copyRef(tb)
        # Force toutes les zones dans une seule base
        #td = C.newPyTree(['Wall',Internal.getZones(tb)])
        for nob in range(len(td[2])):
            b = td[2][nob]
            if b[3] == 'CGNSBase_t':
                zones = Internal.getNodesFromType1(b, 'Zone_t')
                if zones != []:
                    zones = C.convertArray2Tetra(zones)
                    zones = T.join(zones); zones = G.close(zones)
                    b[2] = [zones]

        if extractIBMInfo:
            C._initVars(td,"CoordinateX_PW",0.)
            C._initVars(td,"CoordinateY_PW",0.)
            C._initVars(td,"CoordinateZ_PW",0.)

            C._initVars(td,"CoordinateX_PI",0.)
            C._initVars(td,"CoordinateY_PI",0.)
            C._initVars(td,"CoordinateZ_PI",0.)

            C._initVars(td,"CoordinateX_PC",0.)
            C._initVars(td,"CoordinateY_PC",0.)
            C._initVars(td,"CoordinateZ_PC",0.)

        C._initVars(td,XOD.__PRESSURE__,0.)
        C._initVars(td,XOD.__DENSITY__,0.)
        if utauPresent==1: C._initVars(td,XOD.__UTAU__,0.)
        if yplusPresent==1: C._initVars(td,XOD.__YPLUS__,0.)
        if yplusIPresent==1: C._initVars(td,"yplusIP",0.)
        if vxPresent==1:
            C._initVars(td,XOD.__VELOCITYX__,0.)
            C._initVars(td,XOD.__VELOCITYY__,0.)
            C._initVars(td,XOD.__VELOCITYZ__,0.)
        if kcurvPresent==1:
            C._initVars(td,XOD.__KCURV__,0.)
        if temperaturePresent==1:
            C._initVars(td,XOD.__TEMPERATURE__,0.)            
        if gradxPressurePresent==1:
            C._initVars(td,XOD.__GRADXPRESSURE__,0.)
            C._initVars(td,XOD.__GRADYPRESSURE__,0.)
            C._initVars(td,XOD.__GRADZPRESSURE__,0.)
        if conv1Present==1:
            C._initVars(td,XOD.__CONV1__,0.)
            C._initVars(td,XOD.__CONV2__,0.)
        #print("projectCloudSolution for dim {}".format(dimPb))

        P._projectCloudSolution(z, td, dim=dimPb)

        return td

#=============================================================================
# compute the shear stress using utau
#=============================================================================
def extractShearStress(ts):
    """Computes the shear stress on the surface."""
    ts2 = Internal.copyRef(ts)
    _extractShearStress(ts2)
    return ts2

def _extractShearStress(ts):
    """Computes the shear stress on the surface.""" 
    if Internal.getNodeFromName(ts, Internal.__FlowSolutionCenters__) is None:
        ts = C.node2Center(ts, Internal.__FlowSolutionNodes__)

    # normal vector
    G._getNormalMap(ts)
    C._normalize(ts, ['centers:sx','centers:sy','centers:sz'])

    # Calc tangent vector
    C._initVars(ts, '{centers:tx}={centers:VelocityX}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sx}')
    C._initVars(ts, '{centers:ty}={centers:VelocityY}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sy}')
    C._initVars(ts, '{centers:tz}={centers:VelocityZ}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sz}')
    C._normalize(ts, ['centers:tx','centers:ty','centers:tz'])

    C._initVars(ts, '{centers:tau_wall}={centers:Density}*{centers:utau}*{centers:utau}')
    C._initVars(ts, '{centers:ShearStressXX}=2*{centers:tau_wall}*{centers:tx}*{centers:sx}')
    C._initVars(ts, '{centers:ShearStressYY}=2*{centers:tau_wall}*{centers:ty}*{centers:sy}')
    C._initVars(ts, '{centers:ShearStressZZ}=2*{centers:tau_wall}*{centers:tz}*{centers:sz}')
    C._initVars(ts, '{centers:ShearStressXY}={centers:tau_wall}*({centers:tx}*{centers:sy}+{centers:ty}*{centers:sx})')
    C._initVars(ts, '{centers:ShearStressXZ}={centers:tau_wall}*({centers:tx}*{centers:sz}+{centers:tz}*{centers:sx})')
    C._initVars(ts, '{centers:ShearStressYZ}={centers:tau_wall}*({centers:ty}*{centers:sz}+{centers:tz}*{centers:sy})')
    C._rmVars(ts,['centers:utau','centers:tau_wall','centers:tx','centers:ty','centers:tz','centers:sx','centers:sy','centers:sz'])

    return None

#=============================================================================
# compute the pressure gradients in the wall normal/tangent directions
# needs gradxPressure, gradyPressure and gradzPressure fields
#=============================================================================
def extractLocalPressureGradients(ts):
    """Computes the pressure gradients in the wall normal/tangent direction.""" 
    ts2 = Internal.copyRef(ts)
    _extractLocalPressureGradients(ts2)
    return ts2

def _extractLocalPressureGradients(ts):
    """Computes the pressure gradients in the wall normal/tangent direction.""" 
    if Internal.getNodeFromName(ts, Internal.__FlowSolutionCenters__) is None:
        ts = C.node2Center(ts, Internal.__FlowSolutionNodes__)

    # normal vector
    G._getNormalMap(ts)
    C._normalize(ts, ['centers:sx','centers:sy','centers:sz'])

    # Calc tangent vector
    C._initVars(ts, '{centers:tx}={centers:VelocityX}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sx}')
    C._initVars(ts, '{centers:ty}={centers:VelocityY}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sy}')
    C._initVars(ts, '{centers:tz}={centers:VelocityZ}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sz}')
    C._normalize(ts, ['centers:tx','centers:ty','centers:tz'])

    C._initVars(ts, '{centers:gradnP}={centers:gradxPressure}*{centers:sx}+{centers:gradyPressure}*{centers:sy}+{centers:gradzPressure}*{centers:sz}')
    C._initVars(ts, '{centers:gradtP}={centers:gradxPressure}*{centers:tx}+{centers:gradyPressure}*{centers:ty}+{centers:gradzPressure}*{centers:tz}')
    C._rmVars(ts,['centers:tx','centers:ty','centers:tz','centers:sx','centers:sy','centers:sz'])

    return None

#=============================================================================
# compute yplus_i
#=============================================================================
def extractYplusIP(tc):
    """Computes yplus values at image points and store them in the tc."""
    tp = Internal.copyRef(tc)
    _extractYplusIP(tp)
    return tp

def _extractYplusIP(tc):
    """Computes yplus values at image points and store them in the tc."""
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if (nameSubRegion[:4] == 'IBCD' or nameSubRegion[:4] == '2_IB'):
                yplus = Internal.getNodeFromName(zsr, 'yplus')[1]
                yplusI = numpy.zeros_like(yplus)

                XPW = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                YPW = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                ZPW = Internal.getNodeFromName(zsr, 'CoordinateZ_PW')[1]

                XPC = Internal.getNodeFromName(zsr, 'CoordinateX_PC')[1]
                YPC = Internal.getNodeFromName(zsr, 'CoordinateY_PC')[1]
                ZPC = Internal.getNodeFromName(zsr, 'CoordinateZ_PC')[1]

                XPI = Internal.getNodeFromName(zsr, 'CoordinateX_PI')[1]
                YPI = Internal.getNodeFromName(zsr, 'CoordinateY_PI')[1]
                ZPI = Internal.getNodeFromName(zsr, 'CoordinateZ_PI')[1]

                distCW = numpy.sqrt( (XPW-XPC)*(XPW-XPC) + (YPW-YPC)*(YPW-YPC) + (ZPW-ZPC)*(ZPW-ZPC))
                distIW = numpy.sqrt( (XPW-XPI)*(XPW-XPI) + (YPW-YPI)*(YPW-YPI) + (ZPW-ZPI)*(ZPW-ZPI))

                yplusI = yplus/distCW*distIW

                zsr[2].append(['yplusIP', yplusI, [], 'DataArray_t'])
    return None

#=============================================================================
# Extrapolate the wall pressure (1st order) using pressure gradient info
#=============================================================================
def extractPressureHO(tc, extractDensity=False):
    """Extrapolates the wall pressure (1st order) at the immersed boundaries and stores the solution in the tc."""
    tp = Internal.copyRef(tc)
    _extractPressureHO(tp, extractDensity=extractDensity)
    return tp

def _extractPressureHO(tc, extractDensity=False):
    """Extrapolates the wall pressure (1st order) at the immersed boundaries and stores the solution in the tc."""
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
                if extractDensity: Density  = Internal.getNodeFromName(zsr, 'Density')[1]

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
                    
                    if extractDensity: Density[i]  = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta)
                    Pressure[i] = Pressure[i] - nGradP*beta

    return None

#=============================================================================
# Extrapolate the wall pressure (2nd order) using pressure gradient info
#=============================================================================    
def extractPressureHO2(tc, extractDensity=False):
    """Extrapolates the wall pressure (2nd order) at the immersed boundaries and stores the solution in the tc."""
    tp = Internal.copyRef(tc)
    _extractPressureHO2(tp, extractDensity=extractDensity)
    return tp

def _extractPressureHO2(tc, extractDensity=False):
    """Extrapolates the wall pressure (2nd order) at the immersed boundaries and stores the solution in the tc."""
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
                if extractDensity: Density  = Internal.getNodeFromName(zsr, 'Density')[1]

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

                    if extractDensity: Density[i]  = Density[i]/Pressure[i]*(Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2)
                    Pressure[i] = Pressure[i] - nGradP*beta + 0.5*nnGradP*beta**2

    return None

#=============================================================================
# Computes the convective terms required for the full thin boundary layers equations
#=============================================================================
def extractConvectiveTerms(tc):
    """Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc."""
    tp = Internal.copyRef(tc)
    _extractConvectiveTerms(tp)
    return tp

def _extractConvectiveTerms(tc):
    """Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc."""
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

#=============================================================================
# Computes additional variables required for the IBM post-processing.
# IN: ts (tree): geometry tree with solution at the cell centers
# IN: PInf (float): reference pressure
# IN: QInf (float): reference dynamic pressure
# IN: variables (list): list of extra variables to compute
#=============================================================================
def computeExtraVariables(ts, PInf, QInf, 
                          variables=['Cp','Cf','frictionX','frictionY','frictionZ', 'frictionMagnitude','ShearStress']):
    """Computes additional variables required for the IBM post-processing."""
    ts2 = Internal.copyRef(ts)
    _computeExtraVariables(ts2, PInf, QInf,variables=variables)
    return ts2
                          
def _computeExtraVariables(ts, PInf, QInf, 
                           variables=['Cp','Cf','frictionX','frictionY','frictionZ', 'frictionMagnitude','ShearStress']):
    """Computes additional variables required for the IBM post-processing."""
    import Post.ExtraVariables2 as PE

    #-------------------------
    # pressure forces 
    #-------------------------
    if 'Cp' in variables: 
        PE._extractCp(ts, PInf, QInf)

    if any(var in ['gradnP', 'gradtP'] for var in variables):
        _extractLocalPressureGradients(ts)

    #-------------------------
    # viscous forces
    #-------------------------
    if any(var in ['ShearStress', 'frictionX'] for var in variables):
        _extractShearStress(ts)

    if any(var in ['frictionX', 'frictionMagnitude', 'Cf'] for var in variables):
        PE._extractFrictionMagnitude(ts)
    
    if 'Cf' in variables: 
        PE._extractCf(ts, QInf)

    if 'ShearStress' not in variables: 
        FSC = Internal.getNodesFromName(ts,Internal.__FlowSolutionCenters__)
        Internal._rmNodesFromName(FSC, 'ShearStress*')
    return None

#=============================================================================
# Computes the viscous and pressure forces on the immersed boundaries
# IN: tb (tree): geometry tree (IBM bodies) 
# IN: tc (tree): connectivity tree containing IBM information
# IN: tc2 (tree, optional): connectivity tree containing IBM information for the second image point 
# IN: alpha (float): angle of attack (x-y plane)
# IN: beta (float): angle of attack (x-z plane)
# IN: Sref (float): reference surface for the aerodynamical coefficients (CD/CL)
#       if Sref is None, Sref is computed as the integ of the surface
# IN: order (float): order for the extrapolation of the wall pressure when gradP is True
# IN: gradP (boolean): if True, extract the wall pressure using pressure gradient information
#       if order == 1: 1st order extrapolation
#       if order == 2: 2nd order extrapolation
#       if tc2 is not None: use the information at second image points
#       else: use the information at first image points
# IN: famZones (list): if famZones is not empty, only extract some subregion families (['FAM1','FAM2,...])
# OUT: Surface tree with additional solution fields: shear stress, tangential friction vector and modulus and taun, Cp, Cf
# OUT (screen): CD & CL
#=============================================================================
def _addGradxiP__(tc):
    for z in Internal.getZones(tc):
        subRegions = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
        for zsr in subRegions:
            nameSubRegion = zsr[0]
            if nameSubRegion[:4] == "IBCD":
                pressure = Internal.getNodeFromName(zsr, 'Pressure')[1]
                gradxP   = Internal.getNodeFromName(zsr, 'gradxPressure')

                nIBC    = pressure.shape[0]
                Nlength = numpy.zeros((nIBC),numpy.float64)
                if gradxP is None:                        
                    zsr[2].append(['gradxPressure' , Nlength            , [], 'DataArray_t'])
                    zsr[2].append(['gradyPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
                    zsr[2].append(['gradzPressure' , copy.copy(Nlength) , [], 'DataArray_t'])
    return None

def extractPressureFromFront2__(zw, zw2, extractDensity=False):
    zones_zw  = []
    zones_zw2 = []
    for zone in Internal.getZones(zw): zones_zw.append(zone[0])
    for zone in Internal.getZones(zw2): zones_zw2.append(zone[0])
    nbZones = len(zones_zw)

    for i in range(nbZones): # for multi corps
        szw  = Internal.getNodeFromName(zw, zones_zw[i])
        szw2 = Internal.getNodeFromName(zw2, zones_zw2[i])

        Internal.getNodeFromName(szw, 'Pressure')[1] = Internal.getNodeFromName(szw2, 'Pressure')[1]
        if extractDensity: Internal.getNodeFromName(szw, 'Density')[1]  = Internal.getNodeFromName(szw2, 'Density')[1]

        Internal.getNodeFromName(szw, 'gradxPressure')[1] = Internal.getNodeFromName(szw2, 'gradxPressure')[1]
        Internal.getNodeFromName(szw, 'gradyPressure')[1] = Internal.getNodeFromName(szw2, 'gradyPressure')[1]
        Internal.getNodeFromName(szw, 'gradzPressure')[1] = Internal.getNodeFromName(szw2, 'gradzPressure')[1]
    
    return zw

def _loads0(ts, Sref=None, Pref=None, Qref=None, alpha=0., beta=0., dimPb=3, verbose=0, time=-1):
    import Post.ExtraVariables2 as PE

    if not Internal.getZones(ts):
        print('INFO: loads0: no zones in ts. Returning res and res2.')
        res  = PE.integCp(ts)[0]
        res2 = PE.integTaun(ts)
        return [res, res2, [0,0], [0,0]]

    if Sref is None:
        C._initVars(ts, '__ONE__',1.)
        Sref = P.integ(ts, '__ONE__')[0]
        C._rmVars(ts, ['__ONE__', 'centers:vol'])
        
    RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
    PInf     = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
    RoInf    = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
    VxInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
    VyInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
    VzInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
    VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
    q        = 0.5*RoInf*VInf2
    
    alpha  = math.radians(alpha)
    beta   = math.radians(beta)
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)   
    
    if Qref is None: Qref = q
    if Pref is None: Pref = PInf

    variables = ['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress']

    if Internal.getNodeFromName(ts, 'gradxPressure') is not None: variables += ['gradnP', 'gradtP']
    
    _computeExtraVariables(ts, Pref, Qref, variables=variables)
        
    #===================================
    # Compute pressure & friction forces
    #===================================
    Qadim = Qref*Sref

    res  = PE.integCp(ts)[0]
    res2 = PE.integTaun(ts)

    res    = [-i/Sref for i in res]
    res2   = [ i/Qadim for i in res2]
    
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)

    if dimPb == 3:
        cdp = res[0]*calpha*cbeta  + res[1]*salpha*cbeta  - res[2]*sbeta
        clp = res[1]*calpha        - res[0]*salpha
        cdf = res2[0]*calpha*cbeta + res2[1]*salpha*cbeta - res2[2]*sbeta
        clf = res2[1]*calpha       - res2[0]*salpha
    else:
        cdp = res[0]*calpha  + res[1]*salpha
        clp = res[1]*calpha  - res[0]*salpha
        cdf = res2[0]*calpha + res2[1]*salpha
        clf = res2[1]*calpha - res2[0]*salpha

    if verbose and Cmpi.rank == 0:
        print("Normalized pressure drag = %.4e and lift = %.4e"%(cdp, clp))
        print("Vector of pressure loads: (Fx_P,Fy_P,Fz_P) = (%.4e, %.4e, %.4e)"%(res[0],res[1],res[2]))
    
        print("Normalized skin friction drag = %.4e and lift = %.4e"%(cdf, clf))
        print("Vector of skin friction loads: (Fx_f,Fy_f,Fz_f) = (%.4e,%.4e,%.4e)"%(res2[0], res2[1], res2[2]))

        infoTime = ' (time = %.4e)'%time if time >= 0 else ''
        print("******************************************")
        print("Total Drag%s: %.4e"%(infoTime,(cdp+cdf)))
        print("Total Lift%s: %.4e"%(infoTime,(clp+clf)))
        print("******************************************")

    FSC = Internal.getNodesFromName(ts,Internal.__FlowSolutionCenters__)
    Internal._rmNodesFromName(FSC, 'ShearStress*')

    return [res, res2, [clp, cdp], [clf, cdf]]

def loads(tb_in, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., Sref=None, order=1, gradP=False, famZones=[], extractIBMInfo=False):
    """Computes the viscous and pressure forces on the immersed boundaries"""

    if tc_in is not None:
        if isinstance(tc_in, str):
            tc = C.convertFile2PyTree(tc_in)
        else: tc = tc_in

        if tc2_in is not None:
            if isinstance(tc2_in, str):
                tc2 = C.convertFile2PyTree(tc2_in)
            else: tc2 = tc2_in
        else: tc2 = None
    else: 
        tc = None; tc2 = None
        
    if isinstance(tb_in, str): tb = C.convertFile2PyTree(tb_in)
    else: tb = tb_in

    if Sref is None:
        C._initVars(tb, '__ONE__',1.)
        Sref = P.integ(tb, '__ONE__')[0]; print('Info: loads: no data input for Sref. Default value is taken as the geometry area. Sref = %f'%Sref)
        C._rmVars(tb, ['__ONE__', 'centers:vol'])

    #====================================
    # Wall pressure correction
    #====================================
    if tc is not None and gradP:
        # add gradP fields in tc if necessary
        _addGradxiP__(tc)

        if tc2 is None:
            print('Info: loads: pressure gradients come from tc')
            if order < 2: tc = extractPressureHO(tc)
            else: tc = extractPressureHO2(tc)
        else:
            print('Info: loads: pressure gradients come from tc2')
            if order < 2: tc2 = extractPressureHO(tc2)
            else: tc2 = extractPressureHO2(tc2)

    #====================================
    # Extraction des grandeurs a la paroi
    #====================================
    if tc is None: 
        zw = Internal.getZones(tb)
        zw = T.join(zw)
    else:
        _extractYplusIP(tc)
        zw = extractIBMWallFields(tc, tb=tb, coordRef='wall', famZones=famZones, extractIBMInfo=extractIBMInfo)
    
    #====================================
    # Extract pressure info from tc2 to tc
    #====================================
    if tc2 is not None:
        print('Info: loads: pressure info come from tc2')
        zw2 = extractIBMWallFields(tc2, tb=tb, famZones=famZones)
        zw = extractPressureFromFront2__(zw, zw2)

    dimPb = Internal.getValue(Internal.getNodeFromName(tb, 'EquationDimension'))

    if dimPb == 2: T._addkplane(zw)

    # save basenames for the final tree
    baseNames = []
    for b in Internal.getBases(zw):
        baseNames.append(b[0])

    zw = C.convertArray2Tetra(zw)
    zw = T.reorderAll(zw, 1)

    CD, CL = [], []

    if famZones:
        ts = C.newPyTree(['SKIN'])
        ts[2][1][2] = zw
        RefState = Internal.getNodeFromType(tb,'ReferenceState_t')
        ts[2][1][2].append(RefState)
        tp = C.node2Center(ts,'FlowSolution')
        C._rmVars(ts, 'FlowSolution')
        [res, res2, [clp, cdp], [clf, cdf]] = _loads0(ts, Sref=Sref, Pref=None, Qref=None, alpha=alpha, beta=beta, dimPb=dimPb, verbose=True)

        CD.append(clp+cdp)
        CL.append(clf+cdf)

        if dimPb == 2: # reextrait en 2D
            ts = P.isoSurfMC(ts, "CoordinateZ", 0.)
            nodes = Internal.getNodesFromName(ts, 'CoordinateX')
            xmin = numpy.min(nodes[0][1])
            xmax = numpy.max(nodes[0][1])
            dxi = 1./(xmax-xmin)
            C._initVars(ts, 'xsc=({CoordinateX}-%g)*%g'%(xmin, dxi))

    else: # get loads info per base
        ts = C.newPyTree()
        for c, b in enumerate(Internal.getBases(zw)):
            tp = C.newPyTree(['SKIN']);
            tp[2][1][2] = zw[2][c+1][2]
            RefState = Internal.getNodeFromType(tb,'ReferenceState_t')
            tp[2][1][2].append(RefState)
            tp = C.node2Center(tp,'FlowSolution')
            C._rmVars(tp, 'FlowSolution')
            print("Info: loads: get CD/CL on base %s"%(baseNames[c]))
            [res, res2, [clp, cdp], [clf, cdf]] = _loads0(tp, Sref=Sref, Pref=None, Qref=None, alpha=alpha, beta=beta, dimPb=dimPb, verbose=True)

            CD.append(cdp+cdf)
            CL.append(clp+clf)

            if dimPb == 2: # reextrait en 2D
                tp = P.isoSurfMC(tp, "CoordinateZ", 0.)
                nodes = Internal.getNodesFromName(tp, 'CoordinateX')
                xmin = numpy.min(nodes[0][1])
                xmax = numpy.max(nodes[0][1])
                dxi = 1./(xmax-xmin)
                C._initVars(tp, 'xsc=({CoordinateX}-%g)*%g'%(xmin, dxi))

            base = Internal.newCGNSBase(baseNames[c], parent=ts)
            base[2] += Internal.getZones(tp)

    if isinstance(wall_out, str): C.convertPyTree2File(ts, wall_out)
    return ts, CL, CD

##################
# WORK IN PROGRESS
##################

#==============================================================================
# Return ts, massflow 
# A REFAIRE !!!  [AJ] KEEP FOR NOW
def extractMassFlowThroughSurface(tb, t, famZones=[]):
    print("WARNING: fonction a reprendre !!!")
    ts = Internal.copyRef(tb)
    if famZones !=[]:
        for zone in Internal.getZones(ts):
            familyNode = Internal.getNodeFromName(zone, "FamilyName")
            if familyNode is None or Internal.getValue(familyNode) not in famZones:
                Internal._rmNodesByName(ts, zone[0])

    C._initVars(t,'{centers:cellN}=minimum({centers:cellN},1.)')
    P._extractMesh(t, ts, mode='robust')
    C._rmVars(ts, 'cellN')
    C._initVars(ts, '{MomentumX}={Density}*{VelocityX}')
    C._initVars(ts, '{MomentumY}={Density}*{VelocityY}')
    C._initVars(ts, '{MomentumZ}={Density}*{VelocityZ}')
    massflow = P.integNormProduct(ts, vector=['MomentumX', 'MomentumY', 'MomentumZ'])
    return massflow, ts

#==========================================================================================
# IN:  ts            : skin (TRI zones) distributed already (partial tree here)
# IN:  tc            : connectivity tree
# OUT: tl            : NODE-type zones of IBM points to be projected locally on ts
# OUT: graphWPOST    : graph of coms between tc and tl
# OUT: interDictWPOST: dictionary of intersection domains
#==========================================================================================
def _prepareSkinReconstruction(ts, tc,famZones=[]):
    alphah=2.2
    tBBs=Cmpi.createBBoxTree(ts)
    procDictFUS = Cmpi.getProcDict(tBBs)

    basename=Internal.getName(Internal.getBases(ts)[0])
    tl = C.newPyTree([basename])
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            fam_tmp = Internal.getNodeFromType(IBCD,'FamilyName_t')
            if fam_tmp is not None:fam   = Internal.getValue(fam_tmp)
            if famZones and fam not in famZones:continue
            zname = Internal.getValue(IBCD)
            XW = Internal.getNodeFromName(IBCD,'CoordinateX_PW')[1]
            YW = Internal.getNodeFromName(IBCD,'CoordinateY_PW')[1]
            ZW = Internal.getNodeFromName(IBCD,'CoordinateZ_PW')[1]
            zsize = numpy.empty((1,3), numpy.int32, order='F')
            zsize[0,0] = XW.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
            z = Internal.newZone(name='IBW_Wall_%s_%s'%(zc[0],zname),zsize=zsize,
                                 ztype='Unstructured')
            gc = Internal.newGridCoordinates(parent=z)
            coordx = ['CoordinateX',XW,[],'DataArray_t']
            coordy = ['CoordinateY',YW,[],'DataArray_t']
            coordz = ['CoordinateZ',ZW,[],'DataArray_t']
            gc[2] = [coordx,coordy,coordz]
            n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
            Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
            Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
            FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                           gridLocation='Vertex', parent=z)
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
                utauPresent = 1
                FSN[2].append([XOD.__UTAU__,utauNP[0], [],'DataArray_t'])
            if yplusNP != []:
                yplusPresent = 1
                FSN[2].append([XOD.__YPLUS__,yplusNP[0], [],'DataArray_t'])

            if vxNP != []:
                vxPresent = 1
                FSN[2].append([XOD.__VELOCITYX__,vxNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYY__,vyNP[0], [],'DataArray_t'])
                FSN[2].append([XOD.__VELOCITYZ__,vzNP[0], [],'DataArray_t'])

            Cmpi._setProc(z,Cmpi.rank)
            tl[2][1][2].append(z)

    tlBB=Cmpi.createBBoxTree(tl)
    for zbb in Internal.getZones(tlBB):
        if Internal.getNodeByName(zbb,'FlowSolution') is not None:
            G._BB(zbb)
 
    hmin = -1e10
    for zbb in Internal.getZones(tlBB):    
        zc = Internal.getNodeFromName2(tl,zbb[0])
        hloc = -1e-3
        if zc is not None:
            GCnode = Internal.getNodeFromType(zc,"GridCoordinates_t")
            XN = Internal.getNodeFromName(GCnode,'CoordinateX')
            if XN is not None:
                if XN[1].shape[0]>1:
                    hx = C.getValue(zc,'CoordinateX',1)-C.getValue(zc,'CoordinateX',0)
                    hy = C.getValue(zc,'CoordinateY',1)-C.getValue(zc,'CoordinateY',0)
                    hloc = max(abs(hx),abs(hy))
                    hmin = max(hloc,hmin)
                else:
                    hloc = 1e-3
                xmin = C.getMinValue(zbb,'CoordinateX'); xmax = C.getMaxValue(zbb,'CoordinateX')
                ymin = C.getMinValue(zbb,'CoordinateY'); ymax = C.getMaxValue(zbb,'CoordinateY')
                zmin = C.getMinValue(zbb,'CoordinateZ'); zmax = C.getMaxValue(zbb,'CoordinateZ')
                for k in range(2):
                    for j in range(2):
                        for i in range(2):
                            ind = i+j*2+k*4
                            if i==0: xl = xmin-alphah*hloc
                            else: xl = xmax+alphah*hloc
                            if j==0: yl = ymin-alphah*hloc
                            else: yl = ymax+alphah*hloc
                            if k==0: zl = zmin-alphah*hloc
                            else: zl = zmax+alphah*hloc

                            C.setValue(zbb,'CoordinateX',ind,xl)
                            C.setValue(zbb,'CoordinateY',ind,yl)
                            C.setValue(zbb,'CoordinateZ',ind,zl)

    procDictWPOST = Cmpi.getProcDict(tlBB)
    interDictWPOST = X.getIntersectingDomains(tlBB, tBBs)
    graphWPOST = Cmpi.computeGraph(tlBB, type='bbox3',intersectionsDict=interDictWPOST,
                                   procDict=procDictWPOST, procDict2=procDictFUS, t2=tBBs)
 
    C._initVars(ts,"Pressure",0.)
    C._initVars(ts,"Density",0.)
    C._initVars(ts,XOD.__UTAU__,0.)
    C._initVars(ts,XOD.__YPLUS__,0.)
    C._initVars(ts,XOD.__VELOCITYX__,0.)
    C._initVars(ts,XOD.__VELOCITYY__,0.)
    C._initVars(ts,XOD.__VELOCITYZ__,0.)

    RefStateNode = Internal.getNodeFromName(ts,'ReferenceState')
    if RefStateNode is not None:
        tl[2][1][2].append(RefStateNode)
    FES =  Internal.getNodeFromName(ts,'FlowEquationSet')
    if FES is not None:
        tl[2][1][2].append(FES)
    return tl, graphWPOST


def _computeSkinVariables(ts, tc, tl, graphWPOST,famZones=[]):
    for zc in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromType(zc,"ZoneSubRegion_t")
        allIBCD = Internal.getNodesFromName(allIBCD,"IBCD_*")
        for IBCD in allIBCD:
            fam_tmp = Internal.getNodeFromType(IBCD,'FamilyName_t')
            if fam_tmp is not None: fam   = Internal.getValue(fam_tmp)
            if famZones and fam not in famZones: continue
            zname = Internal.getValue(IBCD)
            znamepostw = 'IBW_Wall_%s_%s'%(zc[0],zname)
            zpostw = Internal.getNodeFromName(tl,znamepostw)
            FSP = Internal.getNodeFromType(zpostw,'FlowSolution_t')
            
            PW     = Internal.getNodeFromName1(IBCD,XOD.__PRESSURE__)
            RHOW   = Internal.getNodeFromName1(IBCD,XOD.__DENSITY__)
            UTAUW  = Internal.getNodeFromName1(IBCD,XOD.__UTAU__)
            YPLUSW = Internal.getNodeFromName1(IBCD, XOD.__YPLUS__)
            VXW    = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYX__)
            VYW    = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYY__)
            VZW    = Internal.getNodeFromName1(IBCD, XOD.__VELOCITYZ__)
            
            
            PW2     = Internal.getNodeFromName1(FSP,XOD.__PRESSURE__)
            RHOW2   = Internal.getNodeFromName1(FSP,XOD.__DENSITY__)
            UTAUW2  = Internal.getNodeFromName1(FSP,XOD.__UTAU__)
            YPLUSW2 = Internal.getNodeFromName1(FSP, XOD.__YPLUS__)
            VXW2    = Internal.getNodeFromName1(FSP, XOD.__VELOCITYX__)
            VYW2    = Internal.getNodeFromName1(FSP, XOD.__VELOCITYY__)
            VZW2    = Internal.getNodeFromName1(FSP, XOD.__VELOCITYZ__)
            
            PW2[1]  =PW[1]
            RHOW2[1]=RHOW[1]
            VXW2[1] =VXW[1]
            VYW2[1] =VYW[1]
            VZW2[1] =VZW[1]

            if UTAUW2 is not None:
                UTAUW2[1] =UTAUW[1]
                YPLUSW2[1]=YPLUSW[1]
            
    tdl = Cmpi.addXZones(tl, graphWPOST)
    tdl = Cmpi.convert2PartialTree(tdl)
    for nobs in range(len(ts[2])):
        if Internal.getType(ts[2][nobs])=='CGNSBase_t':
            for nozs in range(len(ts[2][nobs][2])):
                zs = ts[2][nobs][2][nozs]
                if Internal.getType(zs)=='Zone_t':
                    cloud = []
                    for zl in Internal.getZones(tdl):
                        if zl != [] and zl is not None:
                            zl = C.convertArray2Node(zl)
                            cloud.append(zl)
    
                    if cloud != []:
                        cloud = T.join(cloud)
                        ts[2][nobs][2][nozs] = P.projectCloudSolution(cloud, zs, dim=3)
                        
    return None

