"""Toolbox for IBM preprocessing"""
import numpy
from . import PyTree as X
from . import IBM as X_IBM
from . import OversetData as XOD
from . import Connector
from . import connector

try: range = xrange
except: pass

try:
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Generator.IBM as G_IBM
    import Generator.IBMmodelHeight as G_IBM_Height
    import Transform.PyTree as T
    import Converter.Internal as Internal
    import Dist2Walls.PyTree as DTW
    import Post.PyTree as P
    import Post.IBM as P_IBM
    import Converter
    import Generator
    import Transform
    import Converter.GhostCells as CGC
    import KCore
    import numpy
    import math
except:
    raise ImportError("Connector.ToolboxIBM requires Converter, Generator, Transform, Dist2Walls and Post modules.")


TypesOfIBC = XOD.TypesOfIBC

## IMPORTANT NOTE !!
## This PYTHON FILE will become decrepit after Jan. 1 2023
#====================================================================================
def _blankClosestTargetCells(t, cellNName='cellN', depth=3):
    X_IBM._blankClosestTargetCells(t, cellNName=cellNName, depth=depth)
    return None


def _removeBlankedGrids(t, loc='centers'):
    X_IBM._removeBlankedGrids(t, loc=loc)
    return None


def blankByIBCBodies(t, tb, loc, dim, cellNName='cellN'):
    t=X_IBM.blankByIBCBodies(t, tb, loc, dim, cellNName=cellNName)
    return t


def _addBCOverlaps(t, bbox):
    X_IBM._addBCOverlaps(t, bbox)
    return None


def _addExternalBCs(t, bbox, DEPTH=2, externalBCType='BCFarfield', dimPb=3):
    X_IBM._addExternalBCs(t, bbox, DEPTH=DEPTH, externalBCType=externalBCType, dimPb=dimPb)
    return None


def _modifPhysicalBCs__(zp, depth=2, dimPb=3):
    X_IBM._modifPhysicalBCs__(zp, depth=depth, dimPb=dimPb)
    return None


def getIBMFront(tc, frontvar, dim, frontType, isoFront=False, isFront2=False, SHIFTB=0.):
    front=X_IBM.getIBMFront(tc, frontvar, dim, frontType, isoFront=isoFront, isFront2=isFront2, SHIFTB=SHIFTB)
    return front


def getIBMFrontType1(tc, frontvar, dim):
    front=X_IBM.getIBMFrontType1(tc, frontvar, dim)
    return front


def getIBMFrontType0(tc, frontvar, dim, isFront2=False, frontType=0, SHIFTB=0.):
    front=X_IBM.getIBMFrontType0(tc, frontvar, dim, isFront2=isFront2, frontType=frontType, SHIFTB=SHIFTB)
    return front


def getIBMFrontType0_old(tc, frontvar, dim, isFront2=False, frontType=0, SHIFTB=0.):
    front=X_IBM.getIBMFrontType0_old(tc, frontvar, dim, isFront2=isFront2, frontType=frontType, SHIFTB=SHIFTB)
    return front


def _pushBackImageFront2(t, tc, tbb, interpDataType=1):
    X_IBM._pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)
    return None


def _smoothImageFront(t, tc, dimPb=2):
    X_IBM._smoothImageFront(t, tc, dimPb=dimPb)
    return None


def _smoothImageFrontBackward(t, tc, dimPb=2):
    X_IBM._smoothImageFrontBackward(t, tc, dimPb=dimPb)
    return None


def gatherFront(front):
    front=X_IBM.gatherFront(front)
    return front


def doInterp(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1., isLBM=False):    
    tc=X_IBM.doInterp(t, tc, tbb, tb=tb, typeI=typeI, dim=dim, dictOfADT=dictOfADT, front=front, frontType=frontType, depth=depth, IBCType=IBCType, interpDataType=interpDataType, Reynolds=Reynolds, yplus=yplus, Lref=Lref, isLBM=isLBM)
    return tc


def doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1.):
    tc=X_IBM.doInterp2(t, tc, tbb, tb=tb, typeI=typeI, dim=dim, dictOfADT=dictOfADT, front=front, frontType=frontType, depth=depth, IBCType=IBCType, interpDataType=interpDataType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)
    return tc


def doInterp3(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100., Lref=1.):
    tc=X_IBM.doInterp3(t, tc, tbb, tb=tb, typeI=typeI, dim=dim, dictOfADT=dictOfADT, frontType=frontType, depth=depth, IBCType=IBCType, interpDataType=interpDataType, Reynolds=Reynolds, yplus=yplus, Lref=Lref)
    return tc


def _extractIBMInfo_param(t,tc):
    X_IBM._extractIBMInfo_param(t,tc)
    return None


def extractIBMInfo(tc):
    t=X_IBM.extractIBMInfo(tc)
    return t


def extractIBMInfo2(tc):
    t=X_IBM.extractIBMInfo2(tc)
    return t


def getAllIBMPoints(t, loc='nodes', hi=0., he=0., tb=None, tfront=None, tfront2=None, frontType=0, cellNName='cellN', IBCType=1, depth=2, Reynolds=6.e6, yplus=100., Lref=1., hmod=0.1, isLBM=False):
    dictOfCorrectedPtsByIBCType, dictOfWallPtsByIBCType, dictOfInterpPtsByIBCType=X_IBM.getAllIBMPoints(t, loc=loc, hi=hi, he=he, tb=tb, tfront=tfront, tfront2=tfront2, frontType=frontType, cellNName=cellNName, IBCType=IBCType, depth=depth, Reynolds=Reynolds, yplus=yplus, Lref=Lref, hmod=hmod, isLBM=isLBM)
    return dictOfCorrectedPtsByIBCType, dictOfWallPtsByIBCType, dictOfInterpPtsByIBCType


def prepareIBMData(t, tbody, DEPTH=2, loc='centers', frontType=1, interpDataType=0, smoothing=False, yplus=100., Lref=1., wallAdapt=None, blankingF42=False, isLBM=False,LBMQ=False,isPrintDebug=False):
    t,tc=X_IBM.prepareIBMData(t, tbody, DEPTH=DEPTH, loc=loc, frontType=frontType, interpDataType=interpDataType, smoothing=smoothing, yplus=yplus, Lref=Lref, wallAdapt=wallAdapt, blankingF42=blankingF42, isLBM=isLBM,LBMQ=LBMQ,isPrintDebug=isPrintDebug)
    return t, tc


def prepareIBMData2(t, tbody, DEPTH=2, loc='centers', frontType=1, inv=False, interpDataType=1):
    t,tc=X_IBM.prepareIBMData2(t, tbody, DEPTH=DEPTH, loc=loc, frontType=frontType, inv=inv, interpDataType=interpDataType)
    return t, tc


def createWallAdapt(tc):
    t=X_IBM.createWallAdapt(tc)
    return t


def createIBMWZones(tc,variables=[]):
    tw=X_IBM.createIBMWZones(tc,variables=variables)
    return tw


def _computeKcurvParameter(tc, tb):
    X_IBM._computeKcurvParameter(tc, tb)
    return None


def _signDistance(t):
    X_IBM._signDistance(t)
    return None


def generateCartMesh__(o, parento=None, dimPb=3, vmin=11, DEPTH=2, sizeMax=4000000, check=True,
                       symmetry=0, externalBCType='BCFarfield', bbox=None):
    t=G_IBM.generateCartMesh__(o, parento=parento, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH, sizeMax=sizeMax, check=check,
                       symmetry=symmetry, externalBCType=externalBCType, bbox=bbox)
    return t


def adaptIBMMesh(t, tb, vmin, sensor, factor=1.2, DEPTH=2, sizeMax=4000000,
                 variables=None, refineFinestLevel=False, refineNearBodies=False,
                 check=True, symmetry=0, externalBCType='BCFarfield', fileo='octree.cgns'):
    t2=G_IBM.adaptIBMMesh(t, tb, vmin, sensor, factor=factor, DEPTH=DEPTH, sizeMax=sizeMax,
                 variables=variables, refineFinestLevel=refineFinestLevel, refineNearBodies=refineNearBodies,
                 check=check, symmetry=symmetry, externalBCType=externalBCType, fileo=fileo)
    return t2


def generateIBMMesh(tb, vmin=15, snears=None, dfar=10., dfarList=[], DEPTH=2, tbox=None,
                    snearsf=None, check=True, sizeMax=4000000,
                    symmetry=0, externalBCType='BCFarfield', to=None,
                    fileo=None, expand=2, dfarDir=0):
    res = G_IBM.generateIBMMesh(tb, vmin=vmin, snears=snears, dfar=dfar, dfarList=dfarList, DEPTH=DEPTH, tbox=tbox,
                    snearsf=snearsf, check=check, sizeMax=sizeMax,symmetry=symmetry, externalBCType=externalBCType, to=to,
                    fileo=fileo, expand=expand, dfarDir=dfarDir)
    return res


def buildOctree(tb, snears=None, snearFactor=1., dfar=10., dfarList=[], to=None, tbox=None, snearsf=None,
                dimPb=3, vmin=15, balancing=2, symmetry=0, fileout=None, rank=0, expand=2, dfarDir=0):
    o = G_IBM.buildOctree(tb, snears=snears, snearFactor=snearFactor, dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                dimPb=dimPb, vmin=vmin, balancing=balancing, symmetry=symmetry, fileout=fileout, rank=rank, expand=expand, dfarDir=dfarDir)
    return o


def addRefinementZones(o, tb, tbox, snearsf, vmin, dim):
    tout = G_IBM.addRefinementZones(o, tb, tbox, snearsf, vmin, dim)
    return tout


def octree2StructLoc__(o, parento=None, vmin=21, ext=0, optimized=0, sizeMax=4e6):
    zones = G_IBM.octree2StructLoc__(o, parento=parento, vmin=vmin, ext=ext, optimized=optimized, sizeMax=sizeMax)
    return zones


def mergeByParent__(zones, parent, sizeMax):
    res = G_IBM.mergeByParent__(zones, parent, sizeMax)
    return res


def buildParentOctrees__(o, tb, snears=None, snearFactor=4., dfar=10., dfarList=[], to=None, tbox=None, snearsf=None,
                         dimPb=3, vmin=15, symmetry=0, fileout=None, rank=0, dfarDir=0):
    OCTREEPARENTS = G_IBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=snearFactor, dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=fileout, rank=rank, dfarDir=dfarDir)
    return OCTREEPARENTS


def getMinimumCartesianSpacing(t):
    dxmin=G_IBM_Height.getMinimumCartesianSpacing(t)
    return dxmin


def computeYplus(Re, Cf_law='ANSYS', height=0.1, L=1.):
    val=G_IBM_Height.computeYplus(Re, Cf_law=Cf_law, height=height, L=L)
    return val


def computeModelisationHeight(Re, Cf_law='ANSYS', yplus=100., L=1.):
    val=G_IBM_Height.computeModelisationHeight(Re, Cf_law=Cf_law,yplus=yplus,L=L)
    return val


def computeBestModelisationHeight(Re, h, Cf_law='ANSYS', L=1., q=1.2):
    val1,val2=G_IBM_Height.computeBestModelisationHeight(Re, h, Cf_law=Cf_law, L=L, q=q)
    return val1,val2


def extractIBMWallFields(tc, tb=None, coordRef='wall', famZones=[], front=1):
    td=P_IBM.extractIBMWallFields(tc, tb=tb, coordRef=coordRef, famZones=famZones, front=front)
    return td


