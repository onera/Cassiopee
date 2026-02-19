"""Mesh generation for IBM"""
from . import Generator
from . import PyTree as G

import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.IBM as X_IBM
import Geom.PyTree as D
import Post.PyTree as P
import Converter
import Converter.GhostCells as CGC
import Connector.PyTree as X
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import numpy
import math

EPSCART = 1.e-6

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## LEGACY FUNCTIONS (SEQUENTIAL BEHAVIOR)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def generateCartMesh__(o, parento=None, dimPb=3, vmin=11, DEPTH=2, sizeMax=4000000, check=False,
                       externalBCType='BCFarfield', bbox=None):

    # Estimation du nb de pts engendres
    vminv0 = vmin+2*DEPTH
    vminv = vminv0*vminv0
    if dimPb == 3: vminv=vminv*vminv0
    else: vminv = vminv*2
    nzones0 = Internal.getZoneDim(o)[2]
    npts = nzones0*vminv
    sizeMax = int(sizeMax)
    # DEPTH > 2: ghost cells added for better implicit phase process
    if DEPTH > 2: optimized = 0
    else: optimized = 1
    if DEPTH == 0: ext=0
    else: ext = DEPTH+1

    res = octree2StructLoc__(o, vmin=vmin, ext=ext, optimized=optimized, sizeMax=sizeMax,
                             parento=parento)
    t = C.newPyTree(['CARTESIAN', res])

    dz = 0.01
    if dimPb == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    if bbox is None: bbox = G.bbox(o)
    del o
    _addExternalBCs(t, bbox, DEPTH, externalBCType, dimPb)

    nptsTot = 0
    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        nptsTot += niz*njz*nkz
    print('Expected number of points is %d.'%nptsTot)
    return t

def adaptIBMMesh(t, tb, vmin, sensor, factor=1.2, DEPTH=2, sizeMax=4000000,
                 variables=None, refineFinestLevel=False, refineNearBodies=False,
                 check=False, externalBCType='BCFarfield', fileo='octree.cgns',
                 isAMR=False,valMin=0,valMax=1):
    if fileo is None: raise ValueError("adaptIBMMesh: Octree mesh must be specified by a file.")
    try: to = C.convertFile2PyTree(fileo)
    except: raise ValueError("adaptIBMMesh: %s file not found."%fileo)

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('adaptIBMMesh: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    refstate = Internal.getNodeFromName(tb, 'ReferenceState')

    if refineNearBodies: constraintSurfaces = []
    else: constraintSurfaces = Internal.getZones(tb)
    if refineFinestLevel: refineLevelF = 1
    else: refineLevelF = 0

    o = Internal.getZones(to)[0]
    dims = Internal.getZoneDim(o)
    npts = dims[1]
    C._initVars(t,"{%s}={%s}*({centers:cellN}>0.)*({centers:cellN}<2.)"%(sensor,sensor))
    C._initVars(to, "centers:indicator", 1.)
    to = P.computeIndicatorValue(to, t, sensor) #projects values from t onto octree for var (but only the absolute value)
    res = P.computeIndicatorField(to, sensor, nbTargetPts=factor*npts, \
                                  bodies=constraintSurfaces, \
                                  refineFinestLevel=refineLevelF, \
                                  coarsenCoarsestLevel=1,
                                  isAMR=isAMR,valMin=valMin,valMax=valMax)
    # nettoyage : on n interpole pas tout
    if variables is not None:
        for z in Internal.getZones(t):
            varsc = C.getVarNames(z, excludeXYZ=True,loc='centers')[0]
            for v in varsc:
                if v not in variables: C._rmVars(z, v)

    # adaptation
    if len(res) == 3: to = res[0]
    o = Internal.getZones(to)[0]
    o = G.adaptOctree(o, balancing=2)
    if Cmpi.size==1: C.convertPyTree2File(o, fileo)

    t2 = generateCartMesh__(o, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH,
                            sizeMax=sizeMax, check=check, externalBCType=externalBCType)

    #C._initVars(t2,"centers:Density=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))
    #C._initVars(t2,"centers:VelocityX=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityX'))))
    #C._initVars(t2,"centers:VelocityY=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityY'))))
    #C._initVars(t2,"centers:VelocityZ=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityZ'))))
    #C._initVars(t2,"centers:Temperature=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Temperature'))))
    #C._initVars(t2,"centers:TurbulentSANuTilde=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'TurbulentSANuTildeDensity'))/Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))

    # interpolate the solution on the new mesh
    P._extractMesh(t, t2, 3, mode='accurate')
    return t2

def generateIBMMesh_legacy(tb, vmin=15, snears=0.01, dfar=10., dfarList=[], DEPTH=2, tbox=None,
                           snearsf=None, check=False, sizeMax=4000000,
                           externalBCType='BCFarfield', to=None,
                           fileo=None, expand=2, dfarDir=0, octreeMode=0):
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('generateIBMMesh: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('generateIBMMesh: GoverningEquations is missing in input body tree.')
    # check Euler non consistant avec Musker

    if Internal.getValue(model) == 'Euler':
        for z in Internal.getZones(tb):
            ibctype = Internal.getNodeFromName2(z, 'ibctype')
            if ibctype is not None:
                ibctype = Internal.getValue(ibctype)
                if ibctype == 'Musker' or ibctype == 'Log':
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))
    if to is None:
        o = buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                        dimPb=dimPb, vmin=vmin, fileout=fileo, rank=0,
                        expand=expand, dfarDir=dfarDir, octreeMode=octreeMode)
    else:
        o = Internal.getZones(to)[0]
    if check: C.convertPyTree2File(o, "octree.cgns")

    # retourne les 4 quarts (en 2D) de l'octree parent 2 niveaux plus haut
    # et les 8 octants en 3D sous forme de listes de zones non structurees
    parento = buildParentOctrees__(o, tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=4., dfars=dfar, dfarDir=dfarDir, to=to, tbox=tbox, snearsf=snearsf, octreeMode=octreeMode)
    res = generateCartMesh__(o, parento=parento, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH, sizeMax=sizeMax,
                             check=check, externalBCType=externalBCType)
    return res

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## MACRO FUNCTIONS FOR IBMs
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#==============================================================================
# IN: bbox: bbox des frontieres exterieures
#
# _modifPhysicalBCs__: Reduction de la taille des fenetres des BC physiques
#  pour qu'elles soient traitees comme des ghost cells

# on the other side of the sym plane using DEPTH ghost cells
# BCOverlaps must not be defined on the other side of the symmetry plane
#==============================================================================
def _modifPhysicalBCs__(zp, depth=2, dimPb=3):
    dimZone = Internal.getZoneDim(zp)

    # Physical BCs
    bclist = Internal.getNodesFromType2(zp, 'BC_t')
    for bc in bclist:
        prange = Internal.getNodesFromName1(bc, 'PointRange')
        if prange != []:
            direction = CGC.getDirection__(dimPb, prange)
            # change PointRange for extended mesh
            pr = numpy.copy(prange[0][1])
            ijk = int(direction/2)
            minmax = direction%2
            for dirl in range(dimZone[4]):
                if dirl != ijk:
                    if dimPb == 2 and dirl == 2: pass
                    else:
                        if dirl == 0: N = dimZone[1]
                        elif dirl == 1: N = dimZone[2]
                        else: N = dimZone[3]
                        pr[dirl][0] += depth
                        pr[dirl][1] -= depth
            prange[0][1] = pr
    return None

def _addExternalBCs(t, bbox, DEPTH=2, externalBCType='BCFarfield', dimPb=3):
    """Add external BCs of given type to BCs out of or on bbox."""
    dirs = [0,1,2,3,4,5]
    rangeDir=['imin','jmin','kmin','imax','jmax','kmax']
    if dimPb == 2: dirs = [0,1,3,4]
    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        indM = niz-1+(njz-1)*niz+(nkz-1)*niz*njz
        x1 = C.getValue(zp,'CoordinateX',0)
        y1 = C.getValue(zp,'CoordinateY',0)
        z1 = C.getValue(zp,'CoordinateZ',0)
        x2 = C.getValue(zp,'CoordinateX',indM)
        y2 = C.getValue(zp,'CoordinateY',indM)
        z2 = C.getValue(zp,'CoordinateZ',indM)
        bbz=[x1,y1,z1,x2,y2,z2]
        external = False
        for idir in dirs:
            if abs(bbz[idir]-bbox[idir])< 1.e-6:
                C._addBC2Zone(zp, 'external', externalBCType, rangeDir[idir])
                external = True
        if externalBCType != 'BCOverlap' and externalBCType != 'BCDummy':
            if external: _modifPhysicalBCs__(zp, depth=DEPTH, dimPb=dimPb)
    return None

def _addBCOverlaps(t, bbox):
    """Add BCOverlap boundary condition to BCs entirely inside bbox."""
    xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
    xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]
    for z in Internal.getZones(t):
        # [x1,y1,z1,x2,y2,z2] = G.bbox(z)
        dimZ = Internal.getZoneDim(z)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        indM = niz-1+(njz-1)*niz+(nkz-1)*niz*njz
        x1 = C.getValue(z,'CoordinateX',0)
        y1 = C.getValue(z,'CoordinateY',0)
        z1 = C.getValue(z,'CoordinateZ',0)
        x2 = C.getValue(z,'CoordinateX',indM)
        y2 = C.getValue(z,'CoordinateY',indM)
        z2 = C.getValue(z,'CoordinateZ',indM)
        if x1 > xmin+EPSCART: C._addBC2Zone(z,'overlap1','BCOverlap','imin')
        if x2 < xmax-EPSCART: C._addBC2Zone(z,'overlap2','BCOverlap','imax')
        if y1 > ymin+EPSCART: C._addBC2Zone(z,'overlap3','BCOverlap','jmin')
        if y2 < ymax-EPSCART: C._addBC2Zone(z,'overlap4','BCOverlap','jmax')
        if z1 > zmin+EPSCART: C._addBC2Zone(z,'overlap5','BCOverlap','kmin')
        if z2 < zmax-EPSCART: C._addBC2Zone(z,'overlap6','BCOverlap','kmax')
    return None

def _addBCsForSymmetry(t, bbox=None, dimPb=3, dir_sym=0, X_SYM=0., depth=2):
    if bbox is None: return None
    dirs = [0,1,2,3,4,5]
    rangeDir=['imin','jmin','kmin','imax','jmax','kmax']
    if dimPb == 2: dirs = [0,1,3,4]

    xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
    xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]

    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        indM = niz-1+(njz-1)*niz+(nkz-1)*niz*njz
        x1 = C.getValue(zp,'CoordinateX',0)
        y1 = C.getValue(zp,'CoordinateY',0)
        z1 = C.getValue(zp,'CoordinateZ',0)
        x2 = C.getValue(zp,'CoordinateX',indM)
        y2 = C.getValue(zp,'CoordinateY',indM)
        z2 = C.getValue(zp,'CoordinateZ',indM)
        # blocs interieurs - tt en overlap
        if (x1 > xmin-EPSCART and y1 > ymin-EPSCART and z1 > zmin-EPSCART and
                x2 < xmax+EPSCART and y2 < ymax+EPSCART and z2 < zmax+EPSCART):
            C._fillEmptyBCWith(zp, 'overlap','BCOverlap',dim=dimPb)
        else: # blocs frontieres
            irange = None
            isw = 1; jsw = 1; ksw = 1; iew = niz; jew = njz; kew = nkz
            # ajout des CL physiques sur ghost cells
            t_imin = False; t_imax = False; t_jmin = False; t_jmax = False; t_kmin = False; t_kmax = False
            if dir_sym==1:
                if x1<X_SYM:
                    isw = depth+1
                    t_imin = True
                    C._addBC2Zone(zp,'sym','BCSymmetryPlane','imin')
                if x2>xmax-EPSCART:
                    iew = niz-depth
                    t_imax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','imax')
                if y1<ymin-EPSCART:
                    jsw = depth+1
                    t_jmin = True
                    C._addBC2Zone(zp,'nref','BCFarfield','jmin')
                if y2>ymax+EPSCART:
                    jew = njz-depth
                    t_jmax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','jmax')
                if dimPb==3:
                    if z1 < zmin-EPSCART:
                        ksw = depth+1
                        t_kmin = True
                        C._addBC2Zone(zp,'nref','BCFarfield','kmin')
                    if z2 > zmax+EPSCART:
                        kew = nkz-depth
                        t_kmax = True
                        C._addBC2Zone(zp,'nref','BCFarfield','kmax')

            elif dir_sym==2:
                if y1<X_SYM:
                    jsw = depth+1
                    t_jmin = True
                    C._addBC2Zone(zp,'sym','BCSymmetryPlane','jmin')
                if y2>ymax-EPSCART:
                    jew = njz-depth
                    t_jmax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','jmax')
                if x1<xmin-EPSCART:
                    isw = depth+1
                    t_imin = True
                    C._addBC2Zone(zp,'nref','BCFarfield','imin')
                if x2>xmax+EPSCART:
                    iew = niz-depth
                    t_imax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','imax')
                if dimPb==3:
                    if z1 < zmin-EPSCART:
                        ksw = depth+1
                        t_kmin = True
                        C._addBC2Zone(zp,'nref','BCFarfield','kmin')
                    if z2 > zmax+EPSCART:
                        kew = nkz-depth
                        t_kmax = True
                        C._addBC2Zone(zp,'nref','BCFarfield','kmax')

            elif dir_sym==3:
                if z1<X_SYM:
                    ksw = depth+1
                    t_kmin = True
                    C._addBC2Zone(zp,'sym','BCSymmetryPlane','kmin')
                if z2>zmax-EPSCART:
                    kew = nkz-depth
                    t_kmax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','kmax')
                if x1<xmin-EPSCART:
                    isw = depth+1
                    t_imin = True
                    C._addBC2Zone(zp,'nref','BCFarfield','imin')
                if x2>xmax+EPSCART:
                    iew = niz-depth
                    t_imax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','imax')
                if y1<ymin-EPSCART:
                    jsw = depth+1
                    t_jmin = True
                    C._addBC2Zone(zp,'nref','BCFarfield','jmin')
                if y2>ymax+EPSCART:
                    jew = njz-depth
                    t_jmax = True
                    C._addBC2Zone(zp,'nref','BCFarfield','jmax')

            if not t_jmin:
                C._addBC2Zone(zp,'overlap','BCOverlap',[isw, iew, jsw, jsw, ksw, kew])
            if not t_jmax:
                C._addBC2Zone(zp,'overlap','BCOverlap',[isw, iew, jew, jew, ksw, kew])
            if not t_imin:
                C._addBC2Zone(zp,'overlap','BCOverlap',[isw, isw, jsw, jew, ksw, kew])
            if not t_imax:
                C._addBC2Zone(zp,'overlap','BCOverlap',[iew, iew, jsw, jew, ksw, kew])
            if not t_kmin and dimPb==3:
                C._addBC2Zone(zp,'overlap','BCOverlap',[isw, iew, jsw, jew, ksw, ksw])
            if not t_kmax and dimPb==3:
                C._addBC2Zone(zp,'overlap','BCOverlap',[isw, iew, jsw, jew, kew, kew])

    return None

#==============================================================================
# Generate the full Cartesian mesh (octree/quadtree-based) for IBMs
# IN:
#   tb (tree): geometry tree (IBM bodies)
#   dimPb (2 or 3): problem dimension
#   vmin (int): number of points for each octree level
#   snears (float or list of floats): minimum cell spacing(s) near the bodies
#   dfars (float or list of floats): extent of the domain from the bodies
#   dfarDir (int):
#   tbox (tree): refinement bodies
#   snearsf (float or list of floats) cell spacing(s) to impose inside the refinement bodies
#   check (boolean): if True: write octree.cgns locally
#   to (tree): input octree if already created
#   ext (int): grid extent for overlapping
#   expand (0, 1, 2 or 3): expand minimum cell spacing to other blocks near the bodies
#   octreeMode (0 or 1): octree generation octreeMode. If 0: dfar is exact and snear varies. If 1 it's the opposite
# OUT:
#   t (tree): mesh Tree
#==============================================================================
# only in octree2StructLoc__
def mergeByParent__(zones, parent, sizeMax):
    parent = G.bboxOfCells(parent)
    xmint = Internal.getNodeFromName2(parent,"xmin")[1]
    xmaxt = Internal.getNodeFromName2(parent,"xmax")[1]
    ymint = Internal.getNodeFromName2(parent,"ymin")[1]
    ymaxt = Internal.getNodeFromName2(parent,"ymax")[1]
    zmint = Internal.getNodeFromName2(parent,"zmin")[1]
    zmaxt = Internal.getNodeFromName2(parent,"zmax")[1]

    res = []
    xminAll=[]; yminAll=[]; zminAll=[]; xmaxAll=[]; ymaxAll=[]; zmaxAll=[]
    noz = 0
    for z in zones:
        dimZ = Internal.getZoneDim(z)
        npts = dimZ[1]*dimZ[2]*dimZ[3]
        xmin = C.getValue(z,'CoordinateX',0)
        ymin = C.getValue(z,'CoordinateY',0)
        zmin = C.getValue(z,'CoordinateZ',0)
        xmax = C.getValue(z,'CoordinateX',npts-1)
        ymax = C.getValue(z,'CoordinateY',npts-1)
        zmax = C.getValue(z,'CoordinateZ',npts-1)
        xminAll.append(xmin); xmaxAll.append(xmax)
        yminAll.append(ymin); ymaxAll.append(ymax)
        zminAll.append(zmin); zmaxAll.append(zmax)
        noz += 1

    found=[0]*len(zones)
    for no in range(xmint.shape[0]):
        xmin = xmint[no]; xmax = xmaxt[no]
        ymin = ymint[no]; ymax = ymaxt[no]
        zmin = zmint[no]; zmax = zmaxt[no]
        pool=[]
        for noz in range(len(zones)):
            if found[noz]==0:
                xminz = xminAll[noz]; xmaxz = xmaxAll[noz]
                yminz = yminAll[noz]; ymaxz = ymaxAll[noz]
                zminz = zminAll[noz]; zmaxz = zmaxAll[noz]
                if zminz > zmin-EPSCART and zmaxz < zmax+EPSCART:
                    if yminz > ymin-EPSCART and ymaxz < ymax+EPSCART:
                        if xminz > xmin-EPSCART and xmaxz < xmax+EPSCART:
                            pool.append(zones[noz])
                            found[noz]=1
        if len(pool)> 1:
            res += T.mergeCart(pool, sizeMax=sizeMax)
            del pool
        elif len(pool) == 1: res += pool
    return res

# only in generateIBMMesh and generateCartMesh__
def octree2StructLoc__(o, parento=None, vmin=15, ext=0, optimized=0, sizeMax=4e6, tbOneOver=None):
    sizeMax=int(sizeMax)
    dim = Internal.getZoneDim(o)
    if dim[3] == 'QUAD': dimPb = 2
    elif dim[3] == 'HEXA': dimPb = 3

    if ext == 1: ext = 2
    a = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
    zones = Generator.generator.octree2Struct(a, [vmin])
    c = 1
    for noz in range(len(zones)):
        zones[noz] = C.convertArrays2ZoneNode('cartDummy'+str(c), [zones[noz]])
        c += 1

    ## Rectilinear mesh modifications
    listSavetbOneOver     = []
    if tbOneOver:
        tzones = G.BB(zones)
        if dimPb==2:
            T._addkplane(tzones)
            T._contract(tzones, (0,0,0), (1,0,0), (0,1,0), 0.01)
            tbOneOver=T.addkplane(tbOneOver)
            tbOneOver=T.contract(tbOneOver, (0,0,0), (1,0,0), (0,1,0), 0.01)

        ##RECTILINEAR REGION
        tzones2         = C.newPyTree(['Base', tzones])
        NBases          = len(Internal.getBases(tbOneOver))
        for i in Internal.getBases(tbOneOver):
            checkOneOver = Internal.getNodeByName(i,".Solver#define")
            granularityLocal = 0
            if checkOneOver and  Internal.getNodeByName(checkOneOver, 'granularity'):
                granularityLocal = Internal.getNodeByName(checkOneOver, 'granularity')[1]
            if granularityLocal==1:
                listSavetbOneOverTmp = []
                C._initVars(tzones2, 'cellNOneOver', 1.)
                X_IBM._blankByIBCBodies(tzones2, i, 'nodes', dimPb, cellNName='cellNOneOver')

                for z in Internal.getZones(tzones2):
                    if C.getMinValue(z, 'cellNOneOver')<1:listSavetbOneOverTmp.append(z[0])
                listSavetbOneOver.append(listSavetbOneOverTmp)
            elif granularityLocal==0:
                interDict_scale = X.getIntersectingDomains(G.BB(i), tzones)
                listSavetbOneOverTmp = []
                for j in interDict_scale:
                    for z in interDict_scale[j]: listSavetbOneOverTmp.append(z)
                    listSavetbOneOver.append(listSavetbOneOverTmp)

        del tzones2

    if parento is None:
        ## Rectilinear mesh modifications
        ZONEStbOneOver    = None
        ZONEStbOneOverTmp = None
        if tbOneOver:
            listSavetbOneOverZones = []

            ##Add zones to list of TbOneOver list &
            ##Remove zone from orig list of zones
            for sublist in listSavetbOneOver:
                listSavetbOneOverZonesTmp=[]
                for zTmp in sublist:
                    zAdd=Internal.getNodeFromNameAndType(zones,zTmp,'Zone_t')
                    if zAdd is not None:
                        listSavetbOneOverZonesTmp.append(Internal.copyTree(zAdd))
                        Internal._rmNode(zones, zAdd)
                listSavetbOneOverZones.append(listSavetbOneOverZonesTmp)

            zones = T.mergeCart(zones,sizeMax=sizeMax)
            for sublist in listSavetbOneOverZones:
                ZONEStbOneOverTmp=T.mergeCart(sublist,sizeMax=sizeMax)
                zones           +=ZONEStbOneOverTmp
        else:
            zones = T.mergeCart(zones,sizeMax=sizeMax)
    else:
        ## Rectilinear mesh modifications
        ZONEStbOneOver    = None
        ZONEStbOneOverTmp = None
        eps = 1.e-10
        bbo = G.bbox(parento[0])# 1st octant lower left side
        xmeano=bbo[3]; ymeano=bbo[4]; zmeano=bbo[5]
        # gather zones by parent octant
        if dimPb == 2:
            ZONES=[[],[],[],[]]; noct = 4;
        else:
            ZONES = [[],[],[],[],[],[],[],[]]; noct = 8 ;
        ## Rectilinear mesh modifications
        if tbOneOver:
            ZONEStbOneOver    = []
            ZONEStbOneOverTmp = []
            for i in range(NBases):
                if dimPb==2: ZONEStbOneOverTmp.append([[],[],[],[]])
                else: ZONEStbOneOverTmp.append([[],[],[],[],[],[],[],[]])
                ZONEStbOneOver.append([])

        for z in zones:
            xminz = C.getValue(z,'CoordinateX',0)
            yminz = C.getValue(z,'CoordinateY',0)
            zminz = C.getValue(z,'CoordinateZ',0)
            dimZ = Internal.getZoneDim(z)
            ni = dimZ[1]; nj = dimZ[2]; nk = dimZ[3]
            ind = ni-1 + (nj-1)*ni+(nk-1)*ni*nj
            xmaxz = C.getValue(z,'CoordinateX',ind)
            ymaxz = C.getValue(z,'CoordinateY',ind)
            zmaxz = C.getValue(z,'CoordinateZ',ind)
            # bbz = G.bbox(z)
            # xminz=bbz[0]; yminz=bbz[1]; zminz=bbz[2]
            # xmaxz=bbz[3]; ymaxz=bbz[4]; zmaxz=bbz[5]
            noo = -1
            if dimPb == 3:
                if zmaxz < zmeano+eps:
                    if ymaxz < ymeano+eps:
                        if xmaxz < xmeano+eps: noo=0
                        else: noo=1
                    else:
                        if xmaxz < xmeano+eps: noo=2
                        else: noo=3
                else:
                    if ymaxz < ymeano+eps:
                        if xmaxz < xmeano+eps: noo=4
                        else: noo=5
                    else:
                        if xmaxz < xmeano+eps: noo=6
                        else: noo=7
            else:
                if ymaxz < ymeano+eps:
                    if xmaxz < xmeano+eps: noo=0
                    else: noo=1
                else:
                    if xmaxz < xmeano+eps: noo=2
                    else: noo=3
            if noo > -1:
                if not any(z[0] in sublist for sublist in listSavetbOneOver): ZONES[noo].append(z)
                else:
                    for i in range(NBases):
                        if z[0] in listSavetbOneOver[i]: ZONEStbOneOverTmp[i][noo].append(z)
        #-----------------------------------------------------------------------------
        zones=[]
        for noo in range(noct):
            nzones = len(ZONES[noo])
            if nzones > 1:
                print('Merging %d Cartesian zones of subdomain %d.'%(nzones,noo))
                ZONES[noo] = mergeByParent__(ZONES[noo], parento[noo], sizeMax)
                print('Nb of merged zones: %d.'%len(ZONES[noo]))

            ## Rectilinear mesh modifications
            if tbOneOver:
                for i in range(NBases):
                    nzones = len(ZONEStbOneOverTmp[i][noo])
                    if nzones > 1:
                        print('Merging %d Cartesian zones of subdomain %d - OneOverRegion # %d '%(nzones,noo,i))
                        ZONEStbOneOverTmp[i][noo] = mergeByParent__(ZONEStbOneOverTmp[i][noo], parento[noo], sizeMax)
                        #C.convertPyTree2File(ZONEStbOneOverTmp[i][noo],'t_BaseTmp_'+str(i)+'_'+str(noo)+'.cgns')
                        print('Nb of merged zones - OneOverRegion # %d : %d.' %(i,len(ZONEStbOneOverTmp[i][noo])))
        if dimPb == 3:
            ZONES0 = T.mergeCart(ZONES[0]+ZONES[4],sizeMax=sizeMax)# XM
            ZONES1 = T.mergeCart(ZONES[2]+ZONES[6],sizeMax=sizeMax)# XP
            ZONES2 = T.mergeCart(ZONES[1]+ZONES[5],sizeMax=sizeMax)
            ZONES3 = T.mergeCart(ZONES[3]+ZONES[7],sizeMax=sizeMax)
            del ZONES
            ZONES0 = T.mergeCart(ZONES0+ZONES1,sizeMax=sizeMax)
            del ZONES1
            ZONES2 = T.mergeCart(ZONES2+ZONES3,sizeMax=sizeMax)
            del ZONES3
            zones = T.mergeCart(ZONES0+ZONES2,sizeMax=sizeMax)
            del ZONES0
            del ZONES2

            ## Rectilinear mesh modifications
            if tbOneOver:
                for i in range(NBases):
                    ZONEStbOneOver[i] = T.mergeCart(ZONEStbOneOverTmp[i][0]+ZONEStbOneOverTmp[i][1]+ZONEStbOneOverTmp[i][2]+ \
                                                    ZONEStbOneOverTmp[i][3]+ZONEStbOneOverTmp[i][4]+ZONEStbOneOverTmp[i][5]+ \
                                                    ZONEStbOneOverTmp[i][6]+ZONEStbOneOverTmp[i][7], sizeMax=sizeMax)
                    #C.convertPyTree2File(ZONEStbOneOver[i],'t_Base_'+str(i)+'.cgns')
                    zones +=ZONEStbOneOver[i]
                del ZONEStbOneOver
                del ZONEStbOneOverTmp
        else: # dim=2
            ZONES[0] = T.mergeCart(ZONES[0]+ZONES[2],sizeMax=sizeMax)# XM
            ZONES[1] = T.mergeCart(ZONES[1]+ZONES[3],sizeMax=sizeMax)# XP
            ZONES=ZONES[0:2]
            zones = T.mergeCart(ZONES[0]+ZONES[1],sizeMax=sizeMax)
            del ZONES

            ## Rectilinear mesh modifications
            if tbOneOver:
                for i in range(NBases):
                    ZONEStbOneOver[i] = T.mergeCart(ZONEStbOneOverTmp[i][0]+ZONEStbOneOverTmp[i][1]+ \
                                                    ZONEStbOneOverTmp[i][2]+ZONEStbOneOverTmp[i][3], sizeMax=sizeMax)
                    #C.convertPyTree2File(ZONEStbOneOver[i],'t_Base_'+str(i)+'.cgns')
                    zones +=ZONEStbOneOver[i]
            del ZONEStbOneOver
            del ZONEStbOneOverTmp
    print('After merging: nb Cartesian zones=%d (ext. =%d).'%(len(zones),ext))
    # Cas ext=-1, ne fait pas les extensions ni les BCs ou raccords

    if optimized == -1:
        if dimPb == 3: ratios = [[2,2,2],[4,4,4],[8,8,8],[16,16,16]]
        else: ratios = [[2,2,1],[4,4,1],[8,8,1],[16,16,1]]
        zones = X.connectMatch(zones, dim=dimPb)
        for ratio0 in ratios:
            zones = X.connectNearMatch(zones, ratio=ratio0, dim=dimPb)

    if ext == -1: return zones

    if ext > 0:
        coords = C.getFields(Internal.__GridCoordinates__, zones, api=3)
        coords,rinds = Generator.extendCartGrids(coords, ext=ext, optimized=optimized, extBnd=0)
        C.setFields(coords, zones, 'nodes')
        for noz in range(len(zones)):
            Internal.newRind(value=rinds[noz], parent=zones[noz])
    # Creation des zones du pyTree
    for z in zones: z[0] = C.getZoneName('cart')
    if ext == 0:
        if dimPb == 3: ratios = [[2,2,2],[4,4,4],[8,8,8],[16,16,16]]
        else: ratios = [[2,2,1],[4,4,1],[8,8,1],[16,16,1]]
        zones = X.connectMatch(zones, dim=dimPb)
        for ratio0 in ratios:
            zones = X.connectNearMatch(zones,ratio=ratio0,dim=dimPb)
        return zones
    else:
        bbox0 = G.bbox(o)
        _addBCOverlaps(zones, bbox0)
    return zones

# only in generateIBMMesh and generateIBMMesh
def buildParentOctrees__(o, tb, dimPb=3, vmin=15, snears=0.01, snearFactor=1., dfars=10., dfarDir=0,
                         tbox=None, snearsf=None, to=None, octreeMode=0):
    nzones0 = Internal.getZoneDim(o)[2]
    if nzones0 < 1000: return None

    parento = buildOctree(tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=snearFactor, dfars=dfars, dfarDir=dfarDir,
                          tbox=tbox, snearsf=snearsf, to=to, expand=-1, balancing=0, octreeMode=octreeMode)

    bbo = G.bbox(parento)
    xmino=bbo[0]; xmaxo=bbo[3]; xmeano=0.5*(xmino+xmaxo)
    ymino=bbo[1]; ymaxo=bbo[4]; ymeano=0.5*(ymino+ymaxo)
    zmino=bbo[2]; zmaxo=bbo[5]; zmeano=0.5*(zmino+zmaxo)
    dx = xmeano-xmino; dy = ymeano-ymino; dz = zmeano-zmino
    eps = 1.e-10

    OCTREEPARENTS = None

    if dimPb == 2:
        OCTREEPARENTS=[]
        for ym in [ymino,ymeano]:
            for xm in [xmino,xmeano]:
                C._initVars(parento,'centers:tag',1.)
                C._initVars(parento,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))
                parento2 = P.selectCells2(parento,'centers:tag')
                OCTREEPARENTS.append(parento2)
    else:
        OCTREEPARENTS=[]
        for zm in [zmino,zmeano]:
            for ym in [ymino,ymeano]:
                for xm in [xmino,xmeano]:
                    C._initVars(parento,'centers:tag',1.)
                    C._initVars(parento,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                    C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))
                    C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateZ}>%g)*({centers:CoordinateZ}<%g)'%(zm-eps,zm+dz+eps))
                    parento2 = P.selectCells2(parento,'centers:tag')
                    OCTREEPARENTS.append(parento2)
    return OCTREEPARENTS

def generateIBMOctree(tb, dimPb=3, vmin=15, snears=0.01, dfars=10., dfarDir=0,
                      tbox=None, snearsf=None, check=False, to=None,
                      expand=3, octreeMode=0, optimized=1):
    ## tbox = tbox(area refinement)[legacy tbox] + tb(area for one over)[tbOneOver] + tb(zones to keep as F1)[tbF1]
    ## here we divide tbox into tbOneOver (rectilinear region)  & tbF1 (WM F1 approach region) & tbox; tbox will henceforth only consist of the area that will be refined.
    ## Note: tb(zones to keep as F1)[tbF1] is still in development, is experimental, & subject to major/minor changes with time. Please use with a lot of caution & see A.Jost @ DAAA/DEFI [28/08/2024] as
    ##       there is no non-regression test case yet available.
    tbOneOverF1 = None
    tbOneOver   = None
    tbF1        = None
    tbRM        = None
    if tbox:
        tbRM        = Internal.getNodesFromNameAndType(tbox, '*RM*'     , 'CGNSBase_t')
        tbOneOver   = Internal.getNodesFromNameAndType(tbox, '*OneOver*', 'CGNSBase_t')
        tbF1        = Internal.getNodesFromNameAndType(tbox, '*KeepF1*' , 'CGNSBase_t')
        tbOneOverF1 = tbOneOver+tbF1
        tbox        = Internal.rmNodesByName(Internal.rmNodesByName(Internal.rmNodesByName(tbox, '*OneOver*'), '*KeepF1*'), '*RM*')
        if len(Internal.getBases(tbox))==0: tbox=None

    # Octree identical on all procs
    if to is not None:
        if isinstance(to, str):
            o = C.convertFile2PyTree(to)
            o = Internal.getZones(o)[0]
        else:
            o = Internal.getZones(to)[0]
        parento = None
    else:
        o = buildOctree(tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=1., dfars=dfars, dfarDir=dfarDir,
                        tbox=tbox, snearsf=snearsf, to=to, expand=expand, octreeMode=octreeMode, optimized=optimized)

        # build parent octree 3 levels higher
        # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
        parento = buildParentOctrees__(o, tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=4., dfars=dfars, dfarDir=dfarDir,
                                       tbox=tbox, snearsf=snearsf, to=to, octreeMode=octreeMode)

    # adjust the extent of the box defining the symmetry plane if in tb
    baseSYM = Internal.getNodeFromName1(tb,"SYM")
    dir_sym = 0; X_SYM = 0.; coordsym =  None; symmetry=0
    if baseSYM is not None:
        symmetry=1; symplane = []
        for zsym in Internal.getZones(baseSYM):
            if C.getMaxValue(zsym,'centers:cellN')>0.:
                symplane.append(zsym)
        [xmin,ymin,zmin,xmax,ymax,zmax] = G.bbox(symplane)
        if abs(xmax-xmin) < 1e-6:
            coordsym = 'CoordinateX'; dir_sym=1; X_SYM = xmin
        elif abs(ymax-ymin) < 1e-6:
            coordsym = 'CoordinateY'; dir_sym=2; X_SYM = ymin
        elif abs(zmax-zmin)<1e-6:
            coordsym = 'CoordinateZ'; dir_sym=3; X_SYM = zmin

        if octreeMode==1:
            Internal._rmNodesFromType(baseSYM,"Zone_t")
            [xmin,ymin,zmin,xmax,ymax,zmax] = G.bbox(o)
            L = 0.5*(xmax+xmin); eps = 0.2*L
            xmin = xmin-eps; ymin = ymin-eps; zmin = zmin-eps
            xmax = xmax+eps; ymax = ymax+eps; zmax = zmax+eps
            if dir_sym==1: xmax=X_SYM
            elif dir_sym==2: ymax=X_SYM
            elif dir_sym==3: zmax = X_SYM
            a = D.box((xmin,ymin,zmin),(xmax,ymax,zmax))
            C._initVars(a,'{centers:cellN}=({centers:%s}>-1e-8)'%coordsym)
            baseSYM[2]+=a
        else: pass
        if coordsym is not None:
            to = C.newPyTree(["OCTREE",o])
            bodies = [Internal.getZones(baseSYM)]
            BM2 = numpy.ones((2,1),dtype=Internal.E_NpyInt)
            to = X.blankCellsTri(to, bodies, BM2, blankingType='center_in', cellNName='cellN')
            to = P.selectCells(to,'{centers:cellN}>0.')
            o = Internal.getZones(to)[0]

    if tbRM:
        to       = C.newPyTree(["OCTREE",o])
        bodies   = [Internal.getBases(tbRM)]
        nBasesRM = len(bodies)
        BM2      = numpy.ones((2,nBasesRM),dtype=Internal.E_NpyInt)
        if dimPb ==2:
            XRAYDIM1 = 20000; XRAYDIM2 = XRAYDIM1
            to = X.blankCells(to, bodies, BM2, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=dimPb, cellNName='cellN')
        else:
            to = X.blankCellsTri(to, bodies, BM2, blankingType='node_in', cellNName='cellN')
        to = P.selectCells(to,'{cellN}<=0.',strict=1)
        C._initVars(to,'{cellN}=1')
        o = Internal.getZones(to)[0]

    if Cmpi.rank==0 and check: C.convertPyTree2File(o, 'octree.cgns')
    return o, parento, tbOneOver, tbF1, tbOneOverF1, [symmetry, dir_sym, X_SYM]

# main function
def generateIBMMesh(tb, dimPb=3, vmin=15, snears=0.01, dfars=10., dfarDir=0,
                    tbox=None, snearsf=None, check=False, to=None,
                    ext=2, optimized=1, expand=3, octreeMode=0):
    """Generates the full Cartesian mesh for IBMs."""
    import KCore.test as test
    # refinementSurfFile: surface meshes describing refinement zones
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox

    o, parento, tbOneOver, tbF1, tbOneOverF1, infoSym = generateIBMOctree(tb, dimPb=dimPb, vmin=vmin, snears=snears, dfars=dfars, dfarDir=dfarDir,
                                                                          tbox=tbox, snearsf=snearsf, check=check, to=to,
                                                                          expand=expand, octreeMode=octreeMode, optimized=optimized)

    # Split octree
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[Cmpi.rank]
    del o

    # fill vmin + merge in parallel
    if optimized == -1:
        res = octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=optimized, parento=parento, sizeMax=1000000, tbOneOver=tbOneOverF1)
    else:
        res = octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000, tbOneOver=tbOneOverF1)

    del p
    if parento is not None:
        for po in parento: del po

    #Split les zones qui possede des faces avec racord 1:1, 1:2  et 2:1 (zones voisines avec 3 niveau de resolution different)
    if optimized==-1:
        t1 = C.newPyTree(['CARTESIAN', res])
        zones = Internal.getZones(t1)
        t = splitZoneForConservative(t1,dimPb=dimPb )
    else:
        t = C.newPyTree(['CARTESIAN', res])

    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%Cmpi.rank
    Cmpi._setProc(t, Cmpi.rank)

    C._addState(t, 'EquationDimension', dimPb)

    # Keep F1 regions - for F1 & F42 synergy
    if tbF1:
        tbbBTmp = G.BB(tbF1)
        interDict_scale = X.getIntersectingDomains(tbbBTmp, t)
        for kk in interDict_scale:
            for kkk in interDict_scale[kk]:
                z = Internal.getNodeFromName(t, kkk)
                Internal._createUniqueChild(z, '.Solver#defineTMP', 'UserDefinedData_t')
                Internal._createUniqueChild(Internal.getNodeFromName1(z, '.Solver#defineTMP'), 'SaveF1', 'DataArray_t', value=1)

    # Add xzones for ext
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)

    # Turn Cartesian grid into a rectilinear grid
    test.printMem(">>> cart grids --> rectilinear grids [start]")
    listDone = []
    if tbOneOver:
        test.printMem(">>> cart grids --> rectilinear grids [inside]")
        tbb = G.BB(t)

        if dimPb == 2:
            T._addkplane(tbb)
            T._contract(tbb, (0,0,0), (1,0,0), (0,1,0), 0.01)
            tbOneOver = T.addkplane(tbOneOver)
            tbOneOver = T.contract(tbOneOver, (0,0,0), (1,0,0), (0,1,0), 0.01)

        ## RECTILINEAR REGION
        tzones2  = Internal.copyTree(tbb)
        for i in Internal.getBases(tbOneOver):
            checkOneOver = Internal.getNodeByName(i,".Solver#define")
            if checkOneOver:
                granularityLocal = 0
                if Internal.getNodeByName(checkOneOver, 'granularity'):
                    granularityLocal = Internal.getNodeByName(checkOneOver, 'granularity')[1]
                oneoverX = Internal.getNodeByName(checkOneOver, 'dirx')
                oneoverX = Internal.getValue(oneoverX)
                oneoverY = Internal.getNodeByName(checkOneOver, 'diry')
                oneoverY = Internal.getValue(oneoverY)
                oneoverZ = Internal.getNodeByName(checkOneOver, 'dirz')
                oneoverZ = Internal.getValue(oneoverZ)
                ## Select regions that need to be coarsened
                if granularityLocal == 1:
                    C._initVars(tzones2, 'cellNOneOver', 1.)
                    X_IBM._blankByIBCBodies(tzones2, i, 'nodes', dimPb, cellNName='cellNOneOver')

                    for z in Internal.getZones(tzones2):
                        if C.getMinValue(z, 'cellNOneOver')<1:
                            if z[0] not in listDone:
                                zLocal = Internal.getNodeFromName(t,z[0])
                                T._oneovern(zLocal, (oneoverX,oneoverY,oneoverZ));
                                ### Avoid a zone to be coarsened twice
                                listDone.append(z[0])
                elif granularityLocal == 0:
                    interDict_scale = X.getIntersectingDomains(G.BB(i), tbb)
                    ## Avoid a zone to be coarsened twice
                    for j in interDict_scale:
                        for z in interDict_scale[j]:
                            if z not in listDone:
                                zLocal = Internal.getNodeFromName(t,z)
                                T._oneovern(zLocal, (oneoverX,oneoverY,oneoverZ));
                                listDone.append(z)

        del tzones2

    test.printMem(">>> cart grids --> rectilinear grids [end]")

    symmetry, dir_sym, X_SYM = infoSym
    if symmetry == 0: extBnd = 0
    else: extBnd = ext-1 # nb de ghost cells = ext-1
    if ext > 0:
        zones = Internal.getZones(t)
        coords = C.getFields(Internal.__GridCoordinates__, zones, api=3)
        coords, rinds = Generator.extendCartGrids(coords, ext=ext, optimized=optimized, extBnd=extBnd)
        C.setFields(coords, zones, 'nodes')
        for noz in range(len(zones)):
            Internal.newRind(value=rinds[noz], parent=zones[noz])

    Cmpi._rmXZones(t)

    coords = None; zones = None

    if symmetry == 0:
        _addBCOverlaps(t, bbox=bb)
        _addExternalBCs(t, bbox=bb, dimPb=dimPb)
    else:
        _addBCsForSymmetry(t, bbox=bb, dimPb=dimPb, dir_sym=dir_sym, X_SYM=X_SYM, depth=ext-1)

    if dimPb == 2:
        dz = 0.01
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    return t

# alias generateIBMMesh new version
generateIBMMeshPara = generateIBMMesh
#==============================================================================
# Generate the full octree (3D) or quadtree (2D) for IBMs
# IN:
#   tb (tree): geometry tree (IBM bodies)
#   dimPb (2 or 3): problem dimension
#   vmin (int): number of points for each octree level
#   snears (float or list of floats): minimum cell spacing(s) near the bodies
#   snearFactor (float): snear multiplicator
#   dfars (float or list of floats): extent of the domain from the bodies
#   dfarDir (int):
#   tbox (tree): refinement bodies
#   snearsf (float or list of floats) cell spacing(s) to impose inside the refinement bodies
#   to (tree): input octree if already created
#   expand (0, 1, 2 or 3): expand minimum cell spacing to other blocks near the bodies
#   octreeMode (0 or 1): octree generation octreeMode. If 0: dfar is exact and snear varies. If 1 it's the opposite
# OUT:
#   t (tree): mesh Tree
#==============================================================================
# only in buildOctree
def addRefinementZones__(o, tb, tbox, snearsf, vmin, dim):
    tbSolid = Internal.rmNodesByName(tb, 'IBCFil*')

    if dim == 2:
        tbSolid = T.addkplane(tbSolid)
        tbSolid = T.translate(tbSolid, (0,0,-0.5))

        tbox = T.addkplane(tbox)
        tbox = T.translate(tbox, (0,0,-0.5))

    boxes = []
    for b in Internal.getBases(tbox):
        boxes.append(Internal.getNodesFromType1(b, 'Zone_t'))
    if snearsf is not None:
        if not isinstance(snearsf, list): snearsf = len(boxes)*[snearsf]
        if len(boxes) != len(snearsf):
            raise ValueError('addRefinementZones: Number of refinement bodies is not equal to the length of snearsf list.')
        for i in range(len(snearsf)):
            snearsf[i] = snearsf[i]*(vmin-1)
    else:
        snearsf=[]
        for sbox in boxes:
            for s in Internal.getZones(sbox):
                sdd = Internal.getNodeFromName1(s, ".Solver#define")
                if sdd is not None:
                    snearl = Internal.getNodeFromName1(sdd, "snear")
                    if snearl is not None:
                        snearl = Internal.getValue(snearl)
                        snearsf.append(snearl*(vmin-1))

    to = C.newPyTree(['Base', o])
    end = 0
    G._getVolumeMap(to)
    volmin0 = C.getMinValue(to, 'centers:vol')
    # volume minimum au dela duquel on ne peut pas raffiner
    snears = Internal.getNodesFromName(tbSolid, 'snear')
    if snears != []: minSneartb= min(snears, key=lambda x: x[1])[1][0]*(vmin-1)
    else: minSneartb= min(snearsf)
    volmin0 *= (min(1.0,min(snearsf)/minSneartb)**(dim))
    while end == 0:
        # Do not refine inside obstacles
        C._initVars(to,'cellN',1.)
        to = X_IBM.blankByIBCBodies(to, tbSolid, 'nodes', 3)
        to = C.node2Center(to, 'cellN')
        Internal._rmNodesFromName(to, Internal.__FlowSolutionNodes__)
        C._initVars(to, '{centers:cellN}=({centers:cellN}>0.1)')
        C._initVars(to, '{centers:cellNBody}={centers:cellN}')
        nob = 0
        C._initVars(to, 'centers:indicator', 0.)
        for box in boxes:
            volmin2 = 1.09*(snearsf[nob])**(dim)
            C._initVars(to,'cellN',1.)
            tboxl = C.newPyTree(['BOXLOC']); tboxl[2][1][2] = box
            to = X_IBM.blankByIBCBodies(to, tboxl, 'nodes', 3)
            to = C.node2Center(to, 'cellN')
            C._initVars(to, '{centers:cellN}=({centers:cellN}>0.9)')
            Internal._rmNodesFromName(to, Internal.__FlowSolutionNodes__)
            fact = 1.1
            while C.getMinValue(to, 'centers:cellN') == 1 and fact < 10.:
                print("Info: addRefinementZones: tbox too small - increase tbox by fact = %2.1f"%(fact))
                box2 = T.scale(box, fact)
                tboxl[2][1][2] = box2
                to = X_IBM.blankByIBCBodies(to, tboxl, 'centers', 3)
                fact += 0.1

            C._initVars(to,'{centers:indicator}=({centers:indicator}>0.)+({centers:indicator}<1.)*logical_and({centers:cellN}<0.001, {centers:vol}>%g)'%volmin2)
            nob += 1

        end = 1
        C._initVars(to,'{centers:indicator}={centers:indicator}*({centers:cellNBody}>0.)*({centers:vol}>%g)'%volmin0)

        if  C.getMaxValue(to, 'centers:indicator') == 1.:
            end = 0
            # Maintien du niveau de raffinement le plus fin
            o = Internal.getZones(to)[0]
            o = G.adaptOctree(o, 'centers:indicator', balancing=2)
            to[2][1][2] = [o]
            G._getVolumeMap(to)

    return Internal.getNodeFromType2(to, 'Zone_t')

def buildOctree(tb, dimPb=3, vmin=15, snears=0.01, snearFactor=1., dfars=10., dfarDir=0,
                tbox=None, snearsf=None, to=None, balancing=2, expand=2, octreeMode=0, optimized=1):

    """Builds an octree from the surface definitions."""

    surfaces=[]; dfarListL=[]; snearso=[]

    # list of dfars
    bodies = Internal.getZones(tb)
    if not isinstance(dfars, list):
        dfarList = [dfars*1.]*len(bodies)
        for c, z in enumerate(bodies):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfarList[c] = Internal.getValue(n)*1.
    else:
        if len(bodies) != len(dfars): raise ValueError('buildOctree: Number of bodies is not equal to the size of dfars.')
        dfarList = list(dfars)
    dfar = -1

    # List of snears
    if not isinstance(snears, list):
        snears = [snears*1.]*len(bodies)
        for c, z in enumerate(bodies):
            n = Internal.getNodeFromName2(z, 'snear')
            if n is not None: snears[c] = Internal.getValue(n)*1.
    else:
        if len(bodies) != len(snears): raise ValueError('buildOctree: Number of bodies is not equal to the size of snears.')

    dxmin0 = 1.e10
    for c, z in enumerate(bodies):
        if dfarList[c] > -1: #body snear is only considered if dfar_loc > -1
            dhloc = snears[c]*(vmin-1)*snearFactor
            surfaces.append(z)
            snearso.append(dhloc)
            dfarListL.append(dfarList[c])
            dxmin0 = min(dxmin0, dhloc)

    if to is not None:
        o = Internal.getZones(to)[0]
    else:
        o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarListL, balancing=balancing, dfarDir=dfarDir, octreeMode=octreeMode)
        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if dxmin < 0.65*dxmin0:
            snearso = [2.*i for i in snearso]
            o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarListL, balancing=balancing, dfarDir=dfarDir, octreeMode=octreeMode)
        # Adaptation avant expandLayer (pour corriger eventuellement les sauts de maille)
        if tbox is not None:
            o = addRefinementZones__(o, tb, tbox, snearsf, vmin, dimPb)
            C._rmVars(o, ['centers:indicator', 'centers:cellN', 'centers:vol', 'centers:cellNBody'])

        if expand == 0:
            G._expandLayer(o, level=0, corners=1, balancing=1)
        elif expand == 1:
            vmint = 31
            if vmin < vmint:
                if Cmpi.rank==0: print('buildOctree: octree finest level expanded (expandLayer activated).')
                to = C.newPyTree(['Base',o])
                to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
                C._initVars(o, "centers:indicator", 0.)
                cellN = C.getField("centers:cellN", to)[0]
                octreeA = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
                indic = C.getField("centers:indicator", o, api=1)[0]
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, 0, 2)
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 1, 0, 2) # CB
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 2, 0, 2) # CB
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 3, 0, 2) # CB
                indic = Converter.addVars([indic,cellN])
                indic = Converter.initVars(indic, "{indicator}={indicator}*({cellN}>0.)")
                octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
                o = C.convertArrays2ZoneNode(o[0], [octreeA])

            to = C.newPyTree(['Base',o])
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
            indic = C.getField("centers:cellN",to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
            indic = Converter.initVars(indic, 'indicator', 0., api=1)
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0,0,1)
            indic = Converter.extractVars(indic, ["indicator"])
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

        else: #expand = 2, 3 or 4
            corner = 0
            if expand==2 and optimized==-1: corner = 1

            to = C.newPyTree(['Base',o])
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
            indic = C.getField("centers:indicator", o, api=1)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 3)
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

            if expand > 2: # expand minimum + 1 couche propagee
                # passe 2
                to = C.newPyTree(['Base',o])
                to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
                C._initVars(o, "centers:indicator", 0.)
                cellN = C.getField("centers:cellN", to)[0]
                octreeA = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
                indic = C.getField("centers:indicator", o, api=1)[0]
                indic = Converter.addVars([indic,cellN])
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 4)
                octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
                o = C.convertArrays2ZoneNode(o[0], [octreeA])
                # fin passe 2

            if expand > 3:# expand minimum + 2 couche propagee
                # passe 3
                to = C.newPyTree(['Base',o])
                to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
                C._initVars(o, "centers:indicator", 0.)
                cellN = C.getField("centers:cellN", to)[0]
                octreeA = C.getFields(Internal.__GridCoordinates__, o, api=1)[0]
                indic = C.getField("centers:indicator", o, api=1)[0]
                indic = Converter.addVars([indic,cellN])
                indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 6)
                octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
                o = C.convertArrays2ZoneNode(o[0], [octreeA])
                # fin passe 3

        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        C._rmVars(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if Cmpi.rank == 0: print('Minimum spacing of Cartesian mesh= %f (targeted %f)'%(dxmin/(vmin-1),dxmin0/(vmin-1)))

        nelts = Internal.getZoneDim(o)[2]
        if nelts > 20000:
            print('Warning: number of zones (%d) on rank %d is high (block merging might last a long time).'%(nelts, Cmpi.rank))

    #sym = Internal.getNodeFromPath(tb, 'SYM')
    #if sym is not None: # sym plan exists
    #    tbsym = G.BB(sym)
    #    xmin = C.getMinValue(tbsym,'CoordinateX'); xmax = C.getMaxValue(tbsym,'CoordinateX')
    #    ymin = C.getMinValue(tbsym,'CoordinateY'); ymax = C.getMaxValue(tbsym,'CoordinateY')
    #    zmin = C.getMinValue(tbsym,'CoordinateZ'); zmax = C.getMaxValue(tbsym,'CoordinateZ')
    #    o = P.selectCells(o, '({CoordinateX}<=%g)*({CoordinateX}>=%g)'%(xmax,xmin), strict=1)
    #    o = P.selectCells(o, '({CoordinateY}<=%g)*({CoordinateY}>=%g)'%(ymax,ymin), strict=1)
    #    o = P.selectCells(o, '({CoordinateZ}<=%g)*({CoordinateZ}>=%g)'%(zmax,zmin), strict=1)

    return o

#==============================================================================
#
#==============================================================================
def createRefinementBodies(tb, dimPb=3, hmod=0.01, hmax=None, pointsPerUnitLength=None):
    """Creates refinement bodies from the immersed boundaries to extend the finest resolution in the fluid domain."""
    import Geom.IBM as D_IBM
    import Geom.Offset

    if pointsPerUnitLength is None:
        print("Info: createRefinementBodies: pointsPerUnitLength is None, using default values (25 for 3D or 1000 for 2D).")
        pointsPerUnitLength = 25 if dimPb == 3 else 1000

    refinementBodies = []

    tb = Internal.rmNodesFromName(tb, "SYM")
    tb = Internal.rmNodesFromName(tb, "*_sym")

    if hmax is None:
        snears = Internal.getNodesFromName(tb, 'snear')
        hmax   = min(snears, key=lambda x: x[1])[1][0]

    for z in Internal.getZones(tb):
        snear = Internal.getNodeFromName(z, 'snear')[1]
        zname = z[0]
        if snear <= 1.5*hmax:
            if dimPb == 2:
                z2 = D_IBM.closeContour(z)
            else:
                z2 = D_IBM.closeSurface(z)

            a = Geom.Offset.offsetSurface(z2, offset=hmod, pointsPerUnitLength=pointsPerUnitLength, algo=0, dim=dimPb)

            a = T.splitConnexity(a)
            a = max([za for za in Internal.getZones(a)], key=lambda za: len(Internal.getNodeFromName(za, 'CoordinateX')[1]))

            if dimPb == 2:
                a = T.reorder(a, (1,2,3))
                a = D.uniformize(a, N=1000)

            D_IBM._setSnear(a, snear)
            if dimPb == 2:
                # 2D blanking uses xray method : each body must be in a separate zone so
                # that the blanking remains correct in overlapping regions
                refinementBodies = refinementBodies + [zname,Internal.getZones(a)]
            else:
                refinementBodies = refinementBodies + Internal.getZones(a)

    if dimPb == 2:
        refinementBodies = C.newPyTree(refinementBodies)
    else:
        refinementBodies = C.newPyTree(['Base',refinementBodies])
        refinementBodies = T.reorder(refinementBodies, (-1,))
        refinementBodies = C.convertArray2Tetra(refinementBodies)

    return refinementBodies

def _projectMeshSize(t, NPas=10, span=1, dictNz=None, isCartesianExtrude=False):
    """Predicts the final size of the mesh when extruding 2D to 3D in the z-direction.
    Usage: _projectMeshSize(t, NPas, span, dictNz, isCartesianExtrude)"""
    NP             = Cmpi.size
    rank           = Cmpi.rank
    NPTS           = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS         = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NPTS_noghost   = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    NCELLS_noghost = numpy.zeros(NP, dtype=Internal.E_NpyInt)
    if isinstance(t, str):
        h = Filter.Handle(t)
        t = h.loadFromProc(loadVariables=False)
        h._loadVariables(t, var=['CoordinateX'])


    for z in Internal.getZones(t):
        name_zone = z[0]
        if not isCartesianExtrude:
            if dictNz is not None: Nk = int(dictNz[name_zone])
            else: Nk = NPas-1
        else:
            h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
            NPas_local = int(round(span/h)/4)
            if NPas_local < 2:
                print("WARNING:: MPI rank %d || Zone %s has Nz=%d and is being clipped to Nz=2"%(rank,z[0],NPas_local))
                NPas_local = 2
            Nk = NPas_local
        Nk += 1
        NPTS[rank]           += C.getNPts(z)/2*(Nk+4)
        NCELLS[rank]         += C.getNCells(z)*(Nk+3)
        NPTS_noghost[rank]   += C.getNPts(C.rmGhostCells(z, z, 2))*Nk
        NCELLS_noghost[rank] += C.getNCells(C.rmGhostCells(z, z, 2))*(Nk-1)
    NPTS             = Cmpi.allreduce(NPTS  ,op=Cmpi.SUM)
    NCELLS           = Cmpi.allreduce(NCELLS,op=Cmpi.SUM)
    NPTS_noghost     = Cmpi.allreduce(NPTS_noghost  ,op=Cmpi.SUM)
    NCELLS_noghost   = Cmpi.allreduce(NCELLS_noghost,op=Cmpi.SUM)
    if rank ==0:
        print('Projected mesh size with ghost: {} million points & {} million cells'.format(numpy.sum(NPTS)/1.e6,numpy.sum(NCELLS)/1.e6))
        print('Projected mesh size without ghost: {} million points & {} million cells'.format(numpy.sum(NPTS_noghost)/1.e6,numpy.sum(NCELLS_noghost)/1.e6))
    return None

def extrudeCartesianZDir(t, tb, check=False, extrusion="cart", dz=0.01, NPas=10, span=1 , Ntranche=1,
                         dictNz=None, ific=2, isCartesianExtrude=False, isAutoPeriodic=False, nghost=0):
    """Extrudes a 2D IBM grid and geoemtry. The extraction is done in the z-direction.
        Usage: extrudeCartesianZDir(t, tb, check, extrusion, dz, NPas, span, Ntranche, dictNz, ific, isCartesianExtrude, isAutoPeriodic, nghost)"""

    ## isAutoPeriodic=True is prefered over isAutoPeriodic=False when extruding a Cartesian mesh.
    ## nghost should be 2. It is currently 0 (default value) to pass the non-regression tests. This function was initially developed for
    ## an application and its default value of 0 is kept for reproduciability reasons for that original test case. nghost is only used when isAutoPeriodic=False.
    ## isAutoPeridioc=True assumes 2 ghost cells.

    if  extrusion == "cart": perio = span/float(Ntranche)
    else:                    perio = span/180.*math.pi/float(Ntranche)

    ## Found in original function. Commented here as it is related to interpolations
    ## and not geometric extrusion. This functionality is now in Connector/IBM.py - it is has been
    ## deemed to be a more appropriate location.
    #for z in Internal.getZones(t):
    #    cellN = Internal.getNodeFromName(z, "cellN")[1]
    #    sh    = numpy.shape(cellN)
    #    # modif CellN pour filtrer les cellule solide de l'interp chimere
    #    # et autoriser les ghost comme donneuse
    #    for k in range(sh[2]):
    #        for j in range(sh[1]):
    #            for i in range(sh[0]):
    #                if  cellN[i,j,k] != 0:  cellN[i,j,k] =1

    #Dim = 3D
    c = 0
    for tree in [t,tb]:
        for dim in Internal.getNodesFromName(tree,'EquationDimension'): dim[1]=3

        dz_loc={}; Nk={}
        if c ==0:
            Nk_min = int(10e09)
            for z in Internal.getZones(tree):
                name_zone = z[0]
                if not isCartesianExtrude:
                    if dictNz is not None:
                        Nk[z[0]]     = int(dictNz[name_zone])
                        dz_loc[z[0]] = span/float(Ntranche*Nk[z[0]])
                    else:
                        Nk[z[0]]     = NPas-1
                        dz_loc[z[0]] = dz
                else:
                    h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
                    NPas_local = int(round(span/h))
                    if NPas_local<4:
                        print("WARNING:: Zone %s has Nz=%d and is being clipped to Nz=4"%(z[0],NPas_local))
                        NPas_local=4
                    Nk[z[0]]     = NPas_local
                    dz_loc[z[0]] = span/float(Ntranche*Nk[z[0]])
                if Nk[z[0]] < Nk_min: Nk_min =Nk[z[0]]
        else:
            for z in Internal.getZones(tree):
                Nk[z[0]]     =  Nk_min
                dz_loc[z[0]] = span/float(Ntranche*Nk_min)

        ##Adding ghost cells in K direction
        for z in Internal.getZones(tree):
            Nk[z[0]] += 2*ific-1+c  # -1, for the mesh (already added in the mesh generation) || need to remove this as it is not the case for the geometry (tb)

        for z in Internal.getZones(tree):
            yy_2d = Internal.getNodeFromName(z, "CoordinateY")[1]
            zz_2d = Internal.getNodeFromName(z, "CoordinateZ")[1]
            sh_2d = numpy.shape(yy_2d)

            T._addkplane(z   ,N=Nk[z[0]])

            zz = Internal.getNodeFromName(z, "CoordinateZ")[1]
            yy = Internal.getNodeFromName(z, "CoordinateY")[1]
            sh = numpy.shape(yy)
            r  = numpy.sqrt( zz**2+yy**2)

            perio_loc = perio
            if c==1: period_loc= perio*1.5 #is periodic_loc a typo or a new variable that is not used anywhere else?

            if  len(sh_2d) == 1: #NGON
                if  extrusion == "cart":
                    for k in range( sh[1]):
                        zz[:,k] = (k-ific)* dz_loc[z[0]]
                else:
                    theta0 = numpy.arctan(zz_2d/yy_2d)
                    for k in range( sh[1]):
                        shift = perio_loc/float(sh[1]-5.)*(k-ific)
                        zz[:,k] = r[:,0]*numpy.sin(theta0[:] + shift)
                        yy[:,k] = r[:,0]*numpy.cos(theta0[:] + shift)
            else:
                if  extrusion == "cart":
                    for k in range(sh[2]):
                        zz[:,:,k] = (k-ific)* dz_loc[z[0]]
                else:
                    theta0  = numpy.arctan(zz_2d/yy_2d)
                    for k in range(sh[2]):
                        shift = perio_loc/float(sh[2]-5.)*(k-ific)
                        zz[:,:,k] = r[:,:,0]*numpy.sin(theta0[:,:,0] + shift)
                        yy[:,:,k] = r[:,:,0]*numpy.cos(theta0[:,:,0] + shift)

        if c==1: break
        c += 1
        for z in Internal.getZones(tree):
            zdim = Internal.getValue(z)

            # Modif rind cellule GH en kmin et kmax
            rind = Internal.getNodeFromName1(z, "Rind")[1]
            rind[4] = ific; rind[5] = ific
            # Modif range BC
            BCs = Internal.getNodesFromType2(z, "BC_t")
            for bc in BCs:
                ptrg = Internal.getNodeFromName(bc, "PointRange")[1]
                ptrg[2,0] = 3
                ptrg[2,1] = Nk[z[0]]

            ##Periodic boundary conditions in k direction
            if not isAutoPeriodic:
                # Creatioon connectivite perio dans t
                for idir in ['_kmax','_kmin']:
                    if idir == '_kmax':
                        ktg  = zdim[2,0]
                        ktgD = 1
                        if extrusion =='cart':
                            angle = 0
                            trans = -perio
                        else:
                            angle = -perio
                            trans = 0.
                    else:
                        ktgD  = zdim[2,0]
                        ktg   = 1
                        if extrusion =='cart':
                            angle = 0
                            trans = perio
                        else:
                            angle = perio
                            trans = 0.

                    Conn = Internal.getNodeFromName(z, "ZoneGridConnectivity")
                    name = 'match_'+z[0]+idir
                    Internal.createUniqueChild(Conn, name, 'GridConnectivity1to1_t')
                    tmp1     = Internal.getNodeFromName(Conn, name)
                    Internal.setValue(tmp1, z[0])
                    datap = numpy.empty( (3,2) , numpy.int32)
                    datap[0,0]=1+nghost;datap[1,0]=1+nghost;datap[2,0]=ktg
                    datap[0,1]=zdim[0,0]-nghost;datap[1,1]=zdim[1,0]-nghost;datap[2,1]= ktg
                    Internal.createUniqueChild(tmp1 ,'PointRange', 'IndexRange_t',datap)
                    datap = numpy.empty( (3,2) , numpy.int32)
                    datap[0,0]=1+nghost;datap[1,0]=1+nghost;datap[2,0]=ktgD
                    datap[0,1]=zdim[0,0]-nghost;datap[1,1]=zdim[1,0]-nghost;datap[2,1]= ktgD
                    Internal.createUniqueChild(tmp1 ,'PointRangeDonor', 'IndexRange_t',datap)
                    datap = numpy.empty( 3 , numpy.int32)
                    datap[0]=1;datap[1]=2;datap[2]=3
                    Internal.createUniqueChild(tmp1 ,'Transform', '"int[IndexDimension]"',datap)
                    Internal.createUniqueChild(tmp1 ,'GridConnectivityProperty', 'GridConnectivityProperty_t')

                    prop     = Internal.getNodeFromName(tmp1, 'GridConnectivityProperty')
                    Internal.createUniqueChild(prop ,'Periodic', 'Periodic_t')
                    period    = Internal.getNodeFromName(prop, 'Periodic')
                    datap = numpy.zeros( 3 , numpy.float64)
                    Internal.createUniqueChild(period ,'RotationCenter', 'DataArray_t',datap)
                    datap = numpy.zeros( 3 , numpy.float64)
                    datap[0]= angle
                    Internal.createUniqueChild(period ,'RotationAngle' ,'DataArray_t',datap)
                    datap = numpy.zeros( 3 , numpy.float64)
                    datap[2]= trans
                    Internal.createUniqueChild(period ,'Translation', 'DataArray_t',datap)

                    rot     = Internal.getNodeFromName(period, 'RotationAngle')
                    Units=['Kilogram','Meter','Second','Kelvin','Radian']
                    Internal.createUniqueChild(rot ,'DimensionalUnits' ,'DimensionalUnits_t',Units)

        if isAutoPeriodic:
            for node in Internal.getNodesFromName(t,'EquationDimension'): Internal.setValue(node,3)
            for z in Internal.getZones(t):
                C._addBC2Zone(z, z[0]+'periodKmin', 'BCautoperiod', 'kmin')
                C._addBC2Zone(z, z[0]+'periodKmax', 'BCautoperiod', 'kmax')
            BCs = Internal.getNodesFromType(t, "BC_t")
            for bc in BCs:
                if Internal.getValue(bc)=='BCautoperiod':
                    ptrg = Internal.getNodeFromName(bc, "PointRange")[1]
                    ptrg[0,0] = 3
                    ptrg[0,1] = ptrg[0,1]-2

                    ptrg[1,0] = 3
                    ptrg[1,1] = ptrg[1,1]-2

    if extrusion == 'cyl':
        T._cart2Cyl(t, (0,0,0),(1,0,0))
        T._cart2Cyl(tb, (0,0,0),(1,0,0))
    X_IBM._redispatch__(t=t)
    return t, tb

def checkCartesian(t, nghost=0):
    """Checks if the provided mesh in Cartesian or Rectilinear.
     Usage: checkCartesian(t, nghost)"""
    dimPb = Internal.getNodeFromName(t, 'EquationDimension')
    if dimPb is None: raise ValueError('checkCartesian: EquationDimension is missing in input tree.')
    dimPb = Internal.getValue(dimPb)

    nghostZ = 0
    if dimPb == 3: nghostZ = nghost

    isCartesian = 1
    for z in Internal.getZones(t):
        ##X Direction
        coord = Internal.getNodeFromName2(z, 'CoordinateX')[1]
        i = 0
        dx = coord[i+1+nghost, 0+nghost, 0+nghostZ] - coord[i+nghost, 0+nghost, 0+nghostZ]
        for i in range(1,4):
            dx2  = coord[i+1+nghost, 0+nghost, 0+nghostZ] - coord[i+nghost, 0+nghost, 0+nghostZ]
            diff = abs(dx2-dx)
            if diff > 1e-12: isCartesian=0
        if isCartesian < 1: break

        ##Y Direction
        coord = Internal.getNodeFromName2(z, 'CoordinateY')[1]
        i = 0
        dx = coord[0+nghost, i+1+nghost, 0+nghostZ] - coord[0+nghost, i+nghost, 0+nghostZ]
        for i in range(1,4):
            dx2  = coord[0+nghost, i+1+nghost, 0+nghostZ] - coord[0+nghost, i+nghost, 0+nghostZ]
            diff = abs(dx2-dx)
            if diff > 1e-12: isCartesian=0
        if isCartesian < 1: break

        ##Z Direction
        if dimPb==3:
            coord = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
            i     = 0
            dx    = coord[0+nghost,0+nghost,i+1+nghost] - coord[0+nghost,0+nghost,i  +nghost]
            for i in range(1,4):
                dx2  = coord[0+nghost,0+nghost,i+1+nghost] - coord[0+nghost,0+nghost,i  +nghost]
                diff = abs(dx2-dx)
                if diff > 1e-12: isCartesian=0
        if isCartesian < 1: break

    isCartesian = Cmpi.allreduce(isCartesian, op=Cmpi.MIN)
    if isCartesian == 1: cartesian=True
    else:                cartesian=False
    return cartesian

##Files saved here - are not current supported for modified dist2walls but can plug-in if necessary
def _dist2wallNearBody__(t, tb, type='ortho', signed=0, dim=3, loc='centers', model='NSLaminar'):
    import Dist2Walls.PyTree as DTW
    if model == 'NSLaminar':
        list_final_zones=[]
        for z in Internal.getZones(t):
            list_final_zones.append(z[0])

        tBB =G.BB(t)
        tbBB=G.BB(tb)

        ##PRT1 - Zones flagged by the intersection of the bounding boxes of t and tb
        interDict = X.getIntersectingDomains(tBB, tbBB, taabb=tBB, taabb2=tbBB)
        zt       = []
        zt_names = []
        for i in interDict:
            if interDict[i]:
                zt.append(Internal.getNodeByName(t,i))
                zt_names.append(i)

        if zt_names:
            DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)
        del zt

        ###PRT2 - Zones flagged by the intersection of the bounding boxes of t and scaled version of tb
        list_additional_zones = getZonesScaleUpDown__(tbBB,tBB,zt_names,dim=dim)

        ##PRT2
        if list_additional_zones:
            zt=[]
            for i in list_additional_zones:
                zt.append(Internal.getNodeByName(t,i))
            DTW._distance2Walls(zt, tb, type=type, signed=signed, dim=dim, loc=loc)

    else:
        DTW._distance2Walls(t, tb, type=type, signed=signed, dim=dim, loc=loc)

    return None


def getZonesScaleUpDown__(tbBB, tBB, zt_names, diff_percent=0.15, sweep_num=4, scaleDirection=0, dim=2):
    mean_tb=[]
    for bb in Internal.getZones(tbBB): mean_tb.append(getMean__(bb))

    diff_percentz = diff_percent
    if dim == 2: diff_percentz=0

    list_additional_zones=[]
    list_additional_zonesCGNSFile=[]
    for i in range(1, sweep_num+1):
        if scaleDirection >= 0:
            tbBB_scale = T.scale(tbBB, factor=(1.0+i*diff_percent,1.0+i*diff_percent,1.0+i*diff_percentz))
            _add2listAdditionalZones__(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)
        if scaleDirection<=0:
            tbBB_scale = T.scale(tbBB, factor=(1.0-i*diff_percent,1.0-i*diff_percent,1.0-i*diff_percentz))
            _add2listAdditionalZones__(list_additional_zones,tbBB_scale,tBB,mean_tb,zt_names)

    return list_additional_zones


def getMean__(bb):
    mean_local=[]
    for i in ['CoordinateX', 'CoordinateY','CoordinateZ']: mean_local.append(C.getMeanValue(bb, i))
    return mean_local


def _add2listAdditionalZones__(list_additional_zones, tbBB_scale, tBB, mean_tb, zt_names):
    count=0
    for bb in Internal.getZones(tbBB_scale):
        mean_tbscale=getMean__(bb)
        T._translate(bb, (mean_tb[count][0]-mean_tbscale[0],mean_tb[count][1]-mean_tbscale[1],mean_tb[count][2]-mean_tbscale[2]))

        interDict_scale = X.getIntersectingDomains(tBB, bb)
        for i in interDict_scale:
            if interDict_scale[i] and i not in list_additional_zones and i not in zt_names:
                list_additional_zones.append(i)
        count += 1
    return None

###
## Split les zones generere par generate Cartmesh pour empecher qu'une face possede 3 types de raccord: fin/fin, fin/grossier et grossier/fin
## Sinon souci pour flux conservatif
def splitZoneForConservative(t, dimPb=3):
    zone_out=[]
    candidat={}
    for z in Internal.getZones(t):
        matches = Internal.getNodesByType2(z, "GridConnectivity1to1_t")
        nearmatches = Internal.getNodesByType2(z, "GridConnectivity_t")
        totalRac = len(matches)+len(nearmatches)
        #recherche des zones candidate au split
        if totalRac > dimPb*2 and len(nearmatches) != 0:
            #print("nb raccord", totalRac)
            #recherche des faces problematique
            imin=[];jmin=[];kmin=[];imax=[];jmax=[];kmax=[];
            for rac in matches:
                rg = Internal.getNodeByName(rac,'PointRange')[1]
                #CL I
                if rg[0][1] - rg[0][0]==0:
                    if rg[0][1]==1: imin.append([rac, 1])
                    else: imax.append( [rac,1])
                #CL J
                if rg[1][1] - rg[1][0]==0:
                    if rg[1][1]==1: jmin.append([rac,1])
                    else:  jmax.append( [rac,1])
                #CL K
                if dimPb==3:
                    if rg[2][1] - rg[2][0]==0:
                        if rg[2][1]==1: kmin.append( [rac,1])
                        else:  kmax.append( [rac,1])

            for rac in nearmatches:
                rg    = Internal.getNodeByName1(rac,'PointRange')[1]
                ratio = Internal.getNodeByName2(rac,'NMRatio')[1]
                #CL I
                if rg[0][1] - rg[0][0]==0:
                    if rg[0][1]==1: imin.append([rac, ratio[1]])
                    else: imax.append( [rac, ratio[1] ])
                #CL J
                if rg[1][1] - rg[1][0]==0:
                    if rg[1][1]==1: jmin.append([rac, ratio[0]])
                    else:  jmax.append( [rac, ratio[0]])
                #CL K
                if dimPb==3:
                    if rg[2][1] - rg[2][0]==0:
                        if rg[2][1]==1: kmin.append( [rac, ratio[0]])
                        else:  kmax.append( [rac, ratio[0]])

            ladd = True
            for dir in [[imin,'imin'], [imax,'imax'], [jmin,'jmin'],[jmax,'jmax'], [kmin,'kmin'],[kmax,'kmax']]:
                ratio=[]
                for sf in dir[0]:
                    if sf[1] not in ratio: ratio.append(sf[1])

                if len(ratio) >=3:
                    dims = Internal.getZoneDim(z)
                    #print(z[0], ': zone a decouper en', dir[1],  'ratio:', ratio, 'dim:', dims[1:4])
                    for sf in dir[0]:
                        rg    = Internal.getNodeByName1(sf[0],'PointRange')[1]
                        if sf[1]==0.5:
                            print("zone a splitter", z[0], "dir", dir[1], 'ptrange', rg)
                            ladd=False
                            if z[0] in candidat:
                                tmp = candidat[z[0]]
                                tmp.append( [dir[1], rg])
                                candidat[z[0]]= tmp
                            else:
                                candidat[z[0]]= [[dir[1], rg]]
            if ladd:
                zone_out.append(z)
        else:
            zone_out.append(z)

    for key in candidat:
        zone = Internal.getNodeFromName(t,key)
        print("zone candidate au split", zone[0])
        llist = candidat[key]
        for tmp in llist:
            rg2=tmp[1]
            #print("verif init", tmp[0], rg2[0,0], rg2[0,1], rg2[1,0], rg2[1,1], rg2[2,0], rg2[2,1] )
        racFinal=[]
        while len(llist)!=0:
            l = llist[0]
            rg = l[1]

            if len(llist)==1:
                ldel = False
                for rac in racFinal:
                    tmp = rac[1]
                    if rg[1,1] == tmp[1,1] and rg[2,0] == tmp[2,0] and rg[1,0] == tmp[1,0] and  rg[2,1] == tmp[2,1] and  rg[0,1] == tmp[0,1] and rg[0,0] == tmp[0,0]:
                        #print("on skippe: le dernier rac existe deja")
                        ldel=True
                if not ldel:
                    racFinal.append(l)
                llist.pop(0)

            else:
                i =1
                fusion=False
                same = False
                for l2 in llist[1:]:
                    fusion=False
                    rg2= l2[1]
                    if (l[0]=='imax' or l[0]=='imin') and (l2[0]=='imax' or l2[0]=='imin') :
                        if rg[1,0] == rg2[1,1] and (rg[2,0] == rg2[2,0] and rg[2,1] == rg2[2,1]):
                            fusion=True
                            rg[1,0]= rg2[1,0]
                        elif rg[1,1] == rg2[1,0] and (rg[2,0] == rg2[2,0] and rg[2,1] == rg2[2,1]):
                            fusion=True
                            rg[1,1]= rg2[1,1]
                        elif rg[2,0] == rg2[2,1] and (rg[1,0] == rg2[1,0] and rg[1,1] == rg2[1,1]):
                            fusion=True
                            rg[2,0]= rg2[2,0]
                        elif rg[2,1] == rg2[2,0] and (rg[1,0] == rg2[1,0] and rg[1,1] == rg2[1,1]):
                            fusion=True
                            rg[2,1]= rg2[2,1]
                        elif rg[1,1] == rg2[1,1] and rg[2,0] == rg2[2,0] and rg[1,0] == rg2[1,0] and  rg[2,1] == rg2[2,1] :
                            fusion=True
                            same=True

                    elif (l[0]=='jmax' or l[0]=='jmin') and (l2[0]=='jmax' or l2[0]=='jmin') :
                        if rg[0,0] == rg2[0,1] and (rg[2,0] == rg2[2,0] and rg[2,1] == rg2[2,1]):
                            fusion=True
                            rg[0,0]= rg2[0,0]
                        elif rg[0,1] == rg2[0,0] and (rg[2,0] == rg2[2,0] and rg[2,1] == rg2[2,1]):
                            fusion=True
                            rg[0,1]= rg2[0,1]
                        elif rg[2,0] == rg2[2,1] and (rg[0,0] == rg2[0,0] and rg[0,1] == rg2[0,1]):
                            fusion=True
                            rg[2,0]= rg2[2,0]
                        elif rg[2,1] == rg2[2,0] and (rg[0,0] == rg2[0,0] and rg[0,1] == rg2[0,1]):
                            fusion=True
                            rg[2,1]= rg2[2,1]
                        elif rg[0,1] == rg2[0,1] and rg[2,0] == rg2[2,0] and rg[0,0] == rg2[0,0] and  rg[2,1] == rg2[2,1]:
                            fusion=True
                            same=True

                    elif (l[0]=='kmax' or l[0]=='kmin') and (l2[0]=='kmax' or l2[0]=='kmin') :
                        if rg[0,0] == rg2[0,1] and (rg[1,0] == rg2[1,0] and rg[1,1] == rg2[1,1]):
                            fusion=True
                            rg[0,0]= rg2[0,0]
                        elif rg[0,1] == rg2[0,0] and (rg[1,0] == rg2[1,0] and rg[1,1] == rg2[1,1]):
                            fusion=True
                            rg[0,1]= rg2[0,1]
                        elif rg[1,0] == rg2[1,1] and (rg[0,0] == rg2[0,0] and rg[0,1] == rg2[0,1]):
                            fusion=True
                            rg[1,0]= rg2[1,0]
                        elif rg[1,1] == rg2[1,0] and (rg[0,0] == rg2[0,0] and rg[0,1] == rg2[0,1]):
                            fusion=True
                            rg[1,1]= rg2[1,1]
                        elif rg[0,1] == rg2[0,1] and rg[0,0] == rg2[0,0] and rg[1,0] == rg2[1,0] and  rg[1,1] == rg2[1,1]:
                            same=True
                            fusion=True
                    else:
                        print("skip")

                    if fusion:
                        print("dirF", l[0], rg[0,0], rg[0,1], rg[1,0], rg[1,1], rg[2,0], rg[2,1] )
                        llist.pop(i)

                    i +=1

                if not same: racFinal.append(l)
                if not fusion:
                    llist.pop(0)
                    for tmp in llist:
                        rg2=tmp[1]
                    for tmp in racFinal:
                        rg2=tmp[1]
                else:
                    for tmp in llist:
                        rg2=tmp[1]
                    for tmp in racFinal:
                        rg2=tmp[1]

        dims = Internal.getZoneDim(zone)
        for tmp in racFinal:
            rg2=tmp[1]
            #print("verif Final", tmp[0], rg2[0,0], rg2[0,1], rg2[1,0], rg2[1,1], rg2[2,0], rg2[2,1] )
            #passage surface volume
            #print(zone[0], 'dim:', dims[1:4])
            if tmp[0]=='imin': rg2[0,1]= dims[1]
            if tmp[0]=='imax': rg2[0,0]= 1
            if tmp[0]=='jmin': rg2[1,1]= dims[2]
            if tmp[0]=='jmax': rg2[1,0]= 1
            if tmp[0]=='kmin': rg2[2,1]= dims[3]
            if tmp[0]=='kmax': rg2[2,0]= 1
            #print("verif Final1", tmp[0], rg2[0,0], rg2[0,1], rg2[1,0], rg2[1,1], rg2[2,0], rg2[2,1] )

        i=0
        for rac in racFinal:
            rg = rac[1]
            j=0
            for rac2 in racFinal:
                tmp = rac2[1]
                if rg[1,1] == tmp[1,1] and rg[2,0] == tmp[2,0] and rg[1,0] == tmp[1,0] and  rg[2,1] == tmp[2,1] and  rg[0,1] == tmp[0,1] and rg[0,0] == tmp[0,0] and j!=i:
                    #print("on vire le rac: doublon ")
                    racFinal.pop(i)
                    break
                j+=1
            i+=1

        for tmp in racFinal:
            rg2=tmp[1]
            print("Final filtre", zone[0], tmp[0], rg2[0,0], rg2[0,1], rg2[1,0], rg2[1,1], rg2[2,0], rg2[2,1] )

        ni   = dims[1]; nj = dims[2]; nk = dims[3]
        out=[ [zone, (1,1,1), (ni,nj,nk) ]]
        count=0
        for rac in  racFinal:
            rg = rac[1]; ideb = rg[0,0]; jdeb= rg[1,0]; kdeb= rg[2,0]; ifin = rg[0,1]; jfin= rg[1,1]; kfin= rg[2,1]
            start= (ideb, jdeb, kdeb)
            end  = (ifin, jfin, kfin)
            #print("Split", zone[0], tmp[0], rg[0,0], rg[0,1], rg[1,0], rg[1,1], rg[2,0], rg[2,1] )
            #recherche de la zone de decouper
            i=0
            for sub in out:
                start_tg= sub[1]; end_tg = sub[2]
                if start[0] >= start_tg[0] and end[0] <= end_tg[0] and start[1] >= start_tg[1] and end[1] <= end_tg[1] and start[2] >= start_tg[2] and end[2] <= end_tg[2]:
                    zone_tg = sub[0]
                    i_tg=i
                i+=1

            #print('i_tg', i_tg, 'indice_tg', start_tg[0], end_tg[0], start_tg[1],  end_tg[1] ,start_tg[2], end_tg[2])
            dims = Internal.getZoneDim(zone_tg)
            ni   = dims[1]; nj = dims[2]; nk = dims[3]
            ideb = start[0]-start_tg[0]+1; jdeb = start[1]-start_tg[1]+1; kdeb = start[2]-start_tg[2]+1
            ifin = end[0]-start_tg[0]+1; jfin = end[1]-start_tg[1]+1; kfin = end[2]-start_tg[2]+1
            start_loc=( ideb, jdeb, kdeb)
            end_loc=  ( ifin, jfin, kfin)
            #print("start", start)
            #print("end  ",  end)
            #print("start_loc", start_loc)
            #print("end_loc  ",  end_loc)
            #print("dim",  dims[1:4], zone_tg[0])

            #print(' split centre', start_loc, end_loc ,  Internal.getZoneDim(zone_tg)[1:4], flush=True)
            zp  = T.subzone(zone_tg, start_loc, end_loc )
            #print(' split centre OK', flush=True)
            Internal._rmNodesByName(zp,'ZoneGridConnectivity')
            zp[0]=zone[0]+'split'+str(count)
            count +=1
            out.append( [zp, start_loc, end_loc] )
            #decoupe de 6 sous zone potentiel
            if ideb !=1:
                start= (1   ,  1, kdeb)
                end  = (ideb, nj, kfin)
                start_loc=( 1   , 1   , kdeb)
                end_loc=  ( ideb, nj, kfin)

                zp  = T.subzone(zone_tg, start_loc, end_loc )
                zp[0]=zone[0]+'split'+str(count)
                count +=1
                out.append( [zp, start, end] )
            if jfin != nj:
                start= (ideb, jfin, kdeb)
                end  = (ifin, nj  , kfin)
                start_loc=( ideb, jfin, kdeb)
                end_loc=  ( ifin, nj, kfin)
                zp  = T.subzone(zone_tg, start_loc, end_loc )
                zp[0]=zone[0]+'split'+str(count)
                count +=1
                out.append( [zp, start, end] )
            if ifin != ni:
                start= (ifin,  1  , kdeb)
                end  = ( ni , nj  , kfin)
                start_loc=( ifin,   1 , kdeb)
                end_loc=  (  ni ,  nj , kfin)
                zp  = T.subzone(zone_tg, start_loc, end_loc )
                zp[0]=zone[0]+'split'+str(count)
                count +=1
                out.append( [zp, start, end] )
            if jdeb != 1:
                start= (ideb,  1   , kdeb)
                end  = (ifin, jdeb , kfin)
                start_loc=( ideb,  1  , kdeb)
                end_loc=  ( ifin, jdeb, kfin)
                zp  = T.subzone(zone_tg, start_loc, end_loc )
                out.append( [zp, start, end] )
            if kdeb != 1:
                start= (1,  1  , 1)
                end  = (ni, nj , kdeb)
                start_loc=(   1 ,  1  ,  1  )
                end_loc=  (  ni ,  nj , kdeb)
                zp  = T.subzone(zone_tg, start_loc, end_loc )
                zp[0]=zone[0]+'split'+str(count)
                count +=1
                out.append( [zp, start, end] )
            if kfin != nk:
                start= (1,  1  , kfin)
                end  = (ni, nj , nk)
                start_loc=(  1  ,  1  , kfin)
                end_loc=  (  ni ,  nj ,  nk )
                zp  = T.subzone(zone_tg, start_loc, end_loc )
                zp[0]=zone[0]+'split'+str(count)
                count +=1
                out.append( [zp, start, end] )

            out.pop(i_tg)

        for z in out:
            zone_out.append(z[0])

    t = C.newPyTree(['CARTESIAN', zone_out])
    Internal._rmNodesByName(t, 'ZoneGridConnectivity')

    return t
