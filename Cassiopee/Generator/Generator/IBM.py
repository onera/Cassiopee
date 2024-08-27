"""Mesh generation for IBM"""
from . import Generator
from . import generator
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

def generateIBMMesh(tb, vmin=15, snears=0.01, dfar=10., dfarList=[], DEPTH=2, tbox=None,
                    snearsf=None, check=False, sizeMax=4000000,
                    externalBCType='BCFarfield', to=None,
                    fileo=None, expand=2, dfarDir=0, mode=0):
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
                        expand=expand, dfarDir=dfarDir, mode=mode)
    else:
        o = Internal.getZones(to)[0]
    if check: C.convertPyTree2File(o, "octree.cgns")

    # retourne les 4 quarts (en 2D) de l'octree parent 2 niveaux plus haut
    # et les 8 octants en 3D sous forme de listes de zones non structurees
    parento = buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                                   dimPb=dimPb, vmin=vmin, fileout=None, rank=0, dfarDir=dfarDir, mode=mode)
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
#   mode (0 or 1): octree generation mode. If 0: dfar is exact and snear varies. If 1 it's the opposite
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

# only in generateIBMMeshPara and generateCartMesh__
def octree2StructLoc__(o, parento=None, vmin=15, ext=0, optimized=0, sizeMax=4e6, tbOneOver=None):
    sizeMax=int(sizeMax)
    dim = Internal.getZoneDim(o)
    if dim[3] == 'QUAD': dimPb = 2
    elif dim[3] == 'HEXA': dimPb = 3

    if ext == 1: ext = 2
    a = C.getFields(Internal.__GridCoordinates__, o)[0]
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

        ##RECTILINEAR REGION
        tbbB            = G.BB(tbOneOver)
        NBases          = len(Internal.getBases(tbbB))
        interDict_scale = X.getIntersectingDomains(tbbB, tzones)
        for i in interDict_scale:
            listSavetbOneOverTmp = []
            for z in interDict_scale[i]: listSavetbOneOverTmp.append(z)
            listSavetbOneOver.append(listSavetbOneOverTmp)

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
            ## Rectilinear mesh modifications
            if tbOneOver:
                ZONEStbOneOver    = [[],[],[],[]];
                ZONEStbOneOverTmp = []
                for i in range(NBases): ZONEStbOneOverTmp.append([[],[],[],[]])

        else:
            ZONES = [[],[],[],[],[],[],[],[]]; noct = 8 ;
            ## Rectilinear mesh modifications
            if tbOneOver:
                ZONEStbOneOver    = [[],[],[],[],[],[],[],[]];
                ZONEStbOneOverTmp = []
                for i in range(NBases): ZONEStbOneOverTmp.append([[],[],[],[],[],[],[],[]])
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
                        print('Nb of merged zones - OneOverRegion # %d : %d.' %(1,len(ZONEStbOneOverTmp[i][noo])))
                    ZONEStbOneOver[noo].append(ZONEStbOneOverTmp[i][noo])
                ZONEStbOneOver[noo]=sum(ZONEStbOneOver[noo],[])
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
            if ZONEStbOneOver is not None:
                for i in range(NBases):
                    ZONEStbOneOverTmp[i] = T.mergeCart(ZONEStbOneOverTmp[i][0]+ZONEStbOneOverTmp[i][1]+ZONEStbOneOverTmp[i][2]+ \
                                                       ZONEStbOneOverTmp[i][3]+ZONEStbOneOverTmp[i][4]+ZONEStbOneOverTmp[i][5]+ \
                                                       ZONEStbOneOverTmp[i][6]+ZONEStbOneOverTmp[i][7], sizeMax=sizeMax)
                    zones +=ZONEStbOneOverTmp[i]
            if ZONEStbOneOver is not None: del ZONEStbOneOverTmp
            del ZONEStbOneOver
            
        else: # dim=2
            ZONES[0] = T.mergeCart(ZONES[0]+ZONES[2],sizeMax=sizeMax)# XM
            ZONES[1] = T.mergeCart(ZONES[1]+ZONES[3],sizeMax=sizeMax)# XP
            ZONES=ZONES[0:2]
            zones = T.mergeCart(ZONES[0]+ZONES[1],sizeMax=sizeMax)
            del ZONES
            
            ## Rectilinear mesh modifications
            if ZONEStbOneOver is not None:
                for i in range(NBases):
                    ZONEStbOneOverTmp[i] = T.mergeCart(ZONEStbOneOverTmp[i][0]+ZONEStbOneOverTmp[i][1]+ \
                                                       ZONEStbOneOverTmp[i][2]+ZONEStbOneOverTmp[i][3], sizeMax=sizeMax)
                    zones +=ZONEStbOneOverTmp[i]
            del ZONEStbOneOver
            del ZONEStbOneOverTmp
                      
    print('After merging: nb Cartesian zones=%d (ext. =%d).'%(len(zones),ext))
    # Cas ext=-1, ne fait pas les extensions ni les BCs ou raccords
    if ext == -1: return zones
    
    if ext > 0:
        coords = C.getFields(Internal.__GridCoordinates__, zones,api=2)
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

# only in generateIBMMeshPara and generateIBMMesh
def buildParentOctrees__(o, tb, dimPb=3, vmin=15, snears=0.01, snearFactor=1., dfars=10., dfarDir=0, 
                         tbox=None, snearsf=None, to=None, mode=0):
    nzones0 = Internal.getZoneDim(o)[2]
    if nzones0 < 1000: return None

    parento = buildOctree(tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=snearFactor, dfars=dfars, dfarDir=dfarDir, 
                          tbox=tbox, snearsf=snearsf, to=to, expand=-1, balancing=0, mode=mode)

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

# main function
def generateIBMMeshPara(tb, dimPb=3, vmin=15, snears=0.01, dfars=10., dfarDir=0, 
                        tbox=None, snearsf=None, check=False, to=None,
                        ext=2, expand=3, mode=0,
                        tbOneOver=None):
    import KCore.test as test
        # refinementSurfFile: surface meshes describing refinement zones
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
                
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
                        tbox=tbox, snearsf=snearsf, to=to, expand=expand, mode=mode)

    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = buildParentOctrees__(o, tb, dimPb=dimPb, vmin=vmin, snears=snears, snearFactor=4., dfars=dfars, dfarDir=dfarDir, 
                                   tbox=tbox, snearsf=snearsf, to=to, mode=mode)
    
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

        if mode==1:
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

    if Cmpi.rank==0 and check: C.convertPyTree2File(o, 'octree.cgns')

    # Split octree
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[Cmpi.rank]
    del o

    # fill vmin + merge in parallel
    tbOneOverLocal=None
    if tbOneOver:
        tbOneOverLocal=Internal.getNodesFromNameAndType(tbOneOver, '*OneOver*', 'CGNSBase_t')
    res = octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000, tbOneOver=tbOneOverLocal)
    del p
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%Cmpi.rank
    Cmpi._setProc(t, Cmpi.rank)

    C._addState(t, 'EquationDimension', dimPb)

    # Keep F1 regions - for F1 & F42 synergy
    if tbOneOver:
        tbF1 = Internal.getNodesFromNameAndType(tbOneOver, '*KeepF1*', 'CGNSBase_t')
        if tbF1:
            tbbBTmp         = G.BB(tbF1)
            interDict_scale = X.getIntersectingDomains(tbbBTmp, t)
            for kk in interDict_scale:
                for kkk in interDict_scale[kk]:
                    z=Internal.getNodeFromName(t, kkk)
                    Internal._createUniqueChild(z, '.Solver#defineTMP', 'UserDefinedData_t')
                    Internal._createUniqueChild(Internal.getNodeFromName1(z, '.Solver#defineTMP'), 'SaveF1', 'DataArray_t', value=1)
                    node=Internal.getNodeFromName(t, kkk)

    # Add xzones for ext
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)

    # Turn Cartesian grid into a rectilinear grid
    test.printMem(">>> cart grids --> rectilinear grids [start]")        
    listDone = []
    if tbOneOverLocal:
        tbb = G.BB(t)

        if dimPb==2:
            T._addkplane(tbb)
            T._contract(tbb, (0,0,0), (1,0,0), (0,1,0), 0.01)
        
        ## RECTILINEAR REGION
        ## Select regions that need to be coarsened
        tbbB            = G.BB(tbOneOverLocal)
        interDict_scale = X.getIntersectingDomains(tbbB, tbb)
        ## Avoid a zone to be coarsened twice
        for i in interDict_scale:
            (b,btmp) = Internal.getParentOfNode(tbOneOverLocal,Internal.getNodeByName(tbOneOverLocal,i))
            checkOneOver = Internal.getNodeByName(b,".Solver#define") ##Needed for F1 & F42 approach
            if checkOneOver:
                b        = Internal.getNodeByName(b,".Solver#define")
                oneoverX = int(Internal.getNodeByName(b, 'dirx')[1])
                oneoverY = int(Internal.getNodeByName(b, 'diry')[1])
                oneoverZ = int(Internal.getNodeByName(b, 'dirz')[1])
                for z in interDict_scale[i]:
                    if z not in listDone:
                        zLocal = Internal.getNodeFromName(t,z)
                        T._oneovern(zLocal, (oneoverX,oneoverY,oneoverZ));
                        listDone.append(z)
    test.printMem(">>> cart grids --> rectilinear grids [end]")     

    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    if symmetry==0: extBnd = 0
    else: extBnd = ext-1 # nb de ghost cells = ext-1 
    coords, rinds = Generator.extendCartGrids(coords, ext=ext, optimized=1, extBnd=extBnd)
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
#   mode (0 or 1): octree generation mode. If 0: dfar is exact and snear varies. If 1 it's the opposite
# OUT:
#   t (tree): mesh Tree
#==============================================================================
# only in buildOctree
def addRefinementZones__(o, tb, tbox, snearsf, vmin, dim):
    tbSolid = Internal.rmNodesByName(tb, 'IBCFil*')
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
    volmin0 = 1.*volmin0
    while end == 0:
        # Do not refine inside obstacles
        C._initVars(to, 'centers:cellN', 1.)
        to = X_IBM.blankByIBCBodies(to, tbSolid, 'centers', dim)
        C._initVars(to, '{centers:cellNBody}={centers:cellN}')
        nob = 0
        C._initVars(to, 'centers:indicator', 0.)
        for box in boxes:
            volmin2 = 1.09*(snearsf[nob])**(dim)
            C._initVars(to,'centers:cellN',1.)
            tboxl = C.newPyTree(['BOXLOC']); tboxl[2][1][2] = box
            to = X_IBM.blankByIBCBodies(to, tboxl, 'centers', dim)
            fact = 1.1
            while C.getMinValue(to, 'centers:cellN') == 1 and fact < 10.:
                print("Info: addRefinementZones: tbox too small - increase tbox by fact = %2.1f"%(fact))
                box2 = T.scale(box, fact) 
                tboxl[2][1][2] = box2
                to = X_IBM.blankByIBCBodies(to, tboxl, 'centers', dim)
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
            volminloc = C.getMinValue(to, 'centers:vol')
    return Internal.getNodeFromType2(to, 'Zone_t')

def buildOctree(tb, dimPb=3, vmin=15, snears=0.01, snearFactor=1., dfars=10., dfarDir=0, 
                tbox=None, snearsf=None, to=None, balancing=2, expand=2, mode=0):
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
        o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarListL, balancing=balancing, dfarDir=dfarDir, mode=mode)
        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if dxmin < 0.65*dxmin0:
            snearso = [2.*i for i in snearso]
            o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarListL, balancing=balancing, dfarDir=dfarDir, mode=mode)
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
                octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
                indic = C.getField("centers:indicator", o)[0]
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
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = Converter.initVars(indic, 'indicator', 0.)
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0,0,1)
            indic = Converter.extractVars(indic, ["indicator"])
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

        else: #expand = 2, 3 or 4
            corner = 0
            to = C.newPyTree(['Base',o])
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
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
                octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
                indic = C.getField("centers:indicator", o)[0]
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
                octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
                indic = C.getField("centers:indicator", o)[0]
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
def _projectMeshSize(t, NPas=10, span=1, dictNz=None, isCartesianExtrude=False):
    """Predicts the final size of the mesh when extruding 2D to 3D in the z-direction.
    Usage: loads(t, NPas, span, dictNz, isCartesianExtrude)"""
    NP             = Cmpi.size
    rank           = Cmpi.rank
    NPTS           = numpy.zeros(NP)
    NCELLS         = numpy.zeros(NP)
    NPTS_noghost   = numpy.zeros(NP)
    NCELLS_noghost = numpy.zeros(NP)    
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

