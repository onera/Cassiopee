"""Toolbox for IBM preprocessing"""
import numpy
from . import PyTree as X
from . import OversetData as XOD
from . import Connector
from . import connector

try: range = xrange
except: pass

try:
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Transform.PyTree as T
    import Converter.Internal as Internal
    import Dist2Walls.PyTree as DTW
    import Post.PyTree as P
    import Converter
    import Generator
    import Transform
    import Converter.GhostCells as CGC
    import KCore
    import numpy
    import math
except:
    raise ImportError("Connector.ToolboxIBM requires Converter, Generator, Transform, Dist2Walls and Post modules.")

varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
TOLDIST = 1.e-14
SHIFTF = 1.e-10
EPSCART = 1.e-6
TOLCELLN = 0.01

##from param_solver.h
LBM_IBC_NUM        = 113

TypesOfIBC = XOD.TypesOfIBC

# ?
def _blankClosestTargetCells(t, cellNName='cellN', depth=3):
    for z in Internal.getZones(t):
        connector._blankClosestTargetCells(z, depth, cellNName,
                                           Internal.__GridCoordinates__,
                                           Internal.__FlowSolutionNodes__,
                                           Internal.__FlowSolutionCenters__)
    return None
        
# ==============================================================================
# Generates the fully Cartesian IBM mesh
# ==============================================================================
# Reduction de la taille des fenetres des BC physiques pour qu elles soient 
# traitees comme des ghost cells
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

#---------------------------------------------------------------------------------
# INPUT: t:  
#         tb: bodies - to ensure the refinement prescribed       
#         sensor function to be already computed
#         factor: nb of points is roughly multiplied by factor after remeshing
#----------------------------------------------------------------------------------
def adaptIBMMesh(t, tb, vmin, sensor, factor=1.2, DEPTH=2, sizeMax=4000000,
                 variables=None, refineFinestLevel=False, refineNearBodies=False, 
                 check=True, symmetry=0, externalBCType='BCFarfield', fileo='octree.cgns'):
    import Converter.Mpi as Cmpi
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
    to = P.computeIndicatorValue(to, t, sensor)
    res = P.computeIndicatorField(to, sensor, nbTargetPts=factor*npts, \
                                  bodies=constraintSurfaces, \
                                  refineFinestLevel=refineLevelF, \
                                  coarsenCoarsestLevel=1)
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
                            sizeMax=sizeMax, check=check, symmetry=symmetry, 
                            externalBCType=externalBCType)
    
    #C._initVars(t2,"centers:Density=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))
    #C._initVars(t2,"centers:VelocityX=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityX'))))
    #C._initVars(t2,"centers:VelocityY=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityY'))))
    #C._initVars(t2,"centers:VelocityZ=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityZ'))))
    #C._initVars(t2,"centers:Temperature=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Temperature'))))
    #C._initVars(t2,"centers:TurbulentSANuTilde=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'TurbulentSANuTildeDensity'))/Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))

    # interpolate the solution on the new mesh
    P._extractMesh(t, t2, 3, mode='accurate')
    return t2

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
        elif len(pool)==1: res += pool
    return res

def octree2StructLoc__(o, parento=None, vmin=21, ext=0, optimized=0, sizeMax=4e6):
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

    if parento is None:
        zones = T.mergeCart(zones,sizeMax=sizeMax)
    else:     
        eps=1.e-10
        bbo = G.bbox(parento[0])# 1st octant lower left side
        xmeano=bbo[3]; ymeano=bbo[4]; zmeano=bbo[5]
        # gather zones by parent octant
        if dimPb == 2: ZONES=[[],[],[],[]]; noct = 4
        else: ZONES = [[],[],[],[],[],[],[],[]]; noct = 8
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
            if noo > -1: ZONES[noo].append(z)
        #-----------------------------------------------------------------------------
        zones=[]
        for noo in range(noct):
            nzones = len(ZONES[noo])
            if nzones > 1:
                print('Merging %d Cartesian zones of subdomain %d.'%(nzones,noo))                
                ZONES[noo] = mergeByParent__(ZONES[noo], parento[noo], sizeMax)
                print('Nb of merged zones : %d.' %len(ZONES[noo]))

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

        else: # dim=2
            ZONES[0] = T.mergeCart(ZONES[0]+ZONES[2],sizeMax=sizeMax)# XM
            ZONES[1] = T.mergeCart(ZONES[1]+ZONES[3],sizeMax=sizeMax)# XP
            ZONES=ZONES[0:2]
            zones = T.mergeCart(ZONES[0]+ZONES[1],sizeMax=sizeMax)
            del ZONES
    print('After merging: nb Cartesian zones=%d.'%(len(zones)))

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


def buildParentOctrees__(o, tb, snears=None, snearFactor=4., dfar=10., dfarList=[], to=None, tbox=None, snearsf=None, 
                         dimPb=3, vmin=15, symmetry=0, fileout=None, rank=0, dfarDir=0):
    nzones0 = Internal.getZoneDim(o)[2] 
    if nzones0 < 1000: return None

    parento = buildOctree(tb, snears=snears, snearFactor=snearFactor, dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf, 
                          dimPb=dimPb, vmin=vmin, symmetry=symmetry, balancing=0, rank=rank, expand=-1, dfarDir=dfarDir)
    
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

# IN: bbox: bbox des frontieres exterieures
def generateCartMesh__(o, parento=None, dimPb=3, vmin=11, DEPTH=2, sizeMax=4000000, check=True, 
                       symmetry=0, externalBCType='BCFarfield', bbox=None):

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

def _addBCOverlaps(t, bbox):
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

def _addExternalBCs(t, bbox, DEPTH=2, externalBCType='BCFarfield', dimPb=3):
    dirs = [0,1,2,3,4,5]
    rangeDir=['imin','jmin','kmin','imax','jmax','kmax']
    if dimPb == 2: dirs = [0,1,3,4]
    nptsTot = 0
    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        nptsTot += niz*njz*nkz
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

#--------------------------------------------------------------------------
# to : maillage octree, si not None : on le prend comme squelette 
#--------------------------------------------------------------------------
def buildOctree(tb, snears=None, snearFactor=1., dfar=10., dfarList=[], to=None, tbox=None, snearsf=None, 
                dimPb=3, vmin=15, balancing=2, symmetry=0, fileout=None, rank=0, expand=2, dfarDir=0):
    import Converter.Mpi as Cmpi
    i = 0; surfaces=[]; snearso=[] # pas d'espace sur l'octree
    bodies = Internal.getZones(tb)
    if not isinstance(snears, list): snears = len(bodies)*[snears]
    if len(bodies) != len(snears):
        raise ValueError('buildOctree: Number of bodies is not equal to the size of snears.')
    dxmin0 = 1.e10
    for s in bodies:
        sdd = Internal.getNodeFromName1(s, ".Solver#define")
        if sdd is not None:
            snearl = Internal.getNodeFromName1(sdd, "snear")
            if snearl is not None: 
                snearl = Internal.getValue(snearl)
                snears[i] = snearl*snearFactor
        dhloc = snears[i]*(vmin-1)
        # if C.isNamePresent(s,'centers:cellN') != -1:
        #     s2 = P.selectCells(s,'{centers:cellN}>0.')
        #     surfaces.append(s2)
        # else: surfaces += [s] 
        surfaces += [s]
        snearso += [dhloc]
        dxmin0 = min(dxmin0, dhloc)
        i += 1
    if to is not None:
        o = Internal.getZones(to)[0]
    else:
        o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarList, balancing=balancing,dfarDir=dfarDir)
        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if dxmin < 0.65*dxmin0:
            snearso = [2.*i for i in snearso]
            o = G.octree(surfaces, snearList=snearso, dfar=dfar, dfarList=dfarList, balancing=balancing, dfarDir=dfarDir)
        # Adaptation avant expandLayer (pour corriger eventuellement les sauts de maille)
        if tbox is not None and snearsf is not None:
            o = addRefinementZones(o, tb, tbox, snearsf, vmin, dimPb)
            C._rmVars(o, ['centers:indicator', 'centers:cellN', 'centers:vol', 'centers:cellNBody'])

        #if expand > 0: C.convertPyTree2File(o, 'startOctree.cgns')
        if expand == 0:
            G._expandLayer(o, level=0, corners=1, balancing=1)
        elif expand == 1:
            vmint = 31
            if vmin < vmint:
                if rank==0: print('buildOctree: octree finest level expanded (expandLayer activated).')
                to = C.newPyTree(['Base',o])
                to = blankByIBCBodies(to, tb, 'centers', dimPb)
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
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            indic = C.getField("centers:cellN",to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = Converter.initVars(indic, 'indicator', 0.)
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0,0,1)
            indic = Converter.extractVars(indic, ["indicator"])
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

        elif expand == 2: # expand minimum
            corner = 0
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 3)          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])
            
            # Check
            #to = C.newPyTree(['Base',o])
            #to = blankByIBCBodies(to, tb, 'centers', dimPb)
            #C._initVars(o, "centers:indicator", 0.)
            #cellN = C.getField("centers:cellN", to)[0]
            #octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            #indic = C.getField("centers:indicator", o)[0]
            #indic = Converter.addVars([indic,cellN])
            #indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 5)
            # FIN CHECK
            
        elif expand == 3: # expand minimum + 1 couche propagee
            #C.convertPyTree2File(o, 'octree1.cgns')
            corner = 0
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 3)          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])
            
            # passe 2
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 4)                                                                          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])
            # fin passe 2

            # Check
            #to = C.newPyTree(['Base',o])
            #to = blankByIBCBodies(to, tb, 'centers', dimPb)
            #C._initVars(o, "centers:indicator", 0.)
            #cellN = C.getField("centers:cellN", to)[0]
            #octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            #indic = C.getField("centers:indicator", o)[0]
            #indic = Converter.addVars([indic,cellN])
            #indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 5)
            # FIN CHECK

        elif expand == 4: # expand minimum + 2 couche propagee
            corner = 0
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 3)          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])
            
            # passe 2
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 4)                                                                          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

            # passe 3
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o, "centers:indicator", 0.)
            cellN = C.getField("centers:cellN", to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator", o)[0]
            indic = Converter.addVars([indic,cellN])
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic, 0, corner, 6)                                                                          
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0], [octreeA])

        #if expand > 0: C.convertPyTree2File(o, 'endOctree.cgns')
        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        C._rmVars(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if rank == 0: print('Minimum spacing of Cartesian mesh= %f (targeted %f)'%(dxmin/(vmin-1),dxmin0/(vmin-1)))

        nelts = Internal.getZoneDim(o)[2] 
        if nelts > 20000: 
            print('Warning: number of zones (%d) on rank %d is high (block merging might last a long time).'%(nelts, rank))

    if fileout is not None:
        if Cmpi.size==1: C.convertPyTree2File(o, fileout)
    return o

def generateIBMMesh(tb, vmin=15, snears=None, dfar=10., dfarList=[], DEPTH=2, tbox=None, 
                    snearsf=None, check=True, sizeMax=4000000, 
                    symmetry=0, externalBCType='BCFarfield', to=None, 
                    fileo=None, expand=2, dfarDir=0):
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

    o = buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf, 
                    dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=fileo, rank=0,
                    expand=expand, dfarDir=dfarDir)

    if check: C.convertPyTree2File(o, "octree.cgns")
    
    # retourne les 4 quarts (en 2D) de l'octree parent 2 niveaux plus haut 
    # et les 8 octants en 3D sous forme de listes de zones non structurees
    parento = buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf, 
                                   dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=0, dfarDir=dfarDir)
    res = generateCartMesh__(o, parento=parento, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH, sizeMax=sizeMax, 
                             check=check, symmetry=symmetry, externalBCType=externalBCType)
    return res

# Remove fully blanked grids considering cellNNIBC and cellNChim
def _removeBlankedGrids(t, loc='centers'):
    vari = 'cellNIBC'
    varc = 'cellNChim'
    flag = 'flag'
    if loc == 'centers': 
        vari = 'centers:'+vari
        varc = 'centers:'+varc
        flag = 'centers:'+flag
    C._initVars(t,'{%s}=abs(1.-{%s}*{%s})<0.5'%(flag,vari,varc))

    for z in Internal.getZones(t):
        if C.getMaxValue(z,flag) < 0.5: 
            (parent,noz) = Internal.getParentOfNode(t, z)
            del parent[2][noz]
        else:
            C._rmVars(z,[flag])
    return None

# =============================================================================
# create refinement zones inside tbox bodies with spacing snearsf
# snearsf can be a float or a list of floats. In that case, snears length
# and number of boxes must be equal
# =============================================================================
def addRefinementZones(o, tb, tbox, snearsf, vmin, dim):
    boxes = []
    for b in Internal.getBases(tbox):
        boxes.append(Internal.getNodesFromType1(b, 'Zone_t'))
        
    if snearsf != []:
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
        to = blankByIBCBodies(to, tb, 'centers', dim)
        C._initVars(to, '{centers:cellNBody}={centers:cellN}')
        nob = 0
        C._initVars(to, 'centers:indicator', 0.)
        for box in boxes:
            volmin2 = 1.09*(snearsf[nob])**(dim)
            C._initVars(to,'centers:cellN',1.)
            tboxl = C.newPyTree(['BOXLOC']); tboxl[2][1][2] = box
            to = blankByIBCBodies(to, tboxl, 'centers', dim)
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

# =============================================================================
# Calcul des points IBM a corriger, paroi et a interpoler
# =============================================================================
def getAllIBMPoints(t, loc='nodes', hi=0., he=0., tb=None, tfront=None, 
                    frontType=0, cellNName='cellN', IBCType=1, depth=2, Reynolds=6.e6, yplus=100.,isLBM=False):
    if IBCType == -1: signOfDistCorrected = -1
    else: signOfDistCorrected=1 # signe de la distance aux points corriges

    allCorrectedPts = []; allWallPts = []; allInterpPts = []
    #allInterpPts_layer2 = []; # [AJ] Keep for now
    #-------------------------------------------
    # 1. Get the list of IBC corrected pts
    #-------------------------------------------
    listOfSnearsLoc=[]
    listOfModelisationHeightsLoc = []
    if loc == 'nodes':
        for z in Internal.getZones(t):
            an = C.getFields(Internal.__GridCoordinates__,z)[0]
            ac1 = C.getField(cellNName,z)[0]
            ac1[0] = 'cellN'
            ac2 = C.getField('TurbulentDistance',z)[0]
            ac3 = C.getField('gradxTurbulentDistance',z)[0]
            ac4 = C.getField('gradyTurbulentDistance',z)[0]
            ac5 = C.getField('gradzTurbulentDistance',z)[0]
            an = Converter.addVars([an,ac1,ac2,ac3,ac4,ac5])
            ah = C.getField('hi',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            ah = C.getField('he',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            correctedPts = Connector.getInterpolatedPoints__(an) 
            xt = C.getField('CoordinateX',z)[0][1][0]
            snearl = xt[1]-xt[0]
            listOfSnearsLoc.append(snearl)
            allCorrectedPts.append(correctedPts)
            if frontType == 42: listOfModelisationHeightsLoc.append(computeModelisationHeight(Re=Reynolds, yplus=yplus))
            else: listOfModelisationHeightsLoc.append(0.)
    else:
        for z in Internal.getZones(t):            
            an = C.getFields(Internal.__GridCoordinates__,z)[0]
            an = Converter.node2Center(an)
            ac1 = C.getField('centers:'+cellNName,z)[0]
            ac1[0] = 'cellN'
            ac2 = C.getField('centers:TurbulentDistance',z)[0]
            ac3 = C.getField('centers:gradxTurbulentDistance',z)[0]
            ac4 = C.getField('centers:gradyTurbulentDistance',z)[0]
            ac5 = C.getField('centers:gradzTurbulentDistance',z)[0]
            an = Converter.addVars([an,ac1,ac2,ac3,ac4,ac5])
            ah = C.getField('centers:hi',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            ah = C.getField('centers:he',z)[0]
            if ah != []: an = Converter.addVars([an,ah])
            correctedPts = Connector.getInterpolatedPoints__(an) 
            allCorrectedPts.append(correctedPts)
            xt = C.getField('CoordinateX',z)[0][1][0]
            snearl = xt[1]-xt[0]
            listOfSnearsLoc.append(snearl)
            if frontType == 42: listOfModelisationHeightsLoc.append(computeModelisationHeight(Re=Reynolds, yplus=yplus))
            else: listOfModelisationHeightsLoc.append(0.)
    #-------------------------------------------
    # 2. Get the list of IBC wall and interp pts
    #-------------------------------------------        
    indcell = Converter.extractVars(allCorrectedPts,['indcell'])
    if tb is None or tfront is None: # constant hi, he   
        for nozc in range(len(allCorrectedPts)):
            poshe = KCore.isNamePresent(allCorrectedPts[nozc],'he')
            if poshe == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'he',he)
            poshi = KCore.isNamePresent(allCorrectedPts[nozc],'hi')
            if poshi == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'hi',hi)
        res = connector.getIBMPtsBasic(allCorrectedPts, varsn, 'TurbulentDistance')
    else:
        dictOfBodiesByIBCType={}
        for s in Internal.getZones(tb):
            sdd = Internal.getNodeFromName1(s,".Solver#define")
            if sdd is not None: # check consistency of ibc type with flow equations
                ibctype = Internal.getNodeFromName1(sdd, "ibctype")
                if ibctype is not None:
                    ibctype = Internal.getValue(ibctype)
                else: 
                    if IBCType == -1: ibctype = 'slip' 
                    else: ibctype = 'Musker'
            else: # type of IBC not found: Euler -> slip, other : Musker            
                if IBCType == -1: ibctype = 'slip' 
                else: ibctype = 'Musker'
	    #Over ride as LBM only supports no slip at the moment
            if isLBM==True:
                ibctype = 'noslip'
            famName = Internal.getNodeFromType1(s,'FamilyName_t')
            if famName is not None:
               famName = Internal.getValue(famName)

            ibctypeI = TypesOfIBC[ibctype]
            if famName is not None: ibctype2 = str(ibctypeI)+"#"+famName
            else: ibctype2 = str(ibctypeI)  
            if ibctype2 not in dictOfBodiesByIBCType: dictOfBodiesByIBCType[ibctype2]=[s]
            else: dictOfBodiesByIBCType[ibctype2]+=[s]

        # Regroupement des corps par type de BC - optimise les projections ensuite 
        bodies = []; listOfIBCTypes=[]
        for itype in dictOfBodiesByIBCType:
            s = dictOfBodiesByIBCType.get(itype)
            body = C.getFields(Internal.__GridCoordinates__,s)
            body = Converter.convertArray2Tetra(body)
            body = Transform.join(body)
            bodies.append(body)
            listOfIBCTypes.append(itype)

        if frontType == 0:
            dmin = C.getMaxValue(tfront, 'TurbulentDistance')
            allCorrectedPts = Converter.initVars(allCorrectedPts,'dist',dmin)
            res = connector.getIBMPtsWithoutFront(allCorrectedPts, bodies, varsn, 'dist', signOfDistCorrected)
        else:            
            front = C.getFields(Internal.__GridCoordinates__,tfront)
            front = Converter.convertArray2Tetra(front)
            allCorrectedPts = Converter.extractVars(allCorrectedPts,['CoordinateX','CoordinateY','CoordinateZ']+varsn)
            res = connector.getIBMPtsWithFront(allCorrectedPts, listOfSnearsLoc, listOfModelisationHeightsLoc, bodies, 
                                               front, varsn, signOfDistCorrected, depth)
    allWallPts = res[0]
    allWallPts = Converter.extractVars(allWallPts,['CoordinateX','CoordinateY','CoordinateZ'])

    allInterpPts = res[1] 
    allInterpPts = Converter.extractVars(allInterpPts,['CoordinateX','CoordinateY','CoordinateZ'])
    allInterpPts = Converter.addVars([allInterpPts,indcell])
    allCorrectedPts = Converter.extractVars(allCorrectedPts,['CoordinateX','CoordinateY','CoordinateZ'])

    dictOfInterpPtsByIBCType={} 
    dictOfCorrectedPtsByIBCType={}
    dictOfWallPtsByIBCType={}
    nzonesR = len(allInterpPts)
    if len(res)==3: 
        allIndicesByIBCType = res[2]    
        for noz in range(nzonesR):
            indicesByTypeForZone = res[2][noz]
            nbTypes = len(indicesByTypeForZone)
            for nob in range(nbTypes):
                ibcTypeL = listOfIBCTypes[nob]
                indicesByTypeL = indicesByTypeForZone[nob]
                if indicesByTypeL.shape[0] > 0:                
                    correctedPtsL = Transform.subzone(allCorrectedPts[noz], indicesByTypeL)
                    interpPtsL = Transform.subzone(allInterpPts[noz], indicesByTypeL)
                    wallPtsL = Transform.subzone(allWallPts[noz], indicesByTypeL)                
                else: 
                    correctedPtsL=[]; interpPtsL = []; wallPtsL = []
                if noz == 0:
                    dictOfCorrectedPtsByIBCType[ibcTypeL] = [correctedPtsL]
                    dictOfWallPtsByIBCType[ibcTypeL] = [wallPtsL]
                    dictOfInterpPtsByIBCType[ibcTypeL] = [interpPtsL]
                else:
                    dictOfCorrectedPtsByIBCType[ibcTypeL] += [correctedPtsL]
                    dictOfWallPtsByIBCType[ibcTypeL] += [wallPtsL]
                    dictOfInterpPtsByIBCType[ibcTypeL] += [interpPtsL]
    else:
        if IBCType == -1: ibcTypeL = 0
        else: ibcTypeL = 3 
        for noz in range(nzonesR):
            if noz == 0:
                dictOfCorrectedPtsByIBCType[ibcTypeL] = [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] = [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] = [allInterpPts[noz]]
            else:
                dictOfCorrectedPtsByIBCType[ibcTypeL] += [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] += [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] += [allInterpPts[noz]]        

    return dictOfCorrectedPtsByIBCType, dictOfWallPtsByIBCType, dictOfInterpPtsByIBCType
#=============================================================================
# Returns the front defining the image points
#=============================================================================
def getIBMFront(tc, frontvar, dim, frontType):

    # if frontType == 1 or frontType == 2 : front = getIBMFrontType1(tc,frontvar,dim)
    if frontType > 0 : front = getIBMFrontType1(tc,frontvar,dim)
    else: front = getIBMFrontType0(tc,frontvar,dim)
    front = C.deleteEmptyZones(front)
    Internal._rmNodesFromName(front,"ID_*")
    if frontType == 0: return front
    
    dxmin = 1.e12
    if frontType>0:
        front = Internal.getZones(front)
        dxmax = 1.e-12
        dht = [[]]*len(front)
        nof = 0
        for f in front:
            subf = T.subzone(f, [0], type='elements')
            dx = C.getMaxValue(subf,"CoordinateX")-C.getMinValue(subf,"CoordinateX")
            dy = C.getMaxValue(subf,"CoordinateY")-C.getMinValue(subf,"CoordinateY")
            if dim == 3:
                dz = C.getMaxValue(subf,"CoordinateZ")-C.getMinValue(subf,"CoordinateZ")
                dht[nof] = max(dx,dy,dz)
            else: 
                dht[nof] = max(dx,dy)

            if dht[nof] < dxmin and dht[nof] > 1.e-12:
                dxmin = dht[nof]
            if dht[nof] > dxmax: 
                dxmax = dht[nof]
            nof+=1

        nlevels = int(math.log(dxmax/dxmin)/math.log(2)+1)

        dictOfLevels={}
        for l in range(nlevels): dictOfLevels[l]=[]
        for nof in range(len(dht)):
            if dht[nof]>1.e-12:
                nolvl = int(math.log(dht[nof]/dxmin)/math.log(2))
                dictOfLevels[nolvl].append(front[nof])
        
        front=[]
        for nol in range(nlevels):
            if len(dictOfLevels[nol])>0:
                front.append(T.join(dictOfLevels[nol]))

    return front

# front of first computed cells - with overlapping
def getIBMFrontType1(tc, frontvar, dim):
    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        if z0 is not None:
            zmean = C.getValue(z0, 'CoordinateZ', 0)
            dz = 2*zmean
        else: zmean=0.5; dz = 1.
    else: dz = 0.
    front = []
    for z in Internal.getZones(tc):
        if C.getMinValue(z,frontvar)<0.2 and C.getMaxValue(z,frontvar)>0.8:
            X._maximizeBlankedCells(z,depth=1,dir=1,loc='nodes', cellNName='cellNChim')
            C._initVars(z,'{cellNChim}=minimum(1.,{cellNChim})')
            f = P.frontFaces(z, frontvar)
            if Internal.getZoneDim(f)[1]>0:  
                Internal._rmNodesByName(f,'ID_*')
                front.append(f)
    C._initVars(front,'{tag}=({cellNChim}>0.5)*({cellNChim}<1.5)')
    front = P.selectCells2(front, 'tag', strict=1)
    Internal._rmNodesByName(front,Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(front,Internal.__FlowSolutionCenters__)
    if dim == 2:
        T._addkplane(front)
        T._translate(front,(0,0,-zmean))
        T._contract(front, (0,0,0), (1,0,0), (0,1,0), 0.9*dz)
    return front

# isosurface of max dist of the first interpolable cells
def getIBMFrontType0(tc, frontvar, dim):
    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        zmean = C.getValue(z0, 'CoordinateZ', 0)
        dz = 2*zmean
    else: dz = 0.

    SHIFTD = 1.+SHIFTF
    front = []
    tf = Internal.copyRef(tc)
    C._initVars(tf,'{%s}={%s}-2.*({%s}>1.5)'%(frontvar,frontvar,frontvar))
    for z in Internal.getZones(tf):
        if C.getMinValue(z,frontvar)==0. and C.getMaxValue(z,frontvar)==1.:
            f = P.frontFaces(z, frontvar)
            if Internal.getZoneDim(f)[1]>0:  
                Internal._rmNodesByName(f,'ID_*')
                front.append(f)
    if dim == 2:
        dmin = C.getMaxValue(front, 'TurbulentDistance')
        # Creation du corps 2D pour le preprocessing IBC
        tcl = T.addkplane(tc)
        T._translate(tcl,(0,0,-zmean))
        T._contract(tcl, (0,0,0), (1,0,0), (0,1,0), dz)
        front = P.isoSurfMC(tcl,'TurbulentDistance',dmin*SHIFTD)
        del tcl
    else:
        dmin = C.getMaxValue(front, 'TurbulentDistance')
        front = P.isoSurfMC(tc, 'TurbulentDistance', dmin*SHIFTD)
    return front

#=============================================================================
def getMinimumCartesianSpacing(t):
    baseC = Internal.getNodeFromName1(t, 'CARTESIAN')
    if baseC is None: return -1.

    zonesC = Internal.getZones(baseC)
    dxmin = 1.e6
    for z in zonesC:
        dx = abs(C.getValue(z,'CoordinateX',1)-C.getValue(z,'CoordinateX',0))
        if dx < dxmin: dxmin = dx

    print('Minimum spacing on Cartesian grids = %f.'%dxmin)
    return dxmin

#==============================================================================
# masquage par les corps IBC
# gridType = single or composite - composite means that an off body grid exists
#==============================================================================
def blankByIBCBodies(t, tb, loc, dim, cellNName='cellN'):
    DIM = dim
    blankalgo='tri'
    if DIM == 2: blankalgo = 'xray'

    bodies = []
    for b in Internal.getBases(tb):
        wallsl = Internal.getNodesFromType1(b, 'Zone_t')
        #soldef = Internal.getNodeFromName(wallsl,'.Solver#define')
        bodies.append(wallsl)
        # if wallsl != []:
        #     try:
        #         wallsl = C.convertArray2Tetra(wallsl)
        #         wallsl = T.join(wallsl)
        #         wallsl = G.close(wallsl)
        #         Internal.addChild(wallsl,soldef)
        #         bodies.append([wallsl])
        #         # if DIM == 3:
        #         #     try: P.exteriorFaces(wallsl)
        #         #     except: pass
        #         #     bodies.append([wallsl])
        #         # else: bodies.append([wallsl])
        #     except: 
        #         wallsl = C.convertArray2Tetra(wallsl)
        #         bodies.append(wallsl)

    nbodies = len(bodies)
    if nbodies == 0: 
        print("Warning: blankByIBCBodies: no body defined.")
        return t

    print('Blanking mesh by %d bodies'%nbodies)
    if loc == 'centers': typeb = 'center_in'
    else: typeb = 'node_in'
    nbases = len(Internal.getBases(t))

    bodiesInv=[]
    for body in bodies:
        inv = Internal.getNodeFromName(body,'inv')
        if inv is not None: inv = Internal.getValue(inv)
        else: inv = 0
        if inv==1:
            bodies.remove(body)
            bodiesInv.append(body)

    if blankalgo == 'xray' or DIM == 2:
        BM = numpy.ones((nbases,nbodies),dtype=numpy.int32)
        dh_min = getMinimumCartesianSpacing(t)
        XRAYDIM1 = 2000; XRAYDIM2 = XRAYDIM1
        if dh_min > 0.:
            bb = G.bbox(tb)
            Lxref = bb[3]-bb[0]
            Lyref = bb[4]-bb[1]
            XRAYDIM1 = max(XRAYDIM1,int(Lxref/(0.15*dh_min)))
            XRAYDIM2 = max(XRAYDIM2,int(Lyref/(0.15*dh_min)))
        if DIM == 2: XRAYDIM2 = 2
      
        if loc == 'centers':
            tc = C.node2Center(t)
            for body in bodiesInv:
                print('Reverse blanking for body')
                tc = X.blankCells(tc, [body], BM, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)
                C._initVars(tc,'{%s}=1.-{%s}'%(cellNName,cellNName)) # ecoulement interne

            for body in bodies:
                tc = X.blankCells(tc, [body], BM, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)
            C._cpVars(tc,'%s'%cellNName,t,'centers:%s'%cellNName)
        else:
            t = X.blankCells(t, bodies, BM, blankingType=typeb, delta=TOLDIST, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=DIM, cellNName=cellNName)
    else:
        BM2 = numpy.ones((nbases,1),dtype=numpy.int32)
        for body in bodiesInv:
            print('Reverse blanking for body') 
            t = X.blankCellsTri(t, [body], BM2, blankingType=typeb, cellNName=cellNName)
            C._initVars(t,'{centers:%s}=1-{centers:%s}'%(cellNName,cellNName)) # ecoulement interne
        for body in bodies: 
            t = X.blankCellsTri(t, [body], BM2, blankingType=typeb, cellNName=cellNName)
    return t

#=============================================================================
# distance signee en fonction du masquage Chimere et IBC
#=============================================================================
def _signDistance(t):
    C._initVars(t,'{centers:TurbulentDistance}=-1.*({centers:cellNIBC}*{centers:cellNChim}<1.)*{centers:TurbulentDistance}+({centers:cellNIBC}*{centers:cellNChim}>0.)*{centers:TurbulentDistance}')
    return None

#=============================================================================
# Gather front 
# Si un front est calcule par morceau sur chaque proc, ramene le meme front
# sur tous les procs
#=============================================================================
def gatherFront(front):
    import Converter.Mpi as Cmpi
    zones = Internal.getNodesFromType1(front, 'Zone_t')
    for z in zones: z[0] += '_'+str(Cmpi.rank)
    
    if Cmpi.KCOMM is not None: 
        #allFront = Cmpi.KCOMM.allgather(front)
        #front = []
        #for f in allFront: front += f
        front = Cmpi.allgatherZones(front)
        return front
    else: return front    
    
#=============================================================================
def doInterp(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, 
             frontType=0, depth=2, IBCType=1, interpDataType=1, Reynolds=6.e6, yplus=100.,isLBM=False):    
    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                if interpDataType==1 and dictOfADT is not None:
                                    hook0.append(dictOfADT[zdnrname])
            if interpDataType == 0: hook0 = None

            X._setInterpData(zrcv, dnrZones, nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera')
            for nod in range(len(dnrZones)):
                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]
                tc[2][nobd][2][nozd] = dnrZones[nod]

    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus,isLBM=isLBM)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1] 
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])
                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
                    if interpDataType == 1: hook0 = []
                    else: hook0 = None
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if interpDataType == 1:
                                            if zdnrname not in dictOfADT: 
                                                HOOKADT = C.createHook(zdnr, 'adt')
                                                dictOfADT[zdnrname] = HOOKADT
                                            hook0.append(dictOfADT[zdnrname])

                                        dnrZones.append(zdnr)
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             nature=1, penalty=1, loc='centers', storage='inverse',  
                                             interpDataType=interpDataType, hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL)
                    nozr += 1
                    for nod in range(len(dnrZones)):
                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]
                        tc[2][nobd][2][nozd] = dnrZones[nod]

        if dictOfADT is not None: 
            for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

#=============================================================================
def doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, front=None, 
             frontType=0, depth=2, IBCType=1, interpDataType=1):    
    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension') 
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]
  
        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])
        
        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if (dimPb == 2):dico_dz[z[0]] = 1
        else : dico_dz[z[0]] = dzz

        if (dimPb == 2):dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if (dx > dxmax):dxmax=dx

    niveaux_temps = {} 
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if (dimPb == 2):dzz=max(dxx,dyy)

        
        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0] 

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N) - 2
        if N < 0 :
            niveaux_temps[z[0]] = 2**0
        else :
            niveaux_temps[z[0]] = 2**N            
        ##if (N < 6):N=0
        ##else:N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2	
        #niveaux_temps['cart.294'] = 2

        #niveaux_temps[z[0]] = 2**N

        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))
  

    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                if interpDataType==1 and dictOfADT is not None:
                                    hook0.append(dictOfADT[zdnrname])
            if interpDataType == 0: hook0 = None

            X._setInterpData(zrcv, dnrZones, nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                             interpDataType=interpDataType, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]
            

            for nod in range(len(dnrZones)):

                print('zdnr= ', dnrZones[nod][0])

                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=numpy.int32)
                prangedonor = numpy.zeros(6,dtype=numpy.int32)
                profondeur=numpy.zeros(1,dtype=numpy.int32)
                dirD=numpy.zeros(1,dtype=numpy.int32)
                dirR=numpy.zeros(1,dtype=numpy.int32)
              
                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]
                tc[2][nobd][2][nozd] = dnrZones[nod]

                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo=numpy.zeros(3,dtype=numpy.int32)#XOD.getTransfo(dnrZones[nod],zrcv)

                connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1]+1,dim__[2]+1,dim__[3]+1)


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1]+1 , dim__[2]+1 , dim__[3]+1)

                NMratio = numpy.zeros(3,dtype=numpy.int32)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1] 
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])
                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
                    if interpDataType == 1: hook0 = []
                    else: hook0 = None
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if interpDataType == 1:
                                            if zdnrname not in dictOfADT: 
                                                HOOKADT = C.createHook(zdnr, 'adt')
                                                dictOfADT[zdnrname] = HOOKADT
                                            hook0.append(dictOfADT[zdnrname])

                                        dnrZones.append(zdnr)
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             nature=1, penalty=1, loc='centers', storage='inverse',  
                                             interpDataType=interpDataType, hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL)

                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=numpy.int32)
                        prangedonor = numpy.zeros(6,dtype=numpy.int32)
                        profondeur=numpy.zeros(1,dtype=numpy.int32)
                        dirD=numpy.zeros(1,dtype=numpy.int32)
                        dirR=numpy.zeros(1,dtype=numpy.int32)
              
                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]

  
                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=numpy.zeros(3,dtype=numpy.int32)#XOD.getTransfo(dnrZones[nod],zrcv)

                        connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1]+1,dim__[2]+1,dim__[3]+1)

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1]+1 , dim__[2]+1 , dim__[3]+1)

                        NMratio = numpy.zeros(3,dtype=numpy.int32)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                        print('LEVELS= ', levelrcv, leveldnr)


        if dictOfADT is not None: 
            for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

#=============================================================================
def doInterp3(t, tc, tbb, tb=None, typeI='ID', dim=3, dictOfADT=None, frontType=0, depth=2, IBCType=1, interpDataType=1):    
    

    ReferenceState = Internal.getNodeFromType2(t,'ReferenceState_t')

    bases  = Internal.getNodesFromType1(t     , 'CGNSBase_t') 
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension') 
    dimPb   = Internal.getValue(dimmm)
    dxmax = 0.0


    zones = Internal.getZones(t)

    dico_dx = {}
    dico_dy = {}
    dico_dz = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]
  
        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])
        
        dico_dx[z[0]] = dxx
        dico_dy[z[0]] = dyy
        if dimPb == 2:dico_dz[z[0]] = 1
        else : dico_dz[z[0]] = dzz

        if dimPb == 2:dzz=max(dxx,dyy)

        dx = min(dxx,dyy,dzz)
        if dx > dxmax:dxmax=dx

    niveaux_temps = {} 
    cx = {}

    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])

        if dimPb == 2:dzz=max(dxx,dyy)
        
        dx = min(dxx,dyy,dzz)

        #cx[z[0]]= coordx[1,0,0] 

        N = math.log(dxmax/dx)/math.log(2.0)
        N = round(N)
        if (N < 6):N=0
        else:N=1
        #if (cx[z[0]] < 0.98): niveaux_temps[z[0]] = 2**N
        #else:  niveaux_temps[z[0]] = 1
        #niveaux_temps['cart.264'] = 2	
        #niveaux_temps['cart.294'] = 2

        niveaux_temps[z[0]] = 2**N


        print(niveaux_temps[z[0]])
        #print(round(dxmax/dx))

    
    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
        

        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
               zrcvname = zrcv[0]; rcvZones.append(zrcv)

        #dico={}
        #for zrcv in Internal.getZones(t):
              # listofjoins = Internal.getNodesFromType2(zrcv, 'GridConnectivity_t')
              # if listofjoins is not None: 
               #    prange_list=[]
               #    dir_list=[]
               #    for join in listofjoins:
               #        prange_ = Internal.getNodeFromName1(join,'PointRange')[1]
                       #dirR = CGC.getDirBorderStruct__(prange_,dimPb)
                       #dir_list.append(dirR)
                       #print 'prange_= ', prange_
                #       for i in range(3):
                #           if prange_[i,1] == prange_[i,0] and prange_[i,1] != 1:
                #               prange_[i,1] =  prange_[i,1]-1
                #               prange_[i,0] =  prange_[i,0]-1
                #           elif prange_[i,1] != prange_[i,0] and prange_[i,1] != 1 :
                 #              prange_[i,1] =  prange_[i,1]-1
                 #      prange_=numpy.reshape(prange_,6)
                 #      prange_list.append(prange_)
                 #  dico[zrcv[0]]=prange_list
                   #dico[zrcv[0]]=dir_list
                   # print prange_, zrcv[0]   


        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in range(nbZonesChim):
            zrcv = rcvZones[nozr]
            dim_ = Internal.getZoneDim(zrcv)
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in range(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in range(len(tc[2][nobd][2])):
                        zdnr = tc[2][nobd][2][nozd]
                        if zdnr[3] == 'Zone_t':
                            zdnrname = zdnr[0]
                            if zdnrname in intersectionsDict[zrcvname]:
                                nobOfDnrBases.append(nobd)
                                nobOfDnrZones.append(nozd)
                                dnrZones.append(zdnr)
                                hook0.append(dictOfADT[zdnrname])

            dnrZones = X.setInterpData(zrcv,dnrZones,nature=1,penalty=1,loc='centers',storage='inverse',sameName=1,\
                                       hook=hook0, itype='chimera')


            levelrcv = niveaux_temps[zrcv[0]]

            for nod in range(len(dnrZones)):

                dim__ = Internal.getZoneDim(dnrZones[nod])
                prange = numpy.zeros(6,dtype=numpy.int32)
                prangedonor = numpy.zeros(6,dtype=numpy.int32)
                profondeur=numpy.zeros(1,dtype=numpy.int32)
                dirD=numpy.zeros(1,dtype=numpy.int32)
                dirR=numpy.zeros(1,dtype=numpy.int32)
              
                plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                leveldnr = niveaux_temps[dnrZones[nod][0]]

                nobd = nobOfDnrBases[nod]
                nozd = nobOfDnrZones[nod]

                tc[2][nobd][2][nozd] = dnrZones[nod]

 
                prangebis=numpy.reshape(prange,6)

                info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                transfo = XOD.getTransfo(dnrZones[nod],zrcv)

                connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1]+1,dim__[2]+1,dim__[3]+1)


                #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1]+1 , dim__[2]+1 , dim__[3]+1)

                NMratio = numpy.zeros(3,dtype=numpy.int32)
                NMratio[0]=1
                NMratio[1]=1
                NMratio[2]=1

                info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])



    elif typeI == 'IBCD':
        # detection des zones IBC
        zonesRIBC = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellNIBC')==2.:
                zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

        if zonesRIBC == []: return tc

        print('Building the IBM front.')

        front = getIBMFront(tc, 'cellNFront', dim, frontType)
        # Sortie du front pour debug
        C.convertPyTree2File(front, 'front.cgns')

        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType, Reynolds=Reynolds, yplus=yplus)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1] 
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])

                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]; hook0 = []
                    for nobd in range(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in range(len(tc[2][nobd][2])):
                                zdnr = tc[2][nobd][2][nozd]
                                if zdnr[3] == 'Zone_t':
                                    zdnrname = zdnr[0]
                                    zbb = tbb[2][nobd][2][nozd]
                                    bba = C.getFields(Internal.__GridCoordinates__,zbb)[0]
                                    if Generator.bboxIntersection(interpPtsBB,bba,isBB=True) == 1:
                                        if zdnrname not in dictOfADT: 
                                            HOOKADT = C.createHook(zdnr, 'adt')
                                            dictOfADT[zdnrname] = HOOKADT
                                        dnrZones.append(zdnr)
                                        hook0.append(dictOfADT[zdnrname])
                                        nobOfDnrBases.append(nobd)
                                        nobOfDnrZones.append(nozd)

                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr], \
                                             loc='centers', storage='inverse',  hook=hook0, dim=dim, ReferenceState=ReferenceState, bcType=ibcTypeL)
                    nozr += 1

                    levelrcv = niveaux_temps[zrcv[0]]

                    for nod in range(len(dnrZones)):

                        dim__ = Internal.getZoneDim(dnrZones[nod])
                        prange = numpy.zeros(6,dtype=numpy.int32)
                        prangedonor = numpy.zeros(6,dtype=numpy.int32)
                        profondeur=numpy.zeros(1,dtype=numpy.int32)
                        dirD=numpy.zeros(1,dtype=numpy.int32)
                        dirR=numpy.zeros(1,dtype=numpy.int32)
              
                        plist = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][2][1]
                        plistdnr = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][3][1]
                        coeff = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][4][1]
                        typ = dnrZones[nod][2][len(dnrZones[nod][2])-1][2][5][1]

                        leveldnr = niveaux_temps[dnrZones[nod][0]]

                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]

                        tc[2][nobd][2][nozd] = dnrZones[nod]

  
                        prangebis=numpy.reshape(prange,6)

                        info = dnrZones[nod][2][len(dnrZones[nod][2])-1]
                        info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                        transfo=XOD.getTransfo(dnrZones[nod],zrcv)

                        connector.indiceToCoord2(plist,prangedonor,transfo,profondeur,dirD,typ,dirR,plist.size,dim__[1]+1,dim__[2]+1,dim__[3]+1)

                        #connector.correctCoeffList(plist, coeff, typ, plist.size , dim__[1]+1 , dim__[2]+1 , dim__[3]+1)

                        NMratio = numpy.zeros(3,dtype=numpy.int32)
                        NMratio[0]=1
                        NMratio[1]=1
                        NMratio[2]=1

                        info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                        info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                        info[2].append(['PointPivot', transfo , [], 'IndexArray_t'])
                        info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                        print('LEVELS= ', levelrcv, leveldnr)



        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])

        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    return tc

#=============================================================================
# Performs the full IBM preprocessing using overlapping Cartesian grids
# interpDataType = 1 : interpolation par preconditionnement par ADT
# interpDataType = 0 : interpolation optimisees sur grilles cartesiennes
# frontType 0, 1, 2
#=============================================================================
def prepareIBMData(t, tbody, DEPTH=2, loc='centers', frontType=1, interpDataType=0, smoothing=False, yplus=100., wallAdapt=None,isLBM=False,LBMQ=False,isPrintDebug=False):
    tb =  Internal.copyRef(tbody)

    # tb: fournit model et dimension
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent, LBMLaminar
    model = Internal.getValue(model)

    if model == 'Euler': IBCType =-1
    elif model == 'LBMLaminar':
        IBCType =-1      
        if LBMQ == True:
            IBCType =1
    else: IBCType = 1 # Points cibles externes
    if loc == 'nodes':
        raise NotImplemented("prepareIBMData: prepareIBMData at nodes not yet implemented.")

    #------------------------
    # Ghost cells (overlaps)
    #------------------------
    X._applyBCOverlaps(t, depth=DEPTH,loc='centers',val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    #------------------------
    # Blanking IBM
    #------------------------
    C._initVars(t,'centers:cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)

    t = blankByIBCBodies(t, tb, 'centers', dimPb)
    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #-----------------------------------------
    # calcul de la normale et distance signee
    #-----------------------------------------
    COMPDIST = False # distance deja calculee ou non
    if Internal.getNodeFromName(t, 'TurbulentDistance') is None: COMPDIST=True
    if COMPDIST:
        print('Computing distance field...')
        DTW._distance2Walls(t, tb, loc='centers', type='ortho', signed=0)
    else: pass
    _signDistance(t)

    #-----------------------------------------
    # Pts IBC
    #-----------------------------------------
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    Reynolds = Internal.getNodeFromName(tb, 'Reynolds')
    if Reynolds is not None: Reynolds = Internal.getValue(Reynolds)
    else: Reynolds = 6.e6
    if Reynolds < 1.e6: frontType = 1
    if frontType != 42:
        if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        elif IBCType == 1:
            X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
            if frontType < 2:
                X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
            else:
                DEPTHL=DEPTH+1
                X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0,loc='centers',cellNName='cellN',addGC=False)
                #cree des pts extrapoles supplementaires
                # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
        else:
            raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    else:
        # determination des pts IBC en fonction de la distance a la paroi
        # cree d'abord un front de type 1 pour assurer que l'adaptation ne descende pas en dessous de la limite...
        C._initVars(t,'{centers:cellNMin}={centers:cellNIBC}')
        if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)
        elif IBCType == 1: X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNMin',addGC=False) # pour les gradients
        X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellNMin',addGC=False)

        for z in Internal.getZones(t):
            h = abs(C.getValue(z,'CoordinateX',0)-C.getValue(z,'CoordinateX',1))
            if yplus > 0.:
                height = computeModelisationHeight(Re=Reynolds, yplus=yplus)
            else:
                height = computeBestModelisationHeight(Re=Reynolds, h=h) # meilleur compromis entre hauteur entre le snear et la hauteur de modelisation
            C._initVars(z,'{centers:cellN}=({centers:TurbulentDistance}>%20.16g)+(2*({centers:TurbulentDistance}<=%20.16g)*({centers:TurbulentDistance}>0))'%(height,height))

        # Si wallAdapt, utilisation de la solution precedente pour ne garder que les pts cibles tq y+PC <= y+ref : rapproche le front de la paroi (utile proche bord d'attaque)
        # Attention, le fichier d'adaptation doit etre un nuage de points dont les coordonnees correspondent aux points cibles (Cf modification dans extractIBMWallFields avec tb=None)
        if wallAdapt is not None:
            C._initVars(t,'{centers:yplus}=100000.')
            w = C.convertFile2PyTree(wallAdapt)
            yplus_w = Internal.getNodeFromName(w, 'yplus')[1]
            total = len(Internal.getZones(t))
            cpt = 1
            for z in Internal.getZones(t):
                print("{} / {}".format(cpt, total))
                cellN = Internal.getNodeFromName(z,'cellN')[1]
                if 2 in cellN:
                    yplus_z = Internal.getNodeFromName(z, 'yplus')[1]
                    original_shape = numpy.shape(yplus_z)
                    yplus_z = yplus_z.ravel(order='K')
                    hook = C.createHook(z, function='elementCenters')
                    nodes = C.identifyNodes(hook, w)
                    for pos, node in enumerate(nodes):
                        if node != -1:
                            yplus_z[node-1]=yplus_w[pos]
                    yplus_z = numpy.reshape(yplus_z, original_shape)
                cpt += 1

            C._initVars(t,'{centers:cellN}=({centers:cellN}>0) * ( (({centers:cellN}) * ({centers:yplus}<=%20.16g)) + ({centers:yplus}>%20.16g) )'%(yplus,yplus))
            
        # Securite finale, on aura au min deux rangees de points cibles
        # def maximum(x1, x2): return max(x1,x2)
        # C._initVars(t,'centers:cellN', maximum, ['centers:cellN','centers:cellNMin'])
        # del maximum 
        C._initVars(t, '{centers:cellN} = maximum({centers:cellN}, {centers:cellNMin})')
        C._rmVars(t,['centers:yplus', 'centers:cellNMin'])


        # permet de ne garder que les deux rangees superieures (prises en compte par le solveur), mais complique l'adaptation suivante et la visualisation
        # X._maximizeBlankedCells(t, depth=2, addGC=False)

    _removeBlankedGrids(t, loc='centers')
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
       dims = Internal.getZoneDim(i)
       npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #------------------------------------------------------------------------
    # Nature des points en fonction de leur nature Chimere et leur nature IBC
    #------------------------------------------------------------------------
    # -3 : agit comme un point masque - non donneur pour le type de point
    #  3  : agit comme donneur uniquement
    # updateNatureForIBM: modifie cellNChim, cellNFront, cellNIBM
    # cellNChim=-3, si cellNIBC=0 (masque)
    if IBCType == 1: # Points corriges IBM externes
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)

    else: # EN 2 PARTIES : NECESSITE LE TRANSFERT DU FRONT PAR INTERPOLATION, QUI EST CALCULEE APRES
        print('Euler: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:dummy}={centers:cellN}') # sauvegarde
        C._initVars(t,'{centers:cellN}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellN}>0.5, {centers:cellN}<1.5)')
        C._cpVars(t,'centers:dummy',t,'centers:cellN')
        C._rmVars(t, ['centers:dummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    #------------------------------------------------------------------------
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tc = C.node2Center(t)
    FSN = Internal.getNodesFromName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN,'cellNFront')
    Internal._rmNodesByName(FSN,'cellNIBC')
    Internal._rmNodesByName(FSN, "TurbulentDistance")

    tbb = G.BB(tc)

    # Creation du dictionnaire des ADT pour les raccords
    if interpDataType == 1:
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
    else: dictOfADT = None
    print('Interpolations Chimere.')
    tc = doInterp(t, tc, tbb, tb=None, typeI='ID', dim=dimPb,
                  interpDataType=interpDataType, dictOfADT=dictOfADT)
    if dictOfADT is not None:
        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    for zc in Internal.getZones(tc):
        cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    if frontType == 2: 
        _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType) 
        if smoothing and dimPb==2: _smoothImageFront(t, tc); _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType) 

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        for zc in Internal.getZones(tc):
            cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    # C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t,'centers:TurbulentDistance')
    print('Building the IBM front.')
    front = getIBMFront(tc, 'cellNFront', dimPb, frontType)
    if isPrintDebug:
        C.convertPyTree2File(front, 'IB_front.cgns')
    C.convertPyTree2File(front, 'front.cgns')
    print('Interpolations IBM')
    tc = doInterp(t,tc,tbb, tb=tb,typeI='IBCD',dim=dimPb, dictOfADT=None, front=front, frontType=frontType, depth=DEPTH, IBCType=IBCType, interpDataType=interpDataType, Reynolds=Reynolds, yplus=yplus,isLBM=isLBM)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBCDnr']
    varsRM += ['centers:cellNChim','centers:cellNIBC']
    C._rmVars(t, varsRM)
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance'])
    #----------
    # SORTIE
    #----------
    return t, tc

#=============================================================================
# Performs the full IBM preprocessing using overlapping Cartesian grids
# interpDataType = 1 : interpolation par preconditionnement par ADT
# interpDataType = 0 : interpolation optimisees sur grilles cartesiennes
# frontType 0, 1, 2
#=============================================================================
def prepareIBMData2(t, tbody, DEPTH=2, loc='centers', frontType=1, inv=False, interpDataType=1):
    tb =  Internal.copyRef(tbody)

    # tb: fournit model et dimension
    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None: raise ValueError('prepareIBMData: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)

    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('prepareIBMData: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if model == 'Euler': IBCType =-1
    else: IBCType = 1 # Points cibles externes
    if loc == 'nodes':
        raise NotImplemented("prepareIBMData: prepareIBMData at nodes not yet implemented.")

    #------------------------
    # Ghost cells (overlaps)
    #------------------------
    X._applyBCOverlaps(t, depth=DEPTH,loc='centers',val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    #------------------------
    # Blanking IBM
    #------------------------
    C._initVars(t,'centers:cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)

    t = blankByIBCBodies(t,tb,'centers',dimPb)
    if not inv: C._initVars(t,'{centers:cellNIBC}={centers:cellN}')
    if inv: C._initVars(t,'{centers:cellNIBC}=1-{centers:cellN}') # ecoulement interne


      
    #-----------------------------------------
    # calcul de la normale et distance signee
    #-----------------------------------------
    COMPDIST = False # distance deja calculee ou non
    if Internal.getNodeFromName(t, 'TurbulentDistance') is None: COMPDIST=True
    if COMPDIST:
        print('Computing distance field...')
        DTW._distance2Walls(t,tb,loc='centers',type='ortho',signed=0)
    else: pass
    _signDistance(t)

    #-----------------------------------------
    # Pts IBC
    #-----------------------------------------
    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
    elif IBCType == 1:
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
        if frontType < 2:
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        else:
            DEPTHL=DEPTH+1
            X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0, loc='centers',cellNName='cellN',addGC=False)
            #cree des pts extrapoles supplementaires
            # _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
    else:
        raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    _removeBlankedGrids(t, loc='centers')
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
       dims = Internal.getZoneDim(i)
       npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    #------------------------------------------------------------------------
    # Nature des points en fonction de leur nature Chimere et leur nature IBC
    #------------------------------------------------------------------------
    # -3 : agit comme un point masque - non donneur pour le type de point
    #  3  : agit comme donneur uniquement
    # updateNatureForIBM: modifie cellNChim, cellNFront, cellNIBM
    # cellNChim=-3, si cellNIBC=0 (masque)
    if IBCType == 1: # Points corriges IBM externes
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)

    else: # EN 2 PARTIES : NECESSITE LE TRANSFERT DU FRONT PAR INTERPOLATION, QUI EST CALCULEE APRES
        print('Euler: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:dummy}={centers:cellN}') # sauvegarde
        C._initVars(t,'{centers:cellN}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellN}>0.5, {centers:cellN}<1.5)')
        C._cpVars(t,'centers:dummy',t,'centers:cellN')
        C._rmVars(t, ['centers:dummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    #------------------------------------------------------------------------
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tc = C.node2Center(t)
    FSN = Internal.getNodesFromName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN,'cellNFront')
    Internal._rmNodesByName(FSN,'cellNIBC')
    Internal._rmNodesByName(FSN, "TurbulentDistance")

    tbb = G.BB(tc)

    # Creation du dictionnaire des ADT pour les raccords
    if interpDataType == 1:
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
    else: dictOfADT = None
    print('Interpolations Chimere.')
    tc = doInterp2(t, tc, tbb, tb=None, typeI='ID', dim=dimPb,
                  interpDataType=interpDataType, dictOfADT=dictOfADT)
    if dictOfADT is not None:
        for dnrname in dictOfADT: C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    for zc in Internal.getZones(tc):
        cellNFront = Internal.getNodeFromName2(zc,'cellNFront')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                XOD._setInterpTransfers(t,zc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    if frontType==2 or frontType==3: _pushBackImageFront2(t, tc, tbb, interpDataType=interpDataType)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t,'centers:TurbulentDistance')
    print('Building the IBM front.')
    front = getIBMFront(tc, 'cellNFront', dimPb, frontType)
    print('Interpolations IBM')
    tc = doInterp2(t,tc,tbb, tb=tb,typeI='IBCD',dim=dimPb, dictOfADT=None, front=front, frontType=frontType, depth=DEPTH, IBCType=IBCType, interpDataType=interpDataType)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBCDnr']
    varsRM += ['centers:cellNChim','centers:cellNIBC']
    C._rmVars(t, varsRM)
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance'])
    #----------
    # SORTIE
    #----------
    return t, tc


#=============================================================================
# Performs specific treatment for frontType=2 or frontType=3
# Performs specific treatment for frontType=2
#=============================================================================
def _pushBackImageFront2(t, tc, tbb, interpDataType=1):
    intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)
    C._initVars(t,'{centers:cellNFront_origin}={centers:cellNFront}')
    C._initVars(t,'{centers:cellNIBC_origin}={centers:cellNIBC}')
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._cpVars(t,'centers:cellNFront_origin',tc,'cellNFront_origin')
    C._cpVars(t,'centers:cellNIBC_origin',tc,'cellNIBC_origin')

    C._initVars(t,'{centers:cellNFront2}=1.-({centers:cellNFront}<1.)*(abs({centers:cellNChim})>1.)')

    for z in Internal.getZones(t):
        cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
        if cellNFront != []:
            cellNFront = cellNFront[1]
            sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
            sizeOne =  int(numpy.sum(cellNFront))
            if sizeOne < sizeTot:
                X._setHoleInterpolatedPoints(z, depth=1, dir=0, loc='centers',cellNName='cellNFront2',addGC=False)
                res = X.getInterpolatedPoints(z,loc='centers', cellNName='cellNFront2') #indices,X,Y,Z
                if res is not None:
                    indicesI = res[0]
                    XI = res[1]; YI = res[2]; ZI = res[3]
                    allInterpFields=[]
                    for zc in Internal.getZones(tc):
                        if zc[0] in intersectionsDict[z[0]]:
                            if interpDataType==1: HOOKADT = C.createHook(zc, 'adt')
                            else: HOOKADT = None
                            fields = X.transferFields(zc, XI, YI, ZI, hook=HOOKADT, variables=['cellNFront_origin','cellNIBC_origin'], interpDataType=interpDataType)
                            allInterpFields.append(fields)
                            if interpDataType == 1: C.freeHook(HOOKADT)
                    if allInterpFields != []:
                        C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol',verbose=False)
                        # C._initVars(z,'{centers:cellNFront}=({centers:cellNFront}>0.5)') #ancienne version
                        C._initVars(z,'{centers:cellNFront}={centers:cellNFront}*({centers:cellNFront_origin}>0.5)') # Modification du Front uniquement lorsque celui-ci est repousse
                        # i.e. if cellNFront_origin == 0 and cellNFront == 1 => cellNfront = 0

                        C._initVars(z,'{centers:cellNIBC}={centers:cellNIBC}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                            + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                        # i.e. if cellNChim == 1 and cellNIBC_origin == 2 => cellNIBC = 2

                        C._initVars(z,'{centers:cellNIBCDnr}={centers:cellNIBCDnr}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                            + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                        # le front donneur doit etre egalement maj sur l'exemple du champ cellNIBC

        C._rmVars(z,['centers:cellNFront2'])
        C._rmVars(z,['centers:cellNFront_origin'])
        C._rmVars(z,['centers:cellNIBC_origin'])

    C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    return None      

#=============================================================================
# Performs smoothing after _pushBackImageFront2 
#=============================================================================
def _smoothImageFront(t, tc, dimPb=2):
    for z in Internal.getZones(t):
        cellNFront  = Internal.getNodeFromName(z,'cellNFront')
        cellNIBC    = Internal.getNodeFromName(z,'cellNIBC')
        cellNIBCDnr = Internal.getNodeFromName(z,'cellNIBCDnr')
        dim         = Internal.getZoneDim(z)
        k           = 0
        if cellNIBC != []:
                cellNIBC = cellNIBC[1]
                sizeTot = cellNIBC.shape[0]*cellNIBC.shape[1]*cellNIBC.shape[2]
                sizeOne =  int(numpy.sum(cellNIBC))
                if sizeOne < sizeTot:
                    cellNFront  = cellNFront[1]
                    cellNIBCDnr = cellNIBCDnr[1]
                    for i in range(1, int(dim[1])-2):
                        for j in range(1, int(dim[2])-2):
                            if cellNIBC[i,j,k] == 1 and cellNIBC[i,j-1,k] == 2:
                                if cellNIBC[i-1,j+1,k] == 2:
                                    cellNFront[i,j,k]    = 0
                                    cellNIBC[i,j,k]      = 2
                                    cellNIBCDnr[i,j,k]   = 2

                        for j in range(int(dim[2])-3, 0, -1):
                            if cellNIBC[i,j,k] == 1 and cellNIBC[i,j+1,k] == 2:
                                if cellNIBC[i-1,j-1,k] == 2:
                                    cellNFront[i,j,k]    = 0
                                    cellNIBC[i,j,k]      = 2
                                    cellNIBCDnr[i,j,k]   = 2

                    for i in range(int(dim[1])-3, 0, -1):
                        for j in range(1, int(dim[2])-2):
                            if cellNIBC[i,j,k] == 1 and cellNIBC[i,j-1,k] == 2:
                                if cellNIBC[i+1,j+1,k] == 2:
                                    cellNFront[i,j,k]    = 0
                                    cellNIBC[i,j,k]      = 2
                                    cellNIBCDnr[i,j,k]   = 2

                        for j in range(int(dim[2])-3, 0, -1):
                            if cellNIBC[i,j,k] == 1 and cellNIBC[i,j+1,k] == 2:
                                if cellNIBC[i+1,j-1,k] == 2:
                                    cellNFront[i,j,k]    = 0
                                    cellNIBC[i,j,k]      = 2
                                    cellNIBCDnr[i,j,k]   = 2

    C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    return None  


#=============================================================================
# Compute the height of the modelisation zone for the turbulent boundary layer
#=============================================================================
def computeModelisationHeight(Re, Cf_law='ANSYS', yplus=100., L=1.):
    if Cf_law == 'ANSYS':
        def compute_Cf(Re):
            return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        def compute_Cf(Re):
            return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        def compute_Cf(Re):
            return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        def compute_Cf(Re):
            return 1.46*Re**(-0.5)

    return (yplus*L*math.sqrt(2))/(Re*math.sqrt(compute_Cf(Re)))

#=============================================================================
# Compute the best height of the modelisation zone for the turbulent boundary 
# layer for a given snear
#=============================================================================
def computeBestModelisationHeight(Re, h, Cf_law='ANSYS', L=1., q=1.2):
    if Cf_law == 'ANSYS':
        def compute_Cf(Re):
            return 0.058*Re**(-0.2)
    elif Cf_law == 'PW':
        def compute_Cf(Re):
            return 0.026*Re**(-1/7.)
    elif Cf_law == 'PipeDiameter':
        def compute_Cf(Re):
            return 0.079*Re**(-0.25)
    elif Cf_law == 'Laminar':
        def compute_Cf(Re):
            return 1.46*Re**(-0.5)

    h0 = (L*math.sqrt(2))/(Re*math.sqrt(compute_Cf(Re)))

    return (h0-q*h)/(1.-q) 

#=============================================================================
# Extraction des infos pour le post traitement
# if td=None: return the cloud of points
# else interpolate on td
# input famZones: family of subregions to extract as a list ['FAM1','FAM2,...]
#=============================================================================
def extractIBMWallFields(tc, tb=None, coordRef='wall', famZones=[]):
    xwNP = []; ywNP = []; zwNP = []
    xiNP = []; yiNP = []; ziNP = []
    xcNP = []; ycNP = []; zcNP = []
    pressNP = []; utauNP = []; yplusNP = []; densNP = []
    vxNP = []; vyNP = []; vzNP = []
    KCurvNP = []
    dictOfFamilies={}
    if famZones != []:
        out = []
        allIBCD = Internal.getNodesFromName(tc,"IBCD_*")
        for zsr in allIBCD:
            fam = Internal.getNodeFromType(zsr,'FamilyName_t')
            if fam is not None:
                famName=Internal.getValue(fam)
                if famName in famZones:
                    if famName not in dictOfFamilies: dictOfFamilies[famName]=[zsr]
                    else: dictOfFamilies[famName]+=[zsr]


        # Creation of a single zone
        zsize = numpy.empty((1,3), numpy.int32, order='F')
        zsize[0,0] = 1; zsize[0,1] = 0; zsize[0,2] = 0
        dictOfZoneFamilies={}
        for z in Internal.getZones(tb):
            famName = Internal.getNodeFromType(z,'FamilyName_t')
            if famName is not None:
                famName = Internal.getValue(famName)
                if famName in famZones:
                    if famName not in dictOfZoneFamilies: dictOfZoneFamilies[famName]=[z]
                    else: dictOfZoneFamilies[famName]+=[z]
        for famName in dictOfFamilies:
            zd = Internal.newZone(name='ZIBC_%s'%famName,zsize=zsize,ztype='Unstructured')
            zd[2] += dictOfFamilies[famName]
            tb2 = None
            if tb is not None:
                zones = dictOfZoneFamilies[famName]
                if zones!=[]: tb2=C.newPyTree(['Base']); tb2[2][1][2]=zones

            zd = extractIBMWallFields(zd, tb=tb2, coordRef=coordRef, famZones=[])
            out+=Internal.getZones(zd)
        return out
    
    else:
        allZSR = Internal.getNodesFromType(tc,'ZoneSubRegion_t')
        if allZSR != []:
            allIBCD = Internal.getNodesFromName(allZSR,"IBCD_*")
            for IBCD in allIBCD:
                xPW = Internal.getNodeFromName1(IBCD,"CoordinateX_PW")[1]
                yPW = Internal.getNodeFromName1(IBCD,"CoordinateY_PW")[1]
                zPW = Internal.getNodeFromName1(IBCD,"CoordinateZ_PW")[1]
                xwNP.append(xPW); ywNP.append(yPW); zwNP.append(zPW)
                
                xPI = Internal.getNodeFromName1(IBCD,"CoordinateX_PI")[1]
                yPI = Internal.getNodeFromName1(IBCD,"CoordinateY_PI")[1]
                zPI = Internal.getNodeFromName1(IBCD,"CoordinateZ_PI")[1]
                xiNP.append(xPI); yiNP.append(yPI); ziNP.append(zPI)
                
                xPC = Internal.getNodeFromName1(IBCD,"CoordinateX_PC")[1]
                yPC = Internal.getNodeFromName1(IBCD,"CoordinateY_PC")[1]
                zPC = Internal.getNodeFromName1(IBCD,"CoordinateZ_PC")[1]
                xcNP.append(xPC); ycNP.append(yPC); zcNP.append(zPC)
                                
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

                KCURVW = Internal.getNodeFromName1(IBCD, XOD.__KCURV__)
                if KCURVW is not None: KCurvNP.append(KCURVW[1])

    if pressNP == []: return None
    else:
        pressNP = numpy.concatenate(pressNP)
        densNP = numpy.concatenate(densNP)
        if utauNP != []: utauNP = numpy.concatenate(utauNP)
        if yplusNP != []: yplusNP = numpy.concatenate(yplusNP)
        if vxNP != []: vxNP = numpy.concatenate(vxNP)
        if vyNP != []: vyNP = numpy.concatenate(vyNP)
        if vzNP != []: vzNP = numpy.concatenate(vzNP)
        if KCurvNP != []: KCurvNP = numpy.concatenate(KCurvNP)

        xwNP = numpy.concatenate(xwNP)
        ywNP = numpy.concatenate(ywNP)
        zwNP = numpy.concatenate(zwNP)

        xiNP = numpy.concatenate(xiNP)
        yiNP = numpy.concatenate(yiNP)
        ziNP = numpy.concatenate(ziNP)

        xcNP = numpy.concatenate(xcNP)
        ycNP = numpy.concatenate(ycNP)
        zcNP = numpy.concatenate(zcNP)

    # Creation d une seule zone
    zsize = numpy.empty((1,3), numpy.int32, order='F')
    zsize[0,0] = xwNP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=z)
    if coordRef == 'cible':
        coordx = ['CoordinateX',xcNP,[],'DataArray_t']
        coordy = ['CoordinateY',ycNP,[],'DataArray_t']
        coordz = ['CoordinateZ',zcNP,[],'DataArray_t']
    elif coordRef == 'wall':
        coordx = ['CoordinateX',xwNP,[],'DataArray_t']
        coordy = ['CoordinateY',ywNP,[],'DataArray_t']
        coordz = ['CoordinateZ',zwNP,[],'DataArray_t']
    else:
        coordx = ['CoordinateX',xwNP,[],'DataArray_t']
        coordy = ['CoordinateY',ywNP,[],'DataArray_t']
        coordz = ['CoordinateZ',zwNP,[],'DataArray_t']
    gc[2] = [coordx,coordy,coordz]
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                   gridLocation='Vertex', parent=z)
    FSN[2].append([XOD.__PRESSURE__,pressNP, [],'DataArray_t'])
    FSN[2].append([XOD.__DENSITY__,densNP, [],'DataArray_t'])

    FSN[2].append(["CoordinateX_PW",xwNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateY_PW",ywNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateZ_PW",zwNP, [],'DataArray_t'])

    FSN[2].append(["CoordinateX_PI",xiNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateY_PI",yiNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateZ_PI",ziNP, [],'DataArray_t'])

    FSN[2].append(["CoordinateX_PC",xcNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateY_PC",ycNP, [],'DataArray_t'])
    FSN[2].append(["CoordinateZ_PC",zcNP, [],'DataArray_t'])

    utauPresent = 0; yplusPresent = 0
    if utauNP != []:
        utauPresent = 1
        FSN[2].append([XOD.__UTAU__,utauNP, [],'DataArray_t'])
    if yplusNP != []:
        yplusPresent = 1
        FSN[2].append([XOD.__YPLUS__,yplusNP, [],'DataArray_t'])
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
        if vxPresent==1: 
            C._initVars(td,XOD.__VELOCITYX__,0.)
            C._initVars(td,XOD.__VELOCITYY__,0.)
            C._initVars(td,XOD.__VELOCITYZ__,0.)
        if kcurvPresent==1:
            C._initVars(td,XOD.__KCURV__,0.)
        td = P.projectCloudSolution(z, td, dim=dimPb)
        return td

#=============================================================================
# Extraction des pts IBM: retourne un arbre avec les coordonnees des
# pts IBM a corriger, paroi, miroirs
#=============================================================================
def extractIBMInfo(tc):
    XPC={}; YPC={}; ZPC={}
    XPW={}; YPW={}; ZPW={}
    XPI={}; YPI={}; ZPI={}
    Zones = []
    for z in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromName(z, "IBCD_*")
        for IBCD in allIBCD:
            znames = Internal.getValue(IBCD)
            Zones.append(znames)
            xPC = Internal.getNodesFromName(IBCD,"CoordinateX_PC")[0][1]
            yPC = Internal.getNodesFromName(IBCD,"CoordinateY_PC")[0][1]
            zPC = Internal.getNodesFromName(IBCD,"CoordinateZ_PC")[0][1]
            xPI = Internal.getNodesFromName(IBCD,"CoordinateX_PI")[0][1]
            yPI = Internal.getNodesFromName(IBCD,"CoordinateY_PI")[0][1]
            zPI = Internal.getNodesFromName(IBCD,"CoordinateZ_PI")[0][1]
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_PW")[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_PW")[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_PW")[0][1]
            if znames in XPW:
                a = numpy.concatenate((XPW[znames][0],xPW))
                XPW[znames] = [a]
                a = numpy.concatenate((YPW[znames][0],yPW))
                YPW[znames] = [a]
                a = numpy.concatenate((ZPW[znames][0],zPW))
                ZPW[znames] = [a]
                a = numpy.concatenate((XPI[znames][0],xPI))
                XPI[znames] = [a]
                a = numpy.concatenate((YPI[znames][0],yPI))
                YPI[znames] = [a]
                a = numpy.concatenate((ZPI[znames][0],zPI))
                ZPI[znames] = [a]
                a = numpy.concatenate((XPC[znames][0],xPC))
                XPC[znames] = [a]
                a = numpy.concatenate((YPC[znames][0],yPC))
                YPC[znames] = [a]
                a = numpy.concatenate((ZPC[znames][0],zPC))
                ZPC[znames] = [a]

            else:
                XPW[znames] = [xPW]
                YPW[znames] = [yPW]
                ZPW[znames] = [zPW]
                XPI[znames] = [xPI]
                YPI[znames] = [yPI]
                ZPI[znames] = [zPI]
                XPC[znames] = [xPC]
                YPC[znames] = [yPC]
                ZPC[znames] = [zPC]
    Zones = list(set(Zones))
    corrected = []; wall = []; interp = []
    t = C.newPyTree(['IBM','Wall','Image'])
    for zname in Zones:
        xPC = XPC[zname]; yPC = YPC[zname]; zPC = ZPC[zname]
        size = xPC[0].shape[0]
        coordxPC = ['CoordinateX',xPC[0],[],'DataArray_t']
        coordyPC = ['CoordinateY',yPC[0],[],'DataArray_t']
        coordzPC = ['CoordinateZ',zPC[0],[],'DataArray_t']
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'correctedPts_'+zname

        XPC0 = Internal.getNodeFromName(zone,'CoordinateX')
        parent,d = Internal.getParentOfNode(zone, XPC0)
        parent[2][d] = coordxPC

        YPC0 = Internal.getNodeFromName(zone,'CoordinateY')
        parent,d = Internal.getParentOfNode(zone, YPC0)
        parent[2][d] = coordyPC

        ZPC0 = Internal.getNodeFromName(zone,'CoordinateZ')
        parent,d = Internal.getParentOfNode(zone, ZPC0)
        parent[2][d] = coordzPC
        corrected.append(zone)
        #
        xPI = XPI[zname]; yPI = YPI[zname]; zPI = ZPI[zname]
        size = xPI[0].shape[0]
        coordxPI = ['CoordinateX',xPI[0],[],'DataArray_t']
        coordyPI = ['CoordinateY',yPI[0],[],'DataArray_t']
        coordzPI = ['CoordinateZ',zPI[0],[],'DataArray_t']
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'interpPts_'+zname

        XPI0 = Internal.getNodeFromName(zone,'CoordinateX')
        parent,d = Internal.getParentOfNode(zone, XPI0)
        parent[2][d] = coordxPI

        YPI0 = Internal.getNodeFromName(zone,'CoordinateY')
        parent,d = Internal.getParentOfNode(zone, YPI0)
        parent[2][d] = coordyPI

        ZPI0 = Internal.getNodeFromName(zone,'CoordinateZ')
        parent,d = Internal.getParentOfNode(zone, ZPI0)
        parent[2][d] = coordzPI
        interp.append(zone)

        xPW = XPW[zname];yPW = YPW[zname];zPW = ZPW[zname]
        size = xPW[0].shape[0]
        coordxPW = ['CoordinateX',xPW[0],[],'DataArray_t']
        coordyPW = ['CoordinateY',yPW[0],[],'DataArray_t']
        coordzPW = ['CoordinateZ',zPW[0],[],'DataArray_t']
        zone = G.cart((0,0,0),(1,1,1),(size,1,1))
        zone[0] = 'wallPts_'+zname

        XPW0 = Internal.getNodeFromName(zone,'CoordinateX')
        parent,d = Internal.getParentOfNode(zone, XPW0)
        parent[2][d] = coordxPW

        YPW0 = Internal.getNodeFromName(zone,'CoordinateY')
        parent,d = Internal.getParentOfNode(zone, YPW0)
        parent[2][d] = coordyPW

        ZPW0 = Internal.getNodeFromName(zone,'CoordinateZ')
        parent,d = Internal.getParentOfNode(zone, ZPW0)
        parent[2][d] = coordzPW
        wall.append(zone)

    t[2][1][2] = corrected; t[2][2][2] = wall; t[2][3][2] = interp
    t = C.convertArray2Node(t)
    return t
# --------------------------------------------------------------------------------
# Creation des zones 1D de nom 'Zone#IBCD_*' pour pouvoir etre identifie en post
# retourne tw avec une zone 1D par IBCD de depart
# Sert au post traitement paroi en parallele et instationnaire
def createIBMWZones(tc,variables=[]):
    import Converter.Mpi as Cmpi
    tw = C.newPyTree(['IBM_WALL'])
    for z in Internal.getZones(tc):
        ZSR = Internal.getNodesFromType2(z,'ZoneSubRegion_t')
        procNode = Internal.getNodeFromName(z,'proc')
        proc = -1
        if procNode is not None:
            proc = Internal.getValue(procNode)

        for IBCD in Internal.getNodesFromName(ZSR,"IBCD_*"):
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_PW")[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_PW")[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_PW")[0][1]
            nptsW = xPW.shape[0]
            zw = G.cart((0,0,0),(1,1,1),(nptsW,1,1))
            COORDX = Internal.getNodeFromName2(zw,'CoordinateX'); COORDX[1]=xPW
            COORDY = Internal.getNodeFromName2(zw,'CoordinateY'); COORDY[1]=yPW
            COORDZ = Internal.getNodeFromName2(zw,'CoordinateZ'); COORDZ[1]=zPW
            if variables != [] and variables != None:
                FSN = Internal.newFlowSolution(parent=zw)
                for varo in variables:
                    fieldV = Internal.getNodeFromName2(IBCD,varo)
                    if fieldV is not None: 
                        C._initVars(zw,varo,0.)
                        fieldW = Internal.getNodeFromName2(FSN,varo)
                        fieldW[1] = fieldV[1]

            zw[0]=z[0]+"#"+IBCD[0]
            if proc != -1: Cmpi._setProc(zw,proc)
            tw[2][1][2].append(zw)
    return tw

def _extractIBMInfo_param(t,tc):
    XPC   ={}; 
    Zones = []
    for z in Internal.getZones(tc):
        allIBCD = Internal.getNodesFromName(z, "IBCD_*")
        for IBCD in allIBCD:
            znames = Internal.getValue(IBCD)
            Zones.append(znames)
            xPC = Internal.getNodesFromName(IBCD,"CoordinateX_PC")[0][1]

            if znames in XPC:
                a = numpy.concatenate((XPC[znames][0],xPC))
                XPC[znames] = [a]
            else:
                XPC[znames] = [xPC]
    Zones = list(set(Zones))
    for zname in Zones:
        xPC = XPC[zname];
        size = xPC[0].shape[0]
        z          = Internal.getNodeFromName ( t, zname)
        o          = Internal.getNodeFromName1( z, '.Solver#ownData')
        param_int  = Internal.getNodeFromName1( o, 'Parameter_int')
        param_real = Internal.getNodeFromName1( o, 'Parameter_real')
        param_int[1][LBM_IBC_NUM] = size
    return None


#=============================================================================
# Compute curvature parameter "K" from geom in 2D
#=============================================================================
def _computeKcurvParameter(tc, tb):

    ################################
    ## EXTRADOS ##
    ################################

    tb_extrados = P.selectCells(tb, '{CoordinateY}>0.')
    x_extrados  = Internal.getNodeFromName(tb_extrados, 'CoordinateX')[1]
    y_extrados  = Internal.getNodeFromName(tb_extrados, 'CoordinateY')[1]

    inds_extrados = x_extrados.argsort()
    x_extrados    = x_extrados[inds_extrados]
    y_extrados    = y_extrados[inds_extrados]

    # evite les divisions par 0
    nPts = len(x_extrados)
    listOfPos = []
    for i in range(1,nPts):
        if x_extrados[i] == x_extrados[i-1]:
            listOfPos.append(i)

    for pos in listOfPos:
        x_extrados = numpy.delete(x_extrados, pos)
        y_extrados = numpy.delete(y_extrados, pos)

    nPts = len(x_extrados)
    y1_extrados = numpy.zeros(nPts)
    y2_extrados = numpy.zeros(nPts)

    # derivee premiere y1
    # schema decentre pour CL
    y1_extrados[0] = (y_extrados[1]-y_extrados[0])/((x_extrados[1]-x_extrados[0]))
    y1_extrados[-1] = (y_extrados[-2]-y_extrados[-1])/((x_extrados[-2]-x_extrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y1_extrados[i] = (y_extrados[i+1]-y_extrados[i-1])/((x_extrados[i+1]-x_extrados[i-1]))

    # derivee seconde y2
    # schema decentre pour CL
    y2_extrados[0] = (y1_extrados[1]-y1_extrados[0])/((x_extrados[1]-x_extrados[0]))
    y2_extrados[-1] = (y1_extrados[-2]-y1_extrados[-1])/((x_extrados[-2]-x_extrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y2_extrados[i] = (y1_extrados[i+1]-y1_extrados[i-1])/((x_extrados[i+1]-x_extrados[i-1]))

    k_extrados = y2_extrados/(numpy.power(1 + numpy.power(y1_extrados, 2), 1.5))
    ka_extrados = numpy.abs(y2_extrados)/(numpy.power(1 + numpy.power(y1_extrados, 2), 1.5))

    ################################
    ## INTRADOS ##
    ################################

    tb_intrados = P.selectCells(tb, '{CoordinateY}<0.')
    x_intrados = Internal.getNodeFromName(tb_intrados, 'CoordinateX')[1]
    y_intrados = numpy.abs(Internal.getNodeFromName(tb_intrados, 'CoordinateY')[1]) # abs pour avoir vraie valeur k

    inds_intrados = x_intrados.argsort()
    x_intrados = x_intrados[inds_intrados]
    y_intrados = y_intrados[inds_intrados]

    # evite les divisions par 0
    nPts = len(x_intrados)
    listOfPos = []
    for i in range(1,nPts):
        if x_intrados[i] == x_intrados[i-1]:
            listOfPos.append(i)

    for pos in listOfPos:
        x_intrados = numpy.delete(x_intrados, pos)
        y_intrados = numpy.delete(y_intrados, pos)

    nPts = len(x_intrados)
    y1_intrados = numpy.zeros(nPts)
    y2_intrados = numpy.zeros(nPts)

    # derivee premiere y1
    # schema decentre pour CL
    y1_intrados[0] = (y_intrados[1]-y_intrados[0])/((x_intrados[1]-x_intrados[0]))
    y1_intrados[-1] = (y_intrados[-2]-y_intrados[-1])/((x_intrados[-2]-x_intrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y1_intrados[i] = (y_intrados[i]-y_intrados[i-1])/((x_intrados[i]-x_intrados[i-1]))

    # derivee seconde y2
    # schema decentre pour CL
    y2_intrados[0] = (y1_intrados[1]-y1_intrados[0])/((x_intrados[1]-x_intrados[0]))
    y2_intrados[-1] = (y1_intrados[-2]-y1_intrados[-1])/((x_intrados[-2]-x_intrados[-1]))

    # schema centre
    for i in range(1, nPts-1):
        y2_intrados[i] = (y1_intrados[i]-y1_intrados[i-1])/((x_intrados[i]-x_intrados[i-1]))

    k_intrados = y2_intrados/(numpy.power(1 + numpy.power(y1_intrados, 2), 1.5))
    ka_intrados = numpy.abs(y2_intrados)/(numpy.power(1 + numpy.power(y1_intrados, 2), 1.5))

    ################################
    ## MaJ tc ##
    ################################
    for z in Internal.getZones(tc):
        for zsr in Internal.getNodesFromType1(z, 'ZoneSubRegion_t'):
            nameSubRegion = zsr[0]
            if nameSubRegion[0:4]=='IBCD':
                ibctypeCR = nameSubRegion.split('_')[1]
                print("ibcTypeCR ", ibctypeCR)
                if ibctypeCR=='100':
                    print("coucou", XOD.__KCURV__)
                    KCurv = Internal.getNodeFromName(zsr, XOD.__KCURV__)[1]
                    coordX = Internal.getNodeFromName(zsr, 'CoordinateX_PW')[1]
                    coordY = Internal.getNodeFromName(zsr, 'CoordinateY_PW')[1]
                    nIBC = numpy.shape(coordX)[0]

                    for i in range(nIBC):
                        if coordY[i] > 0:
                            j = 0
                            while(coordX[i] > x_extrados[j] and j < numpy.shape(x_extrados)[0]-1):
                                j += 1
                            KCurv[i] = k_extrados[j]
                        else:
                            j = 0
                            while(coordX[i] > x_intrados[j] and j < numpy.shape(x_intrados)[0]-1):
                                j += 1
                            KCurv[i] = k_intrados[j]
                        KCurv[i] = min(KCurv[i], 100.)
                        KCurv[i] = max(KCurv[i], -100.)

                    Internal.getNodeFromName(zsr, XOD.__KCURV__)[1] = KCurv

    return None
