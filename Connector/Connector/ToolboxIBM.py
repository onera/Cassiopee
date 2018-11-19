"""Toolbox for IBM preprocessing"""
import numpy
import PyTree as X
import OversetData as XOD
import Connector
import connector

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

TypesOfIBC = XOD.TypesOfIBC


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
            for dirl in xrange(dimZone[4]):
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
def adaptIBMMesh(t, tb, vmin, sensor, factor=1.2, DEPTH=2, NP=0, merged=1, sizeMax=4000000,
                 variables=None, refineFinestLevel=False, refineNearBodies=False, 
                 check=True, symmetry=0, externalBCType='BCFarfield', fileo='octree.cgns'):
    
    if fileo is None: raise ValueError("adaptIBMMesh: Octree mesh must be specified by a file.")
    try: to = C.convertFile2PyTree(fileo)
    except: raise ValueError("adaptIBMMesh: %s file not found."%fileo)

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('adaptIBMMesh: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)
    #
    refstate = Internal.getNodeFromName(tb,'ReferenceState')
    #
    if refineNearBodies: constraintSurfaces = []
    else: constraintSurfaces = Internal.getZones(tb)
    if refineFinestLevel: refineLevelF = 1
    else: refineLevelF = 0
    #
    o = Internal.getZones(to)[0]
    dims = Internal.getZoneDim(o)
    npts = dims[1]
    C._initVars(t,"{%s}={%s}*({centers:cellN}>0.)*({centers:cellN}<2.)"%(sensor,sensor))
    C._initVars(to,"centers:indicator", 1.)
    to = P.computeIndicatorValue(to,t,sensor)
    res = P.computeIndicatorField(to,sensor, nbTargetPts=factor*npts, \
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
    if len(res)==3: to = res[0]
    o = Internal.getZones(to)[0]
    o = G.adaptOctree(o, balancing=2)
    C.convertPyTree2File(o,fileo)

    t2 = generateCartMesh__(o, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH, NP=NP, merged=merged, 
                            sizeMax=sizeMax, check=check, symmetry=symmetry, 
                            externalBCType=externalBCType)
    
    #C._initVars(t2,"centers:Density=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))
    #C._initVars(t2,"centers:VelocityX=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityX'))))
    #C._initVars(t2,"centers:VelocityY=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityY'))))
    #C._initVars(t2,"centers:VelocityZ=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'VelocityZ'))))
    #C._initVars(t2,"centers:Temperature=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'Temperature'))))
    #C._initVars(t2,"centers:TurbulentSANuTilde=%g"%(Internal.getValue(Internal.getNodeFromName(refstate,'TurbulentSANuTildeDensity'))/Internal.getValue(Internal.getNodeFromName(refstate,'Density'))))

    # interpolate the solution on the new mesh
    t2 = P.extractMesh(t,t2,3, mode='accurate')
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
    for no in xrange(xmint.shape[0]):
        xmin = xmint[no]; xmax = xmaxt[no]
        ymin = ymint[no]; ymax = ymaxt[no]
        zmin = zmint[no]; zmax = zmaxt[no]
        pool=[]
        for noz in xrange(len(zones)):
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
            res+= T.mergeCart(pool, sizeMax=sizeMax)
            del pool
        elif len(pool)==1: res+=pool
    return res

def octree2StructLoc__(o, vmin=21, ext=0, optimized=0, merged=0, sizeMax=4e6, listOfParents=None):
    sizeMax=int(sizeMax)
    dim = Internal.getZoneDim(o)
    if dim[3] == 'QUAD': dimPb = 2
    elif dim[3] == 'HEXA': dimPb = 3

    if ext == 1: ext = 2

    a = C.getFields(Internal.__GridCoordinates__, o)[0]
    zones = Generator.generator.octree2Struct(a, [vmin])
    c = 1
    for noz in xrange(len(zones)):
        zones[noz] = C.convertArrays2ZoneNode('cartDummy'+str(c), [zones[noz]])        
        c += 1

    if listOfParents is not None:
        bbo = G.bbox(o)
        xmino=bbo[0]; xmaxo=bbo[3]; xmeano=0.5*(xmino+xmaxo)
        ymino=bbo[1]; ymaxo=bbo[4]; ymeano=0.5*(ymino+ymaxo)
        zmino=bbo[2]; zmaxo=bbo[5]; zmeano=0.5*(zmino+zmaxo)
        dx = xmeano-xmino; dy = ymeano-ymino; dz = zmeano-zmino
        eps=1.e-10
        if dimPb == 3: 
            noct = 8
            OCTREEPARENTS=[[],[],[],[],[],[],[],[]]
            ZONES=[[],[],[],[],[],[],[],[]]
            nob = 0
            for zm in [zmino,zmeano]:
                for ym in [ymino,ymeano]:
                    for xm in [xmino,xmeano]:
                        C._initVars(o,'centers:tag',1.)
                        C._initVars(o,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                        C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))
                        C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateZ}>%g)*({centers:CoordinateZ}<%g)'%(zm-eps,zm+dz+eps))
                        o2 = P.selectCells2(o,'centers:tag')
                        if listOfParents is None:
                            OCTREEPARENTS[nob].append(o2)
                        else:
                            for parento in listOfParents:
                                C._initVars(parento,'centers:tag',1.)
                                C._initVars(parento,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                                C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))     
                                C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateZ}>%g)*({centers:CoordinateZ}<%g)'%(zm-eps,zm+dz+eps))
                                parento2 = P.selectCells2(parento,'centers:tag')
                                OCTREEPARENTS[nob].append(parento2)
                        nob+=1
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
                ZONES[noo].append(z)

        else: 
            noct = 4

            OCTREEPARENTS=[[],[],[],[]]
            ZONES=[[],[],[],[]]
            nob = 0

            for ym in [ymino,ymeano]:
                for xm in [xmino,xmeano]:
                    C._initVars(o,'centers:tag',1.)
                    C._initVars(o,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                    C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))
                    o2 = P.selectCells2(o,'centers:tag')

                    if listOfParents is None:
                        OCTREEPARENTS[nob].append(o2)
                    else:
                        for parento in listOfParents:
                            C._initVars(parento,'centers:tag',1.)
                            C._initVars(parento,'{centers:tag}=({centers:CoordinateX}>%g)*({centers:CoordinateX}<%g)'%(xm-eps,xm+dx+eps))
                            C._initVars(parento,'{centers:tag}={centers:tag}*({centers:CoordinateY}>%g)*({centers:CoordinateY}<%g)'%(ym-eps,ym+dy+eps))     
                            parento2 = P.selectCells2(parento,'centers:tag')
                            OCTREEPARENTS[nob].append(parento2)
                    nob+=1
            for z in zones:
                xminz = C.getValue(z,'CoordinateX',0)
                yminz = C.getValue(z,'CoordinateY',0)
                dimZ = Internal.getZoneDim(z)
                ni = dimZ[1]; nj = dimZ[2]; nk = dimZ[3]
                ind = ni-1 + (nj-1)*ni+(nk-1)*ni*nj
                xmaxz = C.getValue(z,'CoordinateX',ind)
                ymaxz = C.getValue(z,'CoordinateY',ind)
                # bbz = G.bbox(z)  
                # xminz=bbz[0]; yminz=bbz[1]
                # xmaxz=bbz[3]; ymaxz=bbz[4]
                if ymaxz < ymeano+eps:
                    if xmaxz < xmeano+eps: noo=0
                    else: noo=1
                else:
                    if xmaxz < xmeano+eps: noo=2
                    else: noo=3

                ZONES[noo].append(z)
        #-----------------------------------------------------------------------------
        #FIN 2D
        #-----------------------------------------------------------------------------
        zones=[]
        for noo in xrange(noct):
            nzones = len(ZONES[noo])
            if nzones > 1:
                nop=0
                print 'Merging %d Cartesian zones of subdomain %d'%(nzones,noo)
                if listOfParents is None:
                    ZONES[noo]=T.mergeCart(ZONES[noo],sizeMax=sizeMax)
                else:
                    for parento in OCTREEPARENTS[noo]:             
                        ZONES[noo]=mergeByParent__(ZONES[noo], parento, sizeMax)
                        nop += 1
                print 'nb of merged zones : %d' %len(ZONES[noo])

        if dimPb==3:
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
    else:
        zones = T.mergeCart(zones,sizeMax=sizeMax)

    print 'After merging: nb Cartesian zones=%d'%(len(zones))

    if ext > 0:
        coords = C.getFields(Internal.__GridCoordinates__,zones)
        coords = Generator.extendOctreeGrids__(coords,ext=ext, optimized=optimized)
        C.setFields(coords, zones, 'nodes')
    # Creation des zones du pyTree
    for z in zones: z[0] = C.getZoneName('cart')
    if ext==0:
        if dimPb == 3: ratios = [[2,2,2],[4,4,4],[8,8,8],[16,16,16]]
        else: ratios = [[2,2,1],[4,4,1],[8,8,1],[16,16,1]]
        zones = X.connectMatch(zones, dim=dimPb)
        for ratio0 in ratios:
            zones = X.connectNearMatch(zones,ratio=ratio0,dim=dimPb)
        return zones
    else:
        bbox0 = G.bbox(o)
        xmin = bbox0[0]; ymin = bbox0[1]; zmin = bbox0[2]
        xmax = bbox0[3]; ymax = bbox0[4]; zmax = bbox0[5]
        noz = 0
        for z in zones:
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
            if (x1 > xmin+EPSCART): C._addBC2Zone(z,'overlap1','BCOverlap','imin')
            if (x2 < xmax-EPSCART): C._addBC2Zone(z,'overlap2','BCOverlap','imax')
            if (y1 > ymin+EPSCART): C._addBC2Zone(z,'overlap3','BCOverlap','jmin')
            if (y2 < ymax-EPSCART): C._addBC2Zone(z,'overlap4','BCOverlap','jmax')
            if (z1 > zmin+EPSCART): C._addBC2Zone(z,'overlap5','BCOverlap','kmin')
            if (z2 < zmax-EPSCART): C._addBC2Zone(z,'overlap6','BCOverlap','kmax')
            zones[noz] = z; noz += 1
    return zones

def generateCartMesh__(o, dimPb=3, vmin=11, DEPTH=2, NP=0, merged=1, sizeMax=4000000, check=True, 
                       symmetry=0, externalBCType='BCFarfield',mergeByParents=True):

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
    if DEPTH==0: ext=0
    else: ext=DEPTH+1

    dimZO = Internal.getZoneDim(o)[2]
    listOfParentOctrees=None
    if mergeByParents and nzones0>1000:
        listOfParentOctrees=[]
        parento = C.initVars(o,'centers:indicator',-1.)
        parento = G.adaptOctree(parento)
        
        for nop in xrange(3):
            C._initVars(parento,'centers:indicator',-1.)
            parento = G.adaptOctree(parento)
            listOfParentOctrees.append(parento)

    res = octree2StructLoc__(o, vmin=vmin, ext=ext, optimized=optimized, merged=merged, sizeMax=sizeMax, 
                             listOfParents=listOfParentOctrees)
    t = C.newPyTree(['CARTESIAN', res])
    dz = 0.01
    if dimPb == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
    bbo = G.bbox(o)        
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
        #bbz = G.bbox(zp)
        external = False
        for idir in dirs:
            if abs(bbz[idir]-bbo[idir])< 1.e-6:                    
                C._addBC2Zone(zp, 'external', externalBCType, rangeDir[idir])
                external = True
        if externalBCType != 'BCOverlap' and externalBCType != 'BCDummy':
            if external: _modifPhysicalBCs__(zp, depth=DEPTH, dimPb=dimPb)
    print 'Expected number of points is %d'%nptsTot
    return t

#--------------------------------------------------------------------------
# to : maillage octree, si not None : on le prend comme squelette 
#--------------------------------------------------------------------------
def generateIBMMesh(tb, vmin=15, snears=None, dfar=10., DEPTH=2, NP=0, tbox=None, 
                    snearsf=None, check=True, merged=1, sizeMax=4000000, 
                    symmetry=0, externalBCType='BCFarfield', to=None, 
                    composite=0, mergeByParents=True, fileo='octree.cgns'):
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    if dimPb is None: raise ValueError('generateIBMMesh: EquationDimension is missing in input body tree.')
    dimPb = Internal.getValue(dimPb)
    
    # type de traitement paroi: pts interieurs ou externes
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('generateIBMMesh: GoverningEquations is missing in input body tree.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)
    if model == 'Euler': IBCType =-1
    else: IBCType = 1 #Turbulent

    # Build octree
    i = 0; surfaces=[]; snearso=[] # pas d'espace sur l'octree
    bodies = Internal.getZones(tb)
    if not isinstance(snears, list): snears = len(bodies)*[snears]
    if len(bodies) != len(snears):
        raise ValueError('generateIBMMesh: Number of bodies is not equal to the size of snears.')
    dxmin0 = 1.e10
    for s in bodies:
        sdd = Internal.getNodeFromName1(s, ".Solver#define")
        if sdd is not None:
            snearl = Internal.getNodeFromName1(sdd, "snear")
            if snearl is not None: 
                snearl = Internal.getValue(snearl)
                snears[i] = snearl
        dhloc = snears[i]*(vmin-1)
        surfaces+=[s]; snearso+=[dhloc]
        dxmin0 = min(dxmin0,dhloc)
        i += 1

    if to is None:
        o = G.octree(surfaces, snearso, dfar=dfar, balancing=2)
        G._getVolumeMap(o); volmin = C.getMinValue(o, 'centers:vol')
        dxmin = (volmin)**(1./dimPb)
        if dxmin < 0.65*dxmin0: 
            snearso = [2.*i for i in snearso]
            o = G.octree(surfaces, snearso, dfar=dfar, balancing=2)

        if composite: # maillages proche corps -> maillage cubique initial decoupe
            tb_bbox = G.bbox(tb)
            xminb = tb_bbox[0];  xmaxb = tb_bbox[3]
            yminb = tb_bbox[1];  ymaxb = tb_bbox[4]
            zminb = tb_bbox[2];  zmaxb = tb_bbox[5]
            bbo = G.bbox(o)
            xmino = bbo[0];  xmaxo = bbo[3]
            ymino = bbo[1];  ymaxo = bbo[4]
            zmino = bbo[2];  zmaxo = bbo[5]
            dfarmin = xmaxo-xmaxb; dirmin = 1
            dfary = ymaxo-ymaxb
            if dfary < dfarmin:
                dfarmin = dfary; dirmin = 2
            if dimPb == 2:
                if dirmin == 1: 
                    C._initVars(o,'{centers:tag}=({centers:CoordinateY}<%g)*({centers:CoordinateY}>%g)'%(ymaxb+dfarmin,yminb-dfarmin))
                elif dirmin == 2:
                    C._initVars(o,'{centers:tag}=({centers:CoordinateX}<%g)*({centers:CoordinateX}>%g)'%(xmaxb+dfarmin,xminb-dfarmin))

            elif dimPb == 3:
                dfarz = zmaxo-zmaxb
                if dfarz < dfarmin:
                    dfarmin = dfarz; dirmin = 3
                C._initVars(o,"centers:tag", 1.)
                if dirmin != 1: 
                    C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateX}<%g)*({centers:CoordinateX}>%g)'%(xmaxb+dfarmin,xminb-dfarmin))
                if dirmin != 2:
                    C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateY}<%g)*({centers:CoordinateY}>%g)'%(ymaxb+dfarmin,yminb-dfarmin))
                if dirmin != 3:
                    C._initVars(o,'{centers:tag}={centers:tag}*({centers:CoordinateZ}<%g)*({centers:CoordinateZ}>%g)'%(zmaxb+dfarmin,zminb-dfarmin))
            
            o = P.selectCells2(o,'centers:tag')
            
        symmetry = 0
        if symmetry != 0:
            bb = G.bbox(o)
            xmoy = 0.5*(bb[3]+bb[0])
            ymoy = 0.5*(bb[4]+bb[1])
            zmoy = 0.5*(bb[5]+bb[2])
            if   symmetry== 1: o = P.selectCells(o,'{CoordinateX}>%g'%(xmoy-TOLDIST))
            elif symmetry==-1: o = P.selectCells(o,'{CoordinateX}<%g'%(xmoy+TOLDIST))
            elif symmetry== 2: o = P.selectCells(o,'{CoordinateY}>%g'%(ymoy-TOLDIST))
            elif symmetry==-2: o = P.selectCells(o,'{CoordinateY}<%g'%(ymoy+TOLDIST))
            elif symmetry== 3: o = P.selectCells(o,'{CoordinateZ}>%g'%(zmoy-TOLDIST))
            elif symmetry==-3: o = P.selectCells(o,'{CoordinateZ}<%g'%(zmoy+TOLDIST))

        vmint = 31
        if vmin < vmint:
            print 'generateIBMMesh: octree finest level expanded (expandLayer activated).'
            to = C.newPyTree(['Base',o])
            to = blankByIBCBodies(to, tb, 'centers', dimPb)
            C._initVars(o,"centers:indicator", 0.)
            cellN = C.getField("centers:cellN",to)[0]
            octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
            indic = C.getField("centers:indicator",o)[0]
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,0,0)
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,1,0) # CB
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,2,0) # CB
            indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,3,0) # CB
                                                                                      
            indic = Converter.addVars([indic,cellN])
            indic = Converter.initVars(indic,"{indicator}={indicator}*({cellN}>0.)")
            octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
            o = C.convertArrays2ZoneNode(o[0],[octreeA])

        to = C.newPyTree(['Base', o])
        G._getVolumeMap(to); volmin = C.getMinValue(to, 'centers:vol')

        dxmin = (volmin)**(1./dimPb)
        print 'Minimum spacing of Cartesian mesh= %f (targeted %f)'%(dxmin/(vmin-1),dxmin0/(vmin-1))

        if tbox is not None and snearsf is not None:
            o = addRefinementZones(o, tb, tbox, snearsf, vmin, dimPb)
        if IBCType == 1 and dimPb == 3: # external points
            to = C.newPyTree('Base'); to[2][1][2]=[o]
            to = blankByIBCBodies(to, tb, 'nodes', dimPb)
            to = X.setHoleInterpolatedPoints(to,depth=-1,loc='nodes')
            # to = P.selectCells2(to,"cellN",strict=0)
            o = Internal.getZones(to)[0]

        nelts = Internal.getZoneDim(o)[2] 
        if nelts > 20000 and merged == 1: 
            print 'Warning: number of zones (%d) might be too big (block merging might last a long time). Try to increase vmin or deactivate merging.'%nelts

    else:
        o = Internal.getZones(to)[0]

    if fileo is not None: C.convertPyTree2File(o,fileo)
    res = generateCartMesh__(o, dimPb=dimPb, vmin=vmin, DEPTH=DEPTH, NP=NP, merged=merged, sizeMax=sizeMax, 
                             check=check, symmetry=symmetry, externalBCType=externalBCType, mergeByParents=mergeByParents)
    return res

def _removeBlankedGrids(t,loc='centers'):
    vari = 'cellNIBC'
    varc = 'cellNChim'
    flag='flag'
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
        
    if not isinstance(snearsf, list): snearsf = len(boxes)*[snearsf]
    if len(boxes) != len(snearsf):
        raise ValueError('addRefinementZones: Number of refinement bodies is not equal to the length of snearsf list.')
    to = C.newPyTree(['Base', o])
    BM = numpy.ones((1,1),numpy.int32)
    end = 0
    G._getVolumeMap(to)
    volmin0 = C.getMinValue(to,'centers:vol')
    #print 'volmin0 = ', volmin0, (volmin0)**(1./3.)
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
            volmin2 = 1.09*(snearsf[nob]*(vmin-1))**(dim)
            C._initVars(to,'centers:cellN',1.)
            to = X.blankCells(to, [box], BM, blankingType='center_in', dim=dim, delta=1.e-10, tol=1.e-8)
            C._initVars(to,'{centers:indicator}=({centers:indicator}>0.)+({centers:indicator}<1.)*logical_and({centers:cellN}<0.001, {centers:vol}>%f)'%volmin2)
            nob += 1

        end = 1
        C._initVars(to,'{centers:indicator}={centers:indicator}*({centers:cellNBody}>0.)*({centers:vol}>%g)'%volmin0)

        if  C.getMaxValue(to,'centers:indicator') == 1.: 
            end = 0
            # Maintien du niveau de raffinement le plus fin
            o = Internal.getZones(to)[0]
            o = G.adaptOctree(o, 'centers:indicator', balancing=2)
            to[2][1][2] = [o]
            G._getVolumeMap(to)
            volminloc = C.getMinValue(to,'centers:vol')
            #print 'volminloc = ', volminloc, (volminloc)**(1./3.)
    return Internal.getNodeFromType2(to, 'Zone_t')

# =============================================================================
# Calcul des points IBM a corriger, paroi et a interpoler
# =============================================================================
def getAllIBMPoints(t, loc='nodes', hi=0., he=0., tb=None, tfront=None, 
                    frontType=0, cellNName='cellN', IBCType=1, depth=2):

    allCorrectedPts = []; allWallPts = []; allInterpPts = []
    #-------------------------------------------
    # 1. Get the list of IBC corrected pts
    #-------------------------------------------
    listOfSnearsLoc=[]
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
    #-------------------------------------------
    # 2. Get the list of IBC wall and interp pts
    #-------------------------------------------        
    indcell = Converter.extractVars(allCorrectedPts,['indcell'])
    if tb is None or tfront is None: # constant hi, he   
        for nozc in xrange(len(allCorrectedPts)):
            poshe = KCore.isNamePresent(allCorrectedPts[nozc],'he')
            if poshe == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'he',he)
            poshi = KCore.isNamePresent(allCorrectedPts[nozc],'hi')
            if poshi == -1: allCorrectedPts[nozc] = Converter.initVars(allCorrectedPts[nozc],'hi',hi)
        res = connector.getIBMPtsBasic(allCorrectedPts, varsn, 'TurbulentDistance')
    else:
        dictOfBodiesByIBCType={}
        ibctypemin = 100
        ibctypemax = -1 # pour les verifications
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
            
            ibctypeI = TypesOfIBC[ibctype]
            ibctypemin = min(ibctypemin,ibctypeI)
            if ibctypeI < 4: ibctypemax = max(ibctypemax,ibctypeI)

            if not dictOfBodiesByIBCType.has_key(ibctype):
                dictOfBodiesByIBCType[ibctype]=[s]
            else: dictOfBodiesByIBCType[ibctype]+=[s]
        
        if ibctypemin < 2 and ibctypemax > 1:
            raise ValueError("prepareIBMData: cannot mix slip/noslip IBC with wall-model IBCs within the same simulation.")

        # Regroupement des corps par type de BC - optimise les projections ensuite 
        bodies = []; listOfIBCTypes=[];
        for itype in dictOfBodiesByIBCType.keys():
            s = dictOfBodiesByIBCType.get(itype)
            body = C.getFields(Internal.__GridCoordinates__,s)
            body = Converter.convertArray2Tetra(body)
            body = Transform.join(body)
            bodies.append(body)
            listOfIBCTypes.append(TypesOfIBC[itype])
        if frontType == 0:
            dmin = C.getMaxValue(tfront, 'TurbulentDistance')
            allCorrectedPts = Converter.initVars(allCorrectedPts,'dist',dmin)
            res = connector.getIBMPtsWithoutFront(allCorrectedPts, bodies, listOfIBCTypes, varsn, 'dist', IBCType)
        else:            
            front = C.getFields(Internal.__GridCoordinates__,tfront)
            front = Converter.convertArray2Tetra(front)
            allCorrectedPts = Converter.extractVars(allCorrectedPts,['CoordinateX','CoordinateY','CoordinateZ']+varsn)
            res = connector.getIBMPtsWithFront(allCorrectedPts, listOfSnearsLoc, bodies, listOfIBCTypes,
                                               front, varsn, IBCType, depth)

    Converter.convertArrays2File(allCorrectedPts, 'correctedPts.plt')
    Converter.convertArrays2File(res[0], 'wallPts.plt')
    Converter.convertArrays2File(res[1], 'imagePts.plt')
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
        for noz in xrange(nzonesR):
            indicesByTypeForZone = res[2][noz]
            nbTypes = len(indicesByTypeForZone)        
            for nob in xrange(nbTypes):
                ibcTypeL = listOfIBCTypes[nob]
                indicesByTypeL = indicesByTypeForZone[nob]
                if indicesByTypeL.shape[0] > 0 :                
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
        for noz in xrange(nzonesR):
            if noz == 0:
                dictOfCorrectedPtsByIBCType[ibcTypeL] = [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] = [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] = [allInterpPts[noz]]
            else:
                dictOfCorrectedPtsByIBCType[ibcTypeL] += [allCorrectedPts[noz]]
                dictOfWallPtsByIBCType[ibcTypeL] += [allWallPts[noz]]
                dictOfInterpPtsByIBCType[ibcTypeL] += [allInterpPts[noz]]        
    # for ibcTypeL in  dictOfCorrectedPtsByIBCType.keys():
    #     Converter.convertArrays2File(dictOfCorrectedPtsByIBCType[ibcTypeL],"correctedPts_%d.plt"%ibcTypeL)
    #     Converter.convertArrays2File(dictOfWallPtsByIBCType[ibcTypeL],"wallPts_%d.plt"%ibcTypeL)
    #     Converter.convertArrays2File(dictOfInterpPtsByIBCType[ibcTypeL],"interpPts_%d.plt"%ibcTypeL)

    #return allCorrectedPts, allWallPts, allInterpPts
    return dictOfCorrectedPtsByIBCType, dictOfWallPtsByIBCType, dictOfInterpPtsByIBCType
#=============================================================================
# Returns the front defining the image points
#=============================================================================
def getIBMFront(tc, frontvar, dim, frontType, shiftIBM):

    if frontType == 1: front = getIBMFrontType1(tc,frontvar,dim)
    else: front = getIBMFrontType0(tc,frontvar,dim)
    front = C.deleteEmptyZones(front)
    Internal._rmNodesFromName(front,"ID_*")
    if shiftIBM is False: return front
    
    dxmin = 1e12
    if frontType==1:
        front = Internal.getZones(front)
        dxmax = 1.e-12
        dht = [[]]*len(front)
        nof = 0
        for f in front:
            subf = T.subzone(f, [0], type='elements')
            dx = C.getMaxValue(subf,"CoordinateX")-C.getMinValue(subf,"CoordinateX")
            if dim == 3:
                dy = C.getMaxValue(subf,"CoordinateY")-C.getMinValue(subf,"CoordinateY")
                dz = C.getMaxValue(subf,"CoordinateZ")-C.getMinValue(subf,"CoordinateZ")
                dht[nof] = max(dx,dy,dz)
            else: dht[nof] = dx

            if dht[nof] < dxmin and dht[nof] > 1.e-12: 
                dxmin = dht[nof]
            if dht[nof] > dxmax: 
                dxmax = dht[nof]
            nof+=1

        nlevels = int(math.log(dxmax/dxmin)/math.log(2)+1)

        dictOfLevels={}
        for l in xrange(nlevels): dictOfLevels[l]=[]
        for nof in xrange(len(dht)):
            if dht[nof]>1.e-12:
                nolvl = int(math.log(dht[nof]/dxmin)/math.log(2))
                dictOfLevels[nolvl].append(front[nof])
        
        front=[]
        for nol in xrange(nlevels):
            if len(dictOfLevels[nol])>0:
                front.append(T.join(dictOfLevels[nol]))

    if shiftIBM and frontType==1:
        # C.convertPyTree2File(front,"front0.cgns")
        G._getNormalMap(front)
        C._normalize(front, ['centers:sx','centers:sy','centers:sz'])
        vshift = min(1.e-6,0.01*dxmin)
        C._initVars(front,"{centers:sx}={centers:sx}*%g"%vshift)
        C._initVars(front,"{centers:sy}={centers:sy}*%g"%vshift)
        C._initVars(front,"{centers:sz}={centers:sz}*%g"%vshift)
        front = C.center2Node(front, Internal.__FlowSolutionCenters__)
        Internal._rmNodesByName(front,Internal.__FlowSolutionCenters__)
        front2 = G.grow(front, ['sx','sy','sz'])
        C._initVars(front2,"tag",2.)
        C._initVars(front,"tag",1.)
        hook = C.createGlobalHook(front, 'nodes')
        front2 = C.identifySolutions(front2, front, hookN=hook, tol=1.e-14)
        C._initVars(front2,'{tag}=abs({tag}-1.)')
        front2 = P.exteriorFaces(front2)
        front = P.selectCells2(front2,'tag',strict=1)
        # C.convertPyTree2File(front,"front1.cgns")
        Internal._rmNodesByName(front,Internal.__FlowSolutionNodes__)
        C.freeHook(hook)
    return front

# front of first computed cells - with overlapping
def getIBMFrontType1(tc,frontvar, dim):
    if dim == 2:
        z0 = Internal.getNodeFromType2(tc, 'Zone_t')
        zmean = C.getValue(z0, 'CoordinateZ', 0)
        dz = 2*zmean
    else: dz = 0.
    front = []
    for z in Internal.getZones(tc):
        if C.getMinValue(z,frontvar)==0. and C.getMaxValue(z,frontvar)==1.:
            X._maximizeBlankedCells(z,depth=1,dir=1,loc='nodes', cellNName='cellNChim')
            C._initVars(z,'{cellNChim}=minimum(1.,{cellNChim})')
            f = P.frontFaces(z, frontvar)
            if Internal.getZoneDim(f)[1]>0:  
                Internal._rmNodesByName(f,'ID_*')
                front.append(f)

    C._initVars(front,'{tag}=({cellNChim}>0.5)*({cellNChim}<1.5)')
    front = P.selectCells2(front,'tag',strict=1)
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

    print 'Minimum spacing on Cartesian grids = ',dxmin
    return dxmin

#==============================================================================
# masquage par les corps IBC
# gridType = single or composite - composite means that an off body grid exists
#==============================================================================
def blankByIBCBodies(t, tb, loc, dim, gridType='single'):
    if gridType == 'composite': blankalgo='xray'
    else: blankalgo='tri'
    DIM = dim
    if DIM == 2: blankalgo = 'xray'
    bodies = []
    for b in Internal.getBases(tb):
        wallsl = Internal.getNodesFromType1(b, 'Zone_t')
        if wallsl != []:
            try:
                wallsl = C.convertArray2Tetra(wallsl)
                wallsl = T.join(wallsl)
                wallsl = G.close(wallsl)
                if DIM == 3:
                    try: P.exteriorFaces(wallsl)
                    except: bodies.append([wallsl])
                else: bodies.append([wallsl])
            except: bodies.append(wallsl)
            wallsl = C.convertArray2Tetra(wallsl)
            #wallsl = T.join(wallsl)
            #wallsl = G.close(wallsl)
            #P.exteriorFaces(wallsl)

    nbodies = len(bodies)
    print 'Blanking mesh by %d immersed bodies'%nbodies
    if loc == 'centers': typeb = 'center_in'
    else: typeb = 'node_in'
    nbases = len(Internal.getBases(t))
    BM = numpy.ones((nbases,nbodies),dtype=numpy.int32)
    if gridType=='composite':
        if nbodies < nbases:
           for nob in xrange(nbodies): BM[nob,nob]=0
    if blankalgo == 'xray' or DIM == 2:
        dh_min = getMinimumCartesianSpacing(t)
        XRAYDIM1 = 2000; XRAYDIM2 = XRAYDIM1
        if dh_min > 0.:
            bb = G.bbox(tb)
            Lxref = bb[3]-bb[0]
            Lyref = bb[4]-bb[1]
            XRAYDIM1 = max(XRAYDIM1,int(Lxref/(0.15*dh_min)))
            XRAYDIM2 = max(XRAYDIM2,int(Lyref/(0.15*dh_min)))
        if DIM == 2:  XRAYDIM2 = 2

        #print 'XRAYDIM=', XRAYDIM1, XRAYDIM2
        if loc == 'centers':
            tc = C.node2Center(t)
            tc = X.blankCells(tc, bodies, BM,blankingType='node_in',XRaydim1=XRAYDIM1,XRaydim2=XRAYDIM2,dim=DIM)
            C._cpVars(tc,'cellN',t,'centers:cellN')
        else:
            t = X.blankCells(t, bodies,BM,blankingType=typeb,delta=TOLDIST,XRaydim1=XRAYDIM1,XRaydim2=XRAYDIM2,dim=DIM)
    else:
        t = X.blankCellsTri(t,bodies,BM,blankingType=typeb)
    return t

#=============================================================================
# distance signee en fonction du masquage Chimere et IBC
#=============================================================================
def _signDistance(t):
    C._initVars(t,'{centers:TurbulentDistance}=-1.*({centers:cellNIBC}*{centers:cellNChim}<1.)*{centers:TurbulentDistance}+({centers:cellNIBC}*{centers:cellNChim}>0.)*{centers:TurbulentDistance}')
    return None

#=============================================================================
def doInterp(t, tc, tbb, shiftIBM, tb=None, typeI='ID', dim=3, dictOfADT=None, frontType=0, depth=2, IBCType=1):    
    ReferenceState = Internal.getNodeFromType2(t,'ReferenceState_t')

    if typeI == 'ID':
        # toutes les zones sont interpolables en Chimere
        intersectionsDict = X.getIntersectingDomains(tbb, method='AABB', taabb=tbb)

        rcvZones = []
        for zrcv in Internal.getZones(t):
            if C.getMaxValue(zrcv,'centers:cellN')==2.:
                zrcvname = zrcv[0]; rcvZones.append(zrcv)
        nozr = 0
        nbZonesChim = len(rcvZones)
        for nozr in xrange(nbZonesChim):
            zrcv = rcvZones[nozr]
            zrcvname = zrcv[0]
            nozr += 1; hook0 = []
            nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]
            for nobd in xrange(len(tc[2])):
                if tc[2][nobd][3] == 'CGNSBase_t':
                    for nozd in xrange(len(tc[2][nobd][2])):
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
            for nod in xrange(len(dnrZones)):
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

        print 'Building the IBM front.'  
        front = getIBMFront(tc, 'cellNFront', dim, frontType, shiftIBM)
        C.convertPyTree2File(front,"front.cgns")
        #
        res = getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType, \
                              cellNName='cellNIBC', depth=depth, IBCType=IBCType)
        nbZonesIBC = len(zonesRIBC)
        dictOfADT = {}
        dictOfCorrectedPtsByIBCType = res[0]
        dictOfWallPtsByIBCType = res[1] 
        dictOfInterpPtsByIBCType = res[2]
        for ibcTypeL in  dictOfCorrectedPtsByIBCType.keys():
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in xrange(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    interpPtsBB=Generator.BB(allInterpPts[nozr])

                    zrcv = zonesRIBC[nozr]
                    zrcvname = zrcv[0]
                    nobOfDnrBases = []; nobOfDnrZones=[]; dnrZones=[]; hook0 = []
                    for nobd in xrange(len(tc[2])):
                        if tc[2][nobd][3] == 'CGNSBase_t':
                            for nozd in xrange(len(tc[2][nobd][2])):
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
                                             nature=1, penalty=1, loc='centers', storage='inverse',  hook=hook0, dim=dim, \
                                             ReferenceState=ReferenceState, bcType=ibcTypeL)
                    nozr += 1
                    for nod in xrange(len(dnrZones)):
                        nobd = nobOfDnrBases[nod]
                        nozd = nobOfDnrZones[nod]
                        tc[2][nobd][2][nozd] = dnrZones[nod]
        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])
    return tc

#=============================================================================
# Performs the full IBM preprocessing using overlapping Cartesian grids
#=============================================================================
def prepareIBMData(t, tbody, DEPTH=2, loc='centers', frontType=1, shiftIBM=False, interp='all', inv=False):
    tb =  Internal.copyRef(tbody)
    if interp == 'all': interpI = 0
    elif interp=='chimera': interpI = 1
    elif interp=='ibm': interpI = 2
    elif interp=='composite': interpI=3# on garde les coordonnees dans tc
    
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
    t = X.applyBCOverlaps(t, depth=DEPTH)
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
        print 'computing distance field...'
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
        if not shiftIBM:
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        else:
            DEPTHL=DEPTH+1
            X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0, loc='centers',cellNName='cellN',addGC=False)         
            _blankClosestTargetCells(t,cellNName='cellN', depth=DEPTHL)
    else:
        raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    _removeBlankedGrids(t, loc='centers')
    print 'Nb of Cartesian grids=', len(Internal.getZones(t))
    npts = 0
    for i in Internal.getZones(t):
       dims = Internal.getZoneDim(i)
       npts += dims[1]*dims[2]*dims[3]
    print 'Final number of points=%5.4f millions.'%(npts/1000000.)

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
        print 'Euler : on repousse le front un peu plus loin.'
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
    Internal._rmNodesByName(tc, "grad*TurbulentDistance")
    Internal._rmNodesByName(tc, Internal.__FlowSolutionCenters__)

    tbb = G.BB(tc)

    if interpI < 2:
        #Creation du dictionnaire des ADT
        #En chimere toutes les zones sont interpolables pour l instant
        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
        print 'Interpolations Chimere'
        tc = doInterp(t, tc, tbb, shiftIBM, tb=None, typeI='ID', dim=dimPb, 
                      dictOfADT=dictOfADT)
        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])

    elif interpI == 3: # maillage mobile : on n interpole pas les blocs externes ici
        for z in Internal.getZones(t):
            if Internal.getNodesFromName(z,"ov_ext*")!=[]:
                C._initVars(z,'centers:cellN', 1.)

        dictOfADT = {}
        for zdnr in Internal.getZones(tc):
            zdnrname = zdnr[0]
            if zdnrname not in dictOfADT:
                HOOKADT = C.createHook(zdnr, 'adt')
                dictOfADT[zdnrname] = HOOKADT
        print 'Interpolations Chimere'
        tc = doInterp(t, tc, tbb, shiftIBM, tb=None, typeI='ID', dim=dimPb, 
                      dictOfADT=dictOfADT)
        for dnrname in dictOfADT.keys(): C.freeHook(dictOfADT[dnrname])

    # setIBCData - IBC
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')    
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    #-----------------------------------------------
    # Transfert du cellNFront
    # A OPTIMISER : NE FAIRE LE TRANSFERT QUE POUR LES ZONES DU FRONT INTERSECTANTES 
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    XOD._setInterpTransfers(t,tc,variables=['cellNFront'],cellNVariable='cellNFront',compact=0)

    ## Fin traitement specifique, vaut 0 ou 1 apres la ligne suivante
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    if interpI != 1:
        print 'Minimum distance: ', C.getMinValue(t,'centers:TurbulentDistance')
        t = P.computeGrad2(t,'centers:TurbulentDistance')
        print 'Interpolations IBM'
        tc = doInterp(t,tc,tbb, shiftIBM, tb=tb,typeI='IBCD',dim=dimPb, dictOfADT=None, frontType=frontType, depth=DEPTH, IBCType=IBCType)

    # cleaning...
    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)

    if interpI==0: Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
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
# Extraction des infos pour le post traitement
# if td=None: return the cloud of points
# else interpolate on td
#=============================================================================
def extractIBMWallFields(tc,tb=None):
    xwNP = []; ywNP = []; zwNP = []
    pressNP = []; utauNP = []; yplusNP = []; densNP = []
    for z in Internal.getZones(tc):
        allZSR = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        if allZSR != []:
            allIBCD = Internal.getNodesFromName(allZSR,"IBCD_*")
            for IBCD in allIBCD:
                xPW = Internal.getNodeFromName1(IBCD,"CoordinateX_PW")[1]
                yPW = Internal.getNodeFromName1(IBCD,"CoordinateY_PW")[1]
                zPW = Internal.getNodeFromName1(IBCD,"CoordinateZ_PW")[1]
                xwNP.append(xPW); ywNP.append(yPW); zwNP.append(zPW)

                PW = Internal.getNodeFromName1(IBCD,X.__PRESSURE__)
                if PW is not None: pressNP.append(PW[1])
                RHOW = Internal.getNodeFromName1(IBCD,X.__DENSITY__)
                if RHOW is not None: densNP.append(RHOW[1])
                UTAUW = Internal.getNodeFromName1(IBCD,X.__UTAU__)
                if UTAUW is not None: utauNP.append(UTAUW[1])

                YPLUSW = Internal.getNodeFromName1(IBCD, X.__YPLUS__)
                if YPLUSW is not None: yplusNP.append(YPLUSW[1])

    if pressNP == []: return None
    else:
        pressNP = numpy.concatenate(pressNP)
        densNP = numpy.concatenate(densNP)
        if utauNP != []: utauNP = numpy.concatenate(utauNP)
        if yplusNP != []: yplusNP = numpy.concatenate(yplusNP)
        xwNP = numpy.concatenate(xwNP)
        ywNP = numpy.concatenate(ywNP)
        zwNP = numpy.concatenate(zwNP)

    # Creation d une seule zone
    zsize = numpy.empty((1,3), numpy.int32, order='Fortran')
    zsize[0,0] = xwNP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=z)
    coordx = ['CoordinateX',xwNP,[],'DataArray_t']
    coordy = ['CoordinateY',ywNP,[],'DataArray_t']
    coordz = ['CoordinateZ',zwNP,[],'DataArray_t']
    gc[2] = [coordx,coordy,coordz]
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                   gridLocation='Vertex', parent=z)
    FSN[2].append([X.__PRESSURE__,pressNP, [],'DataArray_t'])
    FSN[2].append([X.__DENSITY__,densNP, [],'DataArray_t'])
    utauPresent = 0; yplusPresent = 0
    if utauNP != []:
        utauPresent = 1
        FSN[2].append([X.__UTAU__,utauNP, [],'DataArray_t'])
    if yplusNP != []:
        yplusPresent = 1
        FSN[2].append([X.__YPLUS__,yplusNP, [],'DataArray_t'])
    if tb is None: return z
    else:
        dimPb = Internal.getNodeFromName(tb,'EquationDimension')
        if dimPb is None: 
            print 'Warning: extractIBMWallFields: pb dimension is set to 3.'
            dimPb = 3
        else:
            dimPb = Internal.getValue(dimPb)
        td = Internal.copyRef(tb)
        for nob in xrange(len(td[2])):
            b = td[2][nob]
            if b[3] == 'CGNSBase_t':                
                zones = Internal.getNodesFromType1(b, 'Zone_t')
                zones = C.convertArray2Tetra(zones)
                zones = T.join(zones); zones = G.close(zones)
                b[2] = [zones]
        C._initVars(td,X.__PRESSURE__,0.)
        C._initVars(td,X.__DENSITY__,0.)
        if utauPresent==1: C._initVars(td,X.__UTAU__,0.)
        if yplusPresent==1: C._initVars(td,X.__YPLUS__,0.)
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
            znamea = IBCD[1]
            znames = znamea.tostring()
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
            if znames in XPW.keys():
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
