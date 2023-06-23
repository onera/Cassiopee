"""Mesh generation for IBM """
from . import Generator
from . import generator
from . import PyTree as G

import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.IBM as X_IBM
import Post.PyTree as P
import Converter
import Transform
import Converter.GhostCells as CGC
import Connector.PyTree as X
import Converter.Mpi as Cmpi

EPSCART = 1.e-6

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
    X_IBM._addExternalBCs(t, bbox, DEPTH, externalBCType, dimPb)

    nptsTot = 0
    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        nptsTot += niz*njz*nkz
    print('Expected number of points is %d.'%nptsTot)
    return t


def adaptIBMMesh(t, tb, vmin, sensor, factor=1.2, DEPTH=2, sizeMax=4000000,
                 variables=None, refineFinestLevel=False, refineNearBodies=False,
                 check=True, symmetry=0, externalBCType='BCFarfield', fileo='octree.cgns',
                 isAMR=False,valMin=0,valMax=1):
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

def generateIBMMeshPara(tb, vmin=15, snears=None, dimPb=3, dfar=10., dfarList=[], tbox=None,
                    snearsf=None, check=True, symmetry=0, to=None, ext=2,
                    expand=3, dfarDir=0, check_snear=False):    
    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfarList[c] = Internal.getValue(n)*1.

        # refinementSurfFile: surface meshes describing refinement zones
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
        if snearsf is None:
            snearsf = []
            zones = Internal.getZones(tbox)
            for z in zones:
                sn = Internal.getNodeFromName2(z, 'snear')
                if sn is not None: snearsf.append(Internal.getValue(sn))
                else: snearsf.append(1.)
    fileout = None
    if check: fileout = 'octree.cgns'
    # Octree identical on all procs
    if to is not None:
        if isinstance(to, str):
            o = C.convertFile2PyTree(to)
            o = Internal.getZones(o)[0]
        else:
            o = Internal.getZones(to)[0]
        parento = None
    else:
        o = buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList,
                                to=to, tbox=tbox, snearsf=snearsf,
                                dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=Cmpi.rank,
                                expand=expand, dfarDir=dfarDir)

    if Cmpi.rank==0 and check: C.convertPyTree2File(o,fileout)
    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                                        dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=Cmpi.rank)

    if check_snear: exit()
    
    # Split octree
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[Cmpi.rank]
    del o

    # fill vmin + merge in parallel
    res = octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
    del p
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%Cmpi.rank
    Cmpi._setProc(t, Cmpi.rank)

    C._addState(t, 'EquationDimension', dimPb)

    # Add xzones for ext
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)
    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    coords, rinds = Generator.extendCartGrids(coords, ext=ext, optimized=1, extBnd=0)
    C.setFields(coords, zones, 'nodes')
    for noz in range(len(zones)):
        Internal.newRind(value=rinds[noz], parent=zones[noz])
    Cmpi._rmXZones(t)
    coords = None; zones = None
    
    X_IBM._addBCOverlaps(t, bbox=bb)
    X_IBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    if dimPb == 2:
        dz = 0.01
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
               
    return t


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

        elif expand == 2: # expand minimum
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
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
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
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
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
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
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
            to = X_IBM.blankByIBCBodies(to, tb, 'centers', dimPb)
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
        to = X_IBM.blankByIBCBodies(to, tb, 'centers', dim)
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
                print("Info: addRefinementZones: tbox too small - increase tbox by fact = {:2.1f}".format(fact))
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
        X_IBM._addBCOverlaps(zones, bbox0)
    return zones


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






