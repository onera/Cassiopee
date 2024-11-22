import Generator
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Connector.ToolboxIBM as TIBM
import Connector.PyTree as X
import Transform.PyTree as T
import numpy 

# Generates in parallel a Cartesian mesh
# if ext=0, match and nearmatch joins are not computed
def generateCartMesh(t_case, snears=0.01, dfar=10., dfarList=[], vmin=21, check=False, tbox=None, snearsf=None, 
                     ext=2, dimPb=3, sizeMax=1000000, expand=0):

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones): 
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfarList[c] = Internal.getValue(n)*1.
    # a mettre dans la classe ou en parametre de prepare1 ??? 
    to = None
    fileout = None
    if check: fileout = 'octree.cgns'
    DEPTH= ext
    rank = Cmpi.rank

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')
    # Build octree
    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfars=dfarList, to=to,
                         dimPb=dimPb, vmin=vmin, expand=expand)
    # addRefinementZones
    if tbox is not None:
        o = addRefinementZones(o, tbox, snearsf, vmin=vmin, dim=dimPb)


    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfars=dfarList, to=to, tbox=tbox, 
                                        snearsf=snearsf, dimPb=dimPb, vmin=vmin)
    test.printMem(">>> Octree unstruct [end]")

    # Split octree
    test.printMem(">>> Octree unstruct split [start]")
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[rank]
    del o
    test.printMem(">>> Octree unstruct split [end]")

    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=sizeMax)
    del p 
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%rank
    Cmpi._setProc(t, rank)
    C._addState(t, 'EquationDimension', dimPb)
    test.printMem(">>> Octree struct [end]")

    # Add xzones for ext
    if ext>0:
        test.printMem(">>> extended cart grids [start]")
        tbb = Cmpi.createBBoxTree(t)
        interDict = X.getIntersectingDomains(tbb)
        graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
        del tbb
        Cmpi._addXZones(t, graph, variables=[], cartesian=True)
        test.printMem(">>> extended cart grids [after add XZones]")
        zones = Internal.getZones(t)
        coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
        coords, rinds = Generator.extendCartGrids(coords, ext=DEPTH+1, optimized=1, extBnd=0)
        C.setFields(coords, zones, 'nodes')
        #for noz in range(len(zones)):
        #    Internal.newRind(value=rinds[noz], parent=zones[noz])
        Cmpi._rmXZones(t)
        coords = None; zones = None
        test.printMem(">>> extended cart grids (after rmXZones) [end]")

        TIBM._addBCOverlaps(t, bbox=bb)
        TIBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    else:
        if dimPb == 3: ratios = [[2,2,2],[4,4,4],[8,8,8],[16,16,16]]
        else: ratios = [[2,2,1],[4,4,1],[8,8,1],[16,16,1]]
        tbb = Cmpi.createBBoxTree(t)
        interDict = X.getIntersectingDomains(tbb)
        graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
        del tbb
        Cmpi._addXZones(t, graph, variables=[], cartesian=True)
        test.printMem(">>> extended cart grids [after add XZones]")
        t = X.connectMatch(t, dim=dimPb)
        for ratio0 in ratios:
            t = X.connectNearMatch(t,ratio=ratio0,dim=dimPb)
        Cmpi._rmXZones(t)
        C._fillEmptyBCWith(t,"nref","BCFarfield",dim=dimPb)
    return t

def addRefinementZones(o, tbox, snearsf=None, vmin=15, dim=3):
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
    volmin0 = C.getMinValue(to, 'centers:vol')
    # volume minimum au dela duquel on ne peut pas raffiner
    volmin0 = 1.*volmin0
    while end == 0:
        # Do not refine inside obstacles 
        nob = 0
        C._initVars(to, 'centers:indicator', 0.)
        for box in boxes:
            volmin2 = 1.09*(snearsf[nob]*(vmin-1))**(dim)
            C._initVars(to,'centers:cellN',1.)
            to = TIBM.blankByIBCBodies(to, tbox, 'centers', dim)
            C._initVars(to,'{centers:indicator}=({centers:indicator}>0.)+({centers:indicator}<1.)*logical_and({centers:cellN}<0.001, {centers:vol}>%f)'%volmin2)
            nob += 1

        end = 1
        C.convertPyTree2File(to,'octree.cgns')
        C._initVars(to,'{centers:indicator}={centers:indicator}*({centers:vol}>%g)'%volmin0)
        print("max value", C.getMaxValue(to, 'centers:indicator'))
        C.convertPyTree2File(to,'to.cgns')
        if  C.getMaxValue(to, 'centers:indicator') == 1.: 
            end = 0
            # Maintien du niveau de raffinement le plus fin
            o = Internal.getZones(to)[0]
            o = G.adaptOctree(o, 'centers:indicator', balancing=2)
            to[2][1][2] = [o]
            G._getVolumeMap(to)
            volminloc = C.getMinValue(to, 'centers:vol')
    return Internal.getNodeFromType2(to, 'Zone_t')

#====================================================================================
# Prend les snear dans t, les multiplie par factor
def snearFactor(t, factor=1.):
    tp = Internal.copyRef(t)
    _snearFactor(t, factor)
    return tp

def _snearFactor(t, factor=1.):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, factor*Internal.getValue(n))
    return None

# Set snear in zones
def setSnear(t, value):
    tp = Internal.copyRef(t)
    _setSnear(t, value)
    return tp

def _setSnear(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', value)
    return None

# Set dfar in zones 
def setDfar(t, value):
    tp = Internal.copyRef(t)
    _setDfar(t, value)
    return tp

def _setDfar(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', value)
    return None
