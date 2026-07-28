# Check CAD
import OCC.PyTree as OCC

# Return dict of face->edges (topologic)
def getEdgesByFace(hook):
    out = {}
    nbFaces = OCC.getNbFaces(hook)
    for i in range(nbFaces):
        edgeno = OCC.getEdgesByFace(hook, i+1)
        out[i+1] = edgeno
    return out

# Return dict of edge->faces (topologic)
def getFacesByEdge(hook, edgesByFace):
    out = {}
    for i in edgesByFace:
        for j in edgesByFace[i]:
            if j not in out: out[j] = [i]
            else:
                out[j].append(i)
    return out

# Return the list of degenerated edges (geometric)
def getDegeneratedEdges(hook, tol=1.e-12):
    out = []
    nbEdges = OCC.getNbEdges(hook)
    for i in range(nbEdges):
        L = OCC.getEdgeLength(hook, [i+1])
        if L <= tol: out.append(i+1)
    return out

# Return the list of degenerated faces (geometric)
def getDegeneratedFaces(hook, tol=1.e-12):
    out = []
    nbFaces = OCC.getNbFaces(hook)
    for i in range(nbFaces):
        area = OCC.getFaceArea(hook, [i+1])
        if area <= tol: out.append(i+1)
    return out

# Check if faces overlap
def getFaceOverlap(hook, tol=1.e-12, byOCAFLabels=True):
    nbFaces = OCC.getNbFaces(hook)
    # compute bbox of faces
    bb = {}
    for i in range(nbFaces):
        bb[i] = OCC.getBoundingBox(hook, [i+1])

    # check face overlap
    intersectings = []; overlaps = []
    if byOCAFLabels: # by OCAF Label components
        found = {}
        for i in range(nbFaces): found[i+1] = False
        ret = OCC.getFaceNameInOCAF(hook)
        compounds = {}
        for l in range(len(ret)//2):
            label = ret[2*l]
            nos = ret[2*l+1]
            for no in nos: found[no] = True
            compounds[label] = nos
        rest = []
        for no in found:
            if not found[no]: rest.append(no)
        compounds['noLabels'] = rest
        for c in compounds:
            print(f"INFO: overlaps in component {c}...", flush=True)
            nos = compounds[c]; lnos = len(nos)
            for n1 in range(lnos):
                i = nos[n1]-1; bb1 = bb[i]
                for n2 in range(n1+1, lnos):
                    j = nos[n2]-1; bb2 = bb[j]
                    if bb1[3] < bb2[0] or bb1[0] > bb2[3]: continue
                    if bb1[4] < bb2[1] or bb1[1] > bb2[4]: continue
                    if bb1[5] < bb2[2] or bb1[2] > bb2[5]: continue
                    ret = OCC.checkFaceOverlap(hook, i+1, j+1, tol=tol)
                    if ret == -1.:
                        print(f"INFO: Face {i+1} and {j+1} in contact.")
                    elif ret < -1:
                        print(f"INFO: Face {i+1} and {j+1} intersects.")
                        intersectings.append( (i+1,j+1) )
                    elif ret > 0:
                        print(f"INFO: Face {i+1} and {j+1} overlaps.")
                        overlaps.append( (i+1,j+1) )
    else: # for all faces in one go
        for i in range(nbFaces):
            bb1 = bb[i+1]
            for j in range(i+1, nbFaces):
                bb2 = bb[j+1]
                if bb1[3] < bb2[0] or bb1[0] > bb2[3]: continue
                if bb1[4] < bb2[1] or bb1[1] > bb2[4]: continue
                if bb1[5] < bb2[2] or bb1[2] > bb2[5]: continue
                ret = OCC.checkFaceOverlap(hook, i+1, j+1, tol=tol)
                if ret == -1.:
                    print(f"INFO: Face {i+1} and {j+1} in contact.")
                elif ret < -1:
                    print(f"INFO: Face {i+1} and {j+1} intersects.")
                    intersectings.append( (i+1,j+1) )
                elif ret > 0:
                    print(f"INFO: Face {i+1} and {j+1} overlaps.")
                    overlaps.append( (i+1,j+1) )
    return overlaps, intersectings

# Check CAD
def checkCAD(hook, tol=1.e-9, byOCAFLabels=True, repair=False):
    import numpy
    score = 0

    #=======================
    # check global dimension
    #=======================
    bb = OCC.getBoundingBox(hook)
    print(f"INFO: BBox X: {bb[0]} -> {bb[3]}")
    print(f"INFO: BBox Y: {bb[1]} -> {bb[4]}")
    print(f"INFO: BBox Z: {bb[2]} -> {bb[5]}")

    bbm = bb[3]-bb[0]
    bbm = bbm + (bb[4]-bb[1])
    bbm = bbm + (bb[5]-bb[2])
    power = int(numpy.log10(bbm/3))
    tol = tol*numpy.pow(10., power)
    print("INFO: adjusting tol: ", tol)

    #============================
    # check for degenerated edges
    #============================
    print("INFO: checking degenerated egdes...", flush=True)
    degenEdges = getDegeneratedEdges(hook, tol=tol)
    if len(degenEdges) > 0: print("Warning: degenerated edges: ", degenEdges)
    else:
        print("INFO: NONE.")
    if repair:
        print("INFO: degenerated edges are bypassed by mesher. No repair needed.")

    #============================
    # check for degenerated faces
    #============================
    print("INFO: checking degenerated faces...", flush=True)
    degenFaces = getDegeneratedFaces(hook, tol=tol)
    if len(degenFaces) > 0: print("Warning: degenerated faces: ", degenFaces)
    else:
        print("INFO: NONE.")
        score += 1
    if repair:
        print("INFO: removing degen faces.")
        OCC._removeFaces(hook, list(degenFaces))
        OCC._sewing(hook, tol=tol)

    #=======================
    # check for lonely edges
    #=======================
    print("INFO: checking lonely/multiple edges...", flush=True)
    ebf = getEdgesByFace(hook)
    fbe = getFacesByEdge(hook, ebf)
    edges1 = []; edgesm = []
    for e in fbe:
        faces = fbe[e]
        if len(faces) == 1: edges1.append(e)
        if len(faces) > 2: edgesm.append(e)
    blocked = set(degenEdges)
    edges1 = [x for x in edges1 if x not in blocked]
    edgesm = [x for x in edgesm if x not in blocked]

    if len(edges1) > 0:
        print("Warning: lonely edges (connected to one face): ", edges1)
    elif len(edgesm) > 0:
        print("Warning: multiple edges (connected to > 2 faces): ", edgesm)
    else:
        print("INFO: NONE.")
        score += 1
    if repair:
        if len(edges1) > 0:
            print("ERROR: you must fill gaps or holes manually.")
        if len(edgesm) > 0:
            print("ERROR: you must delete non manifold faces manually.")

    #=======================
    # check for face overlap
    #=======================
    print("INFO: checking face overlap...", flush=True)
    overlaps, intersectings = getFaceOverlap(hook, tol=tol, byOCAFLabels=byOCAFLabels)
    if len(intersectings) == 0 and len(overlaps) == 0:
        print("INFO: NONE.")
        score += 1
    if repair:
        if len(intersectings) > 0:
            print("ERROR: you must trim faces manually: ", intersectings)
        if len(overlaps) > 0:
            print("ERROR: you must remove overlap faces manually: ", overlaps)

    print(f"score: {score}/3")
    if score < 3: return False
    else: return True

# check CAD through meshing (coarse)
# generate:
# surface mesh: surface.cgns
# surface components; component.cgns
# internal mesh: mesh.cgns
def checkMesh(hook, tol=1.e-9, byOCAFLabels=True, repair=False):
    import Transform.PyTree as T
    import Converter.PyTree as C
    import Generator.PyTree as G
    import Post.PyTree as P
    import Converter.Internal as Internal
    import numpy

    #========
    # meshing
    #========
    print("INFO: meshing...", flush=True)
    (hmin,hmax,hausd) = OCC.occ.analyseEdges(hook)
    #t = OCC.meshAll(hook, hmin=hmax, hmax=hmax, hausd=hausd) # constant hmax
    t = OCC.meshAll(hook, hmin=hmin, hmax=hmax*2., hausd=hausd*0.1) # variable h
    #t = OCC.meshAllOCC(hook, hausd=hausd*0.1, angularDeflection=10.)
    C.convertPyTree2File(t, 'surface.cgns')

    #==============
    # is watertight
    #==============
    print("INFO: check if CAD is watertight...", flush=True)
    a = OCC.getComponents(t, tol=hmin*1.e-3, byOCAFLabels=byOCAFLabels)
    C.convertPyTree2File(a, 'components.cgns')
    print("INFO: find %d component(s)."%len(a))
    watertight = numpy.empty((len(a)), dtype=numpy.int32)
    watertight[:] = False
    for c, i in enumerate(a):
        leaks = []
        ret = OCC.isWatertight(i, leaks, tol=tol)
        if ret:
            print(f"INFO: component {c} is watertight.")
            watertight[c] = True
        else:
            print(f"Warning: component {c} is NOT watertight ({len(leaks)} holes).")
            watertight[c] = False

    if repair:
        # try automatic hole filling
        for c, i in enumerate(a):
            if not watertight[c]:
                ws = P.exteriorFaces(i)
                ws = T.join(ws)
                ws = T.splitConnexity(ws)
                for w in ws:
                    ext = P.exteriorFaces(w)
                    if ext == []:
                        p = G.fittingPlaster(w)
                        g = G.gapfixer(w, p)
                        i = T.join(i, g, tol=tol)
                pc = P.exteriorFaces(i)
                if pc == []:
                    print("INFO: component {c} has been closed.")

    #===============================
    # reorder stability by component
    #===============================
    orderStable = numpy.empty((len(a)), dtype=numpy.int32)
    orderStable[:] = False
    for c, i in enumerate(a):
        T._reorder(i, (1,))
        e1 = Internal.getElementNodes(i)
        j = T.reorder(i, (1,))
        e2 = Internal.getElementNodes(j)
        #result = e1[0][1] == e2[0][1]
        ret = numpy.array_equal(e1[0][1], e2[0][1])
        if not ret:
            print(f"Warning: component {c} is not reorder stable.")
            orderStable[c] = False
        else: orderStable[c] = True

    #================
    # quality by face
    #================
    FACES = Internal.getNodesFromName1(t, 'FACES')
    zones = Internal.getZones(FACES)
    for c, z in enumerate(zones):
        ret = G.checkMesh(z)
        vmin = ret['vmin']
        if vmin < tol:
            print(f"Warning: face {c+1} contains small cells.")
        rmin = ret['rmin']; rmax = ret['rmax']; rcrit = ret['rcrit'] # regularity
        if rmax > rcrit:
            print(f"Warning: regularity of face {c+1} is critical: {rmax}.")
        omin = ret['omin']; omax = ret['omax']; ocrit = ret['ocrit'] # cell skewness
        if omax > ocrit:
            print(f"Warning: cell skewness of face {c+1} is critical: {omax}.")

    #=====================
    # quality by component
    #=====================
    volmin = numpy.empty((len(a)), dtype=numpy.int32)
    volmin[:] = False
    for c, i in enumerate(a):
        ret = G.checkMesh(i)
        vmin = ret['vmin']
        if vmin < tol: volmin[c] = False
        else: volmin[c] = True

    # inside mesher by component
    out = []
    for c, i in enumerate(a):
        if watertight[c] and orderStable[c] and volmin[c]:
            try:
                b = G.tetraMesher(i)
                if C.getNPts(b) > 0: out.append(b)
            except:
                print(f"Error: component {c}: tetgen fails.")
        else:
            print(f"Warning: component {c} not meshed (", end="")
            if not watertight[c]: print("not watertight ", end="")
            if not orderStable[c]: print("not order stable ", end="")
            if not volmin[c]: print("has too small vol ", end="")
            print(")")
    if len(out) > 0:
        C.convertPyTree2File(out, 'mesh.cgns')
