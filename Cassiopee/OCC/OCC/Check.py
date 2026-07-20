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
            out[j] = i
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
        area = OCC.getFaceArea(hook)
        if area <= tol: out.append(i+1)
    return out

# main function
def check(hook, tol=1.e-12):
    # check for degenerated edges
    ret = getDegeneratedEdges(hook, tol=tol)
    if len(ret) > 0:
        print("Warning: degenerated edges: ", ret)
    # check for degenerated faces
    ret = getDegeneratedFaces(hook, tol=tol)
    if len(ret) > 0:
        print("Warning: degenerated faces: ", ret)
    # check for lonely edges
    ebf = getEdgesByFace(hook)
    fbe = getFacesByEdge(hook, ebf)
    e1 = []; em = []
    for e in fbe:
        faces = fbe[e]
        if len(faces) == 1: e1.append(e)
        if len(faces) > 2: em.append(e)
    if len(e1) > 0:
        print("Warning: lonely edges (connected to one face): ", e1)
    if len(em) > 0:
        print("Warning: multiple edges (connected to > 2 faces): ", em)
    return None
