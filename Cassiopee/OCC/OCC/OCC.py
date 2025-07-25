"""OpenCascade definition module.
"""
__version__ = '4.1'
__author__ = "Sam Landier, Christophe Benoit"

from . import occ

import Converter
import Transform
import Generator
import Geom
import KCore

__all__ = ['convertCAD2Arrays',
           'switch2UV', 'switch2UV2', '_scaleUV', '_unscaleUV',
           'meshSTRUCT', 'meshSTRUCT__', 'meshTRI', 'meshTRI__', 'meshTRIU__',
           'meshTRIHO', 'meshQUAD', 'meshQUAD__', 'meshQUADHO', 'meshQUADHO__',
           'ultimate', 'meshAllEdges', 'meshAllFacesTri', 'meshAllFacesStruct',
           'meshAllFacesTri', 'meshFaceWithMetric', 'identifyTags__',
           '_projectOnFaces',
           'readCAD', 'writeCAD',
           'getNbEdges', 'getNbFaces', 'getFileAndFormat', 'getFaceArea',
           '_translate', '_rotate', '_scale', '_sewing',
           '_splitFaces', '_mergeFaces', '_trimFaces', '_removeFaces',
           'printOCAF', 'getFaceNameInOCAF', 'getEdgeNameInOCAF']

# algo=0: mailleur open cascade (chordal_error)
# algo=1: algorithme T3mesher (h, chordal_error, growth_ratio)
# algo=2: algorithme T3mesher (h, chordal_error, growth_ratio, merge_tol)
def convertCAD2Arrays(fileName, format=None,
                      h=0., chordal_err=0., growth_ratio=0.,
                      merge_tol=-1, algo=1, join=True):
    """Convert a CAD (IGES or STEP) file to arrays.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if format is None: format = Converter.convertExt2Format__(fileName)
    if algo == 0: # pure OCC
        if chordal_err == 0.: chordal_err = 1.
        a = occ.convertCAD2Arrays0(fileName, format, "None", "None", chordal_err)
        a = Generator.close(a)
    elif algo == 1: # OCC+T3Mesher
        a = occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio, join)
    else: # OCC+T3Mesher v2
        a = occ.convertCAD2Arrays2(fileName, format, h, chordal_err, growth_ratio, merge_tol, join)

    # if nothing is read, try to read as edges (suppose draft)
    if Converter.getNPts(a) == 0:
        if h == 0.: h = 1.e-2
        hook = occ.readCAD(fileName, format)
        a = occ.meshGlobalEdges(hook, h)
    return a

# IN: edges: liste d'arrays STRUCT possedant x,y,z,u,v
# OUT: liste d'arrays STRUCT ayant uv dans x,y et z=0
def switch2UV(edges):
    """Switch uv to coordinates."""
    out = []
    for e in edges:
        ni = e[2]; nj = e[3]; nk = e[4]
        uv = Converter.array('x,y,z',ni,nj,nk)
        uv[1][0,:] = e[1][3,:]
        uv[1][1,:] = e[1][4,:]
        uv[1][2,:] = 0.
        out.append(uv)
    return out

def switch2UV2(edges):
    """Switch uv to coordinates keeping uv field."""
    out = []
    for e in edges:
        ni = e[2]; nj = e[3]; nk = e[4]
        uv = Converter.array('x,y,z,u,v',ni,nj,nk)
        uv[1][0,:] = e[1][3,:]
        uv[1][1,:] = e[1][4,:]
        uv[1][2,:] = 0.
        uv[1][3,:] = e[1][3,:]
        uv[1][4,:] = e[1][4,:]
        out.append(uv)
    return out

# Scale u,v in [0,1]
# IN: edges: liste d'arrays
# IN: vu: name for u variable
# IN: vv: name for v variable
# scale entre 0 et 1 les variables vu et vv
def _scaleUV(edges, vu='x', vv='y'):
    """Scale vu and vv in [0,1]."""
    umax = Converter.getMaxValue(edges, vu)
    umin = Converter.getMinValue(edges, vu)
    vmax = Converter.getMaxValue(edges, vv)
    vmin = Converter.getMinValue(edges, vv)
    du = max(umax-umin, 1.e-10); du = 1./du
    dv = max(vmax-vmin, 1.e-10); dv = 1./dv
    for e in edges:
        pu = KCore.isNamePresent(e, vu)
        pv = KCore.isNamePresent(e, vv)
        e[1][pu,:] = (e[1][pu,:]-umin)*du
        e[1][pv,:] = (e[1][pv,:]-vmin)*dv
    return (umin,umax,vmin,vmax)

# unscale u,v back
# IN: edges: liste d'arrays
# IN: T: min max of parameters as returned by scaleUV
def _unscaleUV(edges, T, vu='x', vv='y'):
    """Unscale vu and vv with given minmax."""
    (umin,umax,vmin,vmax) = T
    du = umax-umin
    dv = vmax-vmin
    for e in edges:
        pu = KCore.isNamePresent(e, vu)
        pv = KCore.isNamePresent(e, vv)
        e[1][pu,:] = e[1][pu,:]*du+umin
        e[1][pv,:] = e[1][pv,:]*dv+vmin
    return None

# Mailleur structure de CAD
# IN: N: the number of points for each patch boundary
def meshSTRUCT(fileName, format='fmt_iges', N=11):
    """Return a STRUCT discretisation of CAD."""
    hook = occ.readCAD(fileName, format)
    return meshSTRUCT__(hook, N)

# subfunction
# IN: hook: CAD hook
# IN: faceSubSet: a list of faces to mesh
# OUT: faceNo: keep the CAD face number for each zone
# OUT: one mesh per CAD face
def meshSTRUCT__(hook, N=11, faceSubset=None, faceNo=None):
    """Return a STRUCT discretisation of CAD."""
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N, -1., -1.)
        #print("Face %d has %d edges."%(i+1,len(edges)))
        # edges dans espace uv
        edges = switch2UV(edges)
        # scale uv
        T = _scaleUV(edges)
        # force la fermeture des boucles
        edges = Generator.close(edges, 1.e-6) # the weakness
        # TFI dans espace uv
        try:
            als = Generator.allTFI(edges)
            # unscale uv
            _unscaleUV(als, T)
            for a in als:
                # evaluation sur la CAD
                o = occ.evalFace(hook, a, i+1)
                out.append(o)
                if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, "edgesOfFace%03d.plt"%(i+1))
    return out

# Mailleur CAD non structure TRI
# IN: N: number of points for each face boundary, discarded if < 0
# IN: hmax: mesh step, discarded if < 0
# IN: order: ordre du maillage de sortie
def meshTRI(fileName, format="fmt_step", N=11, hmax=-1., order=1):
    hook = occ.readCAD(fileName, format)
    return meshTRI__(hook, N, hmax)

# subfunction
# IN: hook: CAD hook
# IN: faceSubSet: a list of faces to mesh
# OUT: faceNo: keep the CAD face number for each zone
# OUT: one mesh per CAD face
def meshTRI__(hook, N=11, hmax=-1., hausd=-1., order=1, faceSubset=None, faceNo=None):
    """Return a TRI discretisation of CAD."""
    if hmax > 0 and hausd > 0: out = meshTRIH2__(hook, hmax, hausd, order, faceSubset, faceNo)
    elif hmax > 0 or hausd > 0: out = meshTRIH__(hook, hmax, hausd, order, faceSubset, faceNo)
    else: out = meshTRIN__(hook, N, order, faceSubset, faceNo)
    return out

# reordonne les edges par face pour que le mailleur TRI puisse mailler l'entre deux
# les edges interieurs sont numerotes dans le sens inverse de l'edge exterieur
# limitation : un seul niveau d'edge dans l'edge exterieur
def reorderEdgesByFace__(edges):
    import Post
    import KCore.Vector as Vector
    from operator import itemgetter, attrgetter
    splitEdges = Transform.splitConnexity(edges)
    if len(splitEdges) == 1:
        print("== Single closed curve")
        return edges
    print("== Multiple closed curves")
    # classe les edges par surface englobee
    sortedEdges = []
    for c, e in enumerate(splitEdges):
        #e = Converter.convertBAR2Struct(e)
        #e = Transform.reorder(e, (1,2,3))
        a = Generator.close(e, 1.e-4)
        try:
            a = Generator.T3mesher2D(e, triangulateOnly=1) # must not fail!
            v = Generator.getVolumeMap(a)
            surf = Post.integ([a], [v], [])[0]
            #n = Generator.getNormalMap(a)
            #nz = Converter.getMaxValue(n, 'sz')
            #connect = a[2]
            ## les 3 points du premier triangle
            #ind0 = connect[0,0]
            #ind1 = connect[1,0]
            #ind2 = connect[2,0]
            #print('index', ind0, ind1, ind2)
            #P0 = Converter.getValue(a, ind0-1)
            #P1 = Converter.getValue(a, ind1-1)
            #P2 = Converter.getValue(a, ind2-1)
            #tri = Geom.polyline([P0,P1,P2])
            #print('points', P0,P1,P2)
            #hook = Converter.createHook(es, function='nodes')
            #inds = Converter.identifyNodes(hook, tri)
            #print('in e', inds)
            #order = 1
            #if ind0 == inds[0]:
            #    if ind1 == inds[0]-1: order = -1
            #    elif ind1 == inds[0]+1: order = 1
            #elif ind0 == inds[1]:
            #    if ind1 == inds[1]-1: order = -1
            #    elif ind1 == inds[0]+1: order = 1
            #else: # ind0 = inds[2]
            #    if ind1 == inds[2]-1: order = -1
            #    elif ind1 == inds[2]+1: order = 1
            #print('order', order)
            nz = 1; order = 1
        except:
            print("Warning: finding the exterior edge failed.")
            surf = 0.; nz = 1.; order = 1
        sortedEdges.append((e, surf, nz, order))
    sorted(sortedEdges, key=itemgetter(1), reverse=True)
    for i in sortedEdges: print('surf', i[1], i[2])
    # reorder
    edgesOut = []
    for c, se in enumerate(sortedEdges):
        #b = se[0]
        b = Converter.convertBAR2Struct(se[0])
        PG = Generator.barycenter(b) # must be in curve
        P0 = Converter.getValue(b, 0)
        P1 = Converter.getValue(b, 1)
        P01 = Vector.sub(P1, P0)
        PG0 = Vector.sub(P0, PG)
        cross = Vector.cross(PG0, P01)[2]
        #print('cross', c, cross)
        #cross = se[3] # order
        if cross < 0 and c == 0:  # must be exterior curve
            b = Transform.reorder(b, (-1,2,3))
        elif cross > 0 and c != 0: # must be interior curves
            b = Transform.reorder(b, (-1,2,3))
        edgesOut.append(b)
    edgesOut = Converter.convertArray2Tetra(edgesOut)
    edgesOut = Transform.join(edgesOut)
    edgesOut = Generator.close(edgesOut, 1.e-6)
    return edgesOut

# mesh with constant N
def meshTRIN__(hook, N=11, order=1, faceSubset=None, faceNo=None):
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        print("Meshing Face %d ========================================"%(i+1))
        # maille les edges de la face i avec N pt et parametres
        edges = occ.meshEdgesByFace(hook, i+1, N, -1., -1.)
        # edges dans espace uv
        edges = switch2UV(edges)
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-4) # the weakness
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-6)
        edges = reorderEdgesByFace__(edges)
        try:
            a = Generator.T3mesher2D(edges, grading=1.)
            _unscaleUV([a], T)
            if order > 1: a = Converter.convertLO2HO(a, order=order)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
            if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edgesOfFace%03d.plt'%(i+1))
    return out

# prend le Ue des edges dans globalEdges
def meshTRIU__(hook, globalEdges, order=1, faceSubset=None, faceNo=None):
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # maille les edges de la face i avec le U de l'edge
        edges = occ.meshEdgesByFace2(hook, i+1, globalEdges)
        # edges dans espace uv
        edges = switch2UV(edges)
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-4) # the weakness
        # Delaunay dans espace uv
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-6)
        try:
            a = Generator.T3mesher2D(edges, grading=1.)
            _unscaleUV([a], T)
            if order > 1: a = Converter.convertLO2HO(a, order=order)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
            if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edgesOfFace%03d.plt'%(i+1))
    return out

# mesh with hmax or hausd
def meshTRIH__(hook, hmax=-1., hausd=-1, order=1, faceSubset=None, faceNo=None):
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # edges de la face i mailles avec hmax et parametres
        edges = occ.meshEdgesByFace(hook, i+1, -1, hmax, hausd)
        # edges dans espace uv
        edges = switch2UV(edges)
        T = _scaleUV(edges)
        edges = Generator.close(edges, 1.e-4) # the weakness
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-6)
        edges = reorderEdgesByFace__(edges)
        try:
            a = Generator.T3mesher2D(edges, grading=1.)
            _unscaleUV([a], T)
            if order > 1: a = Converter.convertLO2HO(a, order=order)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
            if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edgesOfFace%03d.plt'%(i+1))

    return out

# using trimesh
def meshTRIH2__(hook, hmax=-1., hausd=-1, order=1, faceSubset=None, faceNo=None):
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        edges = occ.meshEdgesByFace(hook, i+1, -1, hmax, hausd)
        _scaleUV(edges, vu='u', vv='v')
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-6)
        #edges = reorderEdgesByFace__(edges)
        try:
            #edges doit contenir les coords + uv normalises pour entrer dans trimesh
            a = occ.trimesh(hook, edges, i+1, hmax, hausd)
            out.append(a)
            if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edgesOfFace%03d.plt'%(i+1))
    return out

def meshTRIHO(fileName, format="fmt_step", N=11):
    """Return a TRI HO discretisation of CAD."""
    a = convertCAD2Arrays(fileName, format,
                          h=0., chordal_err=0., growth_ratio=0.,
                          merge_tol=-1, algo=2)
    #a = meshTRI(fileName, format, N)
    hook = occ.readCAD(fileName, format)
    out = []
    for c, i in enumerate(a):
        print('Projection %d/%d'%(c,len(a)))
        b = Converter.convertLO2HO(i, order=2)
        occ.projectOnFaces(hook, b, None)
        out.append(b)
    return out

def meshQUADHO(fileName, format="fmt_step", N=11):
    """Return a QUAD HO discretisation of CAD."""
    hook = occ.readCAD(fileName, format)
    return meshQUADHO__(hook, N)

def meshQUADHO__(hook, N=11, faceSubset=None, faceNo=None):
    """Return a QUAD HO discretisation of CAD."""
    if faceNo is None: faceNo = []
    a = meshSTRUCT__(hook, N, faceSubset, faceNo)
    a = Converter.convertArray2Hexa(a)
    out = []
    for c, i in enumerate(a):
        print('Projection %d/%d'%(c,len(a)))
        b = Converter.convertLO2HO(i, order=2)
        occ.projectOnFaces(hook, b, [faceNo[c]])
        out.append(b)
    return out

def meshQUAD(fileName, format="fmt_step", N=11, order=1):
    """Return a QUAD discretisation of CAD."""
    hook = occ.readCAD(fileName, format)
    return meshQUAD__(hook, N, order)

def meshQUAD__(hook, N=11, order=1, faceSubset=None, faceNo=None):
    """Return a QUAD HO discretisation of CAD."""
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N, -1., -1.)
        #print("Face %d has %d edges."%(i+1,len(edges)))
        # edges dans espace uv
        edges = switch2UV(edges)
        # scale uv
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-6) # the weakness
        # TFI dans espace uv
        try:
            als = Generator.allTFI(edges)
            # unscale uv
            _unscaleUV(als, T)
            als = Converter.convertArray2Hexa(als)
            if order > 1: als = Converter.convertLO2HO(als, order=order)
            for a in als:
                # evaluation sur la CAD
                o = occ.evalFace(hook, a, i+1)
                out.append(o)
                if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, "edgesOfFace%03d.plt"%(i+1))
    return out

#===============================================================================
# the ultimate TRI regular mesher
# hmax on all edges, hmax in physical space
#===============================================================================

# enforce edges in face
def _enforceEdgesInFace(a, edges):
    xo = a[1] # array1
    c = 0
    for e in edges:
        xe = e[1]
        npts = xe.shape[1]-1
        xo[0, c:c+npts] = xe[0, 0:npts]
        xo[1, c:c+npts] = xe[1, 0:npts]
        xo[2, c:c+npts] = xe[2, 0:npts]
        c += npts
    return None

# hmax: hmax sur les edges et dans les faces (constant)
# hausd: erreur de corde, pas pris en compte
def ultimate(hook, hmax, hausd=-1, metric=True):
    mesh = []
    FAILED1 = []; FAILED2 = []
    nbFaces = occ.getNbFaces(hook)

    for i in range(nbFaces):
        print('Meshing face %d ======================'%i)

        # maillage des edges de la Face (sortie par wire)
        wires = occ.meshEdgesByFace3(hook, i+1, hmax, hausd)

        # sortie a plat de tous les edges
        #plat = []
        #for w in wires: plat += w
        #Converter.convertArrays2File(plat, '%03d_edgeXY.plt'%i)

        # join des edges par wire (structured)
        edges = []
        for w in wires:
            e = Transform.join(w)
            edges.append(e)

        # sauvegarde des edges
        edgesSav = []
        for e in edges: edgesSav.append(Converter.copy(e))

        # TRIMESH METRIC TRY
        SUCCESS = False
        if metric:
            SUCCESS = meshFaceWithMetric(hook, i, edges, hmax, hausd, 1.1, mesh, FAILED1)

        if not SUCCESS: # TRIMESH sans metric
            edges = edgesSav
            meshFaceInUV(hook, i, edges, 1., mesh, FAILED2)

    FAIL1 = len(FAILED1)
    print("METRICFAILURE = %d / %d"%(FAIL1, nbFaces))
    for f in FAILED1:
        print("METRICFAILED on face = %03d_edgeUV.plt"%f)
    FAIL2 = len(FAILED2)
    print("FINAL FAILURE = %d / %d"%(FAIL2,nbFaces))
    for f in FAILED2:
        print("FINAL FAILED on face = %03d_edgeUV.plt"%f)

    return mesh

#=======================================================================================
# BACK AGAIN
#=======================================================================================

#===============================================================================
# TRI Mesh Face no i of CAD from parametrized edges
# IN: hook: cad hook
# IN: i: no de la face
# IN: edges structured one per wire
# hmax: hmin/hmax/hausd par face
#===============================================================================
def meshFaceWithMetric(hook, i, edges, hmin, hmax, hausd, mesh, FAILED):

    # save edges
    edgesSav = []
    for e in edges: edgesSav.append(Converter.copy(e))

    # must close in uv space
    edges = switch2UV2(edges)
    T = _scaleUV(edges)
    edges = Converter.convertArray2Tetra(edges)
    edges = Transform.join(edges)
    edges = Generator.close(edges, 1.e-10)
    _unscaleUV([edges], T)

    pt = edges[1]
    edges = occ.evalFace(hook, edges, i)
    edges = Converter.addVars(edges, ['u','v'])
    edges[1][3,:] = pt[3,:]
    edges[1][4,:] = pt[4,:]

    if edges[2].shape[1] == 0: return True # pass

    # supprime les edges collapsed
    #edges2 = Generator.close(edges, 1.e-6)

    # Scale UV des edges
    _scaleUV([edges], vu='u', vv='v')
    try:
        a = occ.trimesh(hook, edges, i, hmin, hmax, hausd, 1.1)
        _enforceEdgesInFace(a, edgesSav)
        a = Generator.close(a, 1.e-10) # needed for periodic faces
        if occ.getFaceOrientation(hook, i) == 0:
            a = Transform.reorder(a, (-1,))
        mesh.append(a)
        SUCCESS = True
    except Exception as e:
        SUCCESS = False
        Converter.convertArrays2File(edges, '%03d_edgeUV.plt'%i) # pas vraiment UV
        FAILED.append(i)

    return SUCCESS

# TRI mesh face regular in UV space
def meshFaceInUV(hook, i, edges, grading, mesh, FAILED):

    # save edges
    edgesSav = []
    for e in edges: edgesSav.append(Converter.copy(e))

    # Passage des edges dans espace uv
    edges = switch2UV(edges)
    T = _scaleUV(edges)
    edges = Converter.convertArray2Tetra(edges)
    edges = Transform.join(edges)

    # Maillage de la face
    try:
        a = Generator.T3mesher2D(edges, grading=grading)
        _unscaleUV([a], T)
        o = occ.evalFace(hook, a, i)
        _enforceEdgesInFace(o, edgesSav)
        if occ.getFaceOrientation(hook, i) == 0:
            o = Transform.reorder(o, (-1,))
        mesh.append(o)
        SUCCESS = True
    except Exception as e:
        SUCCESS = False
        Converter.convertArrays2File(edges, '%03d_edgeUV.plt'%i)
        FAILED.append(i)

    return SUCCESS

# pointed hat mesh face in UV space
def meshFaceWithPointedHat(hook, i, edges, mesh):

    # save edges
    #edgesSav = []
    #for e in edges: edgesSav.append(Converter.copy(e))

    # Passage des edges dans espace uv
    edges = switch2UV(edges)
    T = _scaleUV(edges)

    # Maillage l'edge le plus long par pointedhat (no fail)
    lmax = -1.
    for e in edges:
        l = Geom.getLength(e)
        if l > lmax:
            lmax = l
            X = Generator.barycenter(e)
            a = Generator.pointedHat(e, X)

    a = Converter.convertArray2Tetra(a)
    a = Generator.close(a, 1.e-10)
    _unscaleUV([a], T)
    o = occ.evalFace(hook, a, i)
    #_enforceEdgesInFace(o, edgesSav)
    if occ.getFaceOrientation(hook, i) == 0:
        o = Transform.reorder(o, (-1,))
    mesh.append(o)

    return True

# mesh all CAD edges with hmin, hmax, hausd
def meshAllEdges(hook, hmin, hmax, hausd, N, edgeList=None):
    """Mesh all CAD edges with hmin, hmax, hausd."""
    if edgeList is None:
        nbEdges = occ.getNbEdges(hook)
        edgeList = range(1, nbEdges+1)
    dedges = []
    for i in edgeList:
        e = occ.meshOneEdge(hook, i, hmin, hmax, hausd, N, None)
        dedges.append(e)
    dedges = Generator.zip(dedges, tol=hmax/100.) # safe and necessary for corner/seam points
    return dedges

#=================================================================
# mesh TRI given CAD faces from discrete edges U + hList
# IN: hook: hook of cad
# IN: hmax: hmax size
# IN: hausd: hausd size
# IN: dedges: list of all discretized edges (all CAD edges)
# IN: metric: if True use metric else mesh in uv
# IN: faceList: list of faces to mesh (start 1)
# IN: hList: list of (hmin, hmax, hausd) for each face to mesh
#==================================================================
def meshAllFacesTri(hook, dedges, metric=True, faceList=[], hList=[]):
    nbFaces = len(faceList)
    FAILED1 = []; FAILED2 = []; dfaces = []
    for c, i in enumerate(faceList):

        print("========== face %d / %d ==========="%(i,nbFaces))

        wires = occ.meshEdgesOfFace(hook, i, dedges)

        # join des edges par wire (structured)
        edges = []
        for w in wires:
            e = Transform.join(w)
            if e != []: edges.append(e)

        # sauvegarde des edges
        edgesSav = []
        for e in edges: edgesSav.append(Converter.copy(e))

        # TRIMESH METRIC TRY
        SUCCESS = False
        if metric:
            hsize = hList[c]
            SUCCESS = meshFaceWithMetric(hook, i, edges, hsize[0], hsize[1], hsize[2], dfaces, FAILED1)

        if not SUCCESS: # TRIMESH sans metric
            edges = edgesSav
            SUCCESS = meshFaceInUV(hook, i, edges, 1., dfaces, FAILED2)

        if not SUCCESS: # pointed hat
            edges = edgesSav
            #dfaces.append(None)
            SUCCESS = meshFaceWithPointedHat(hook, i, edges, dfaces)

    FAIL1 = len(FAILED1)
    print("METRICFAILURE = %d / %d"%(FAIL1, nbFaces))
    for f in FAILED1:
        print("METRICFAILED on face = %03d_edgeUV.plt"%f)
    FAIL2 = len(FAILED2)
    print("FINAL FAILURE = %d / %d"%(FAIL2, nbFaces))
    for f in FAILED2:
        print("FINAL FAILED on face = %03d_edgeUV.plt"%f)

    return dfaces

#==================================================================
# mesh STRUCT given CAD faces from discrete edges U
#==================================================================
def meshAllFacesStruct(hook, dedges, faceList=[]):
    nbFaces = len(faceList)
    FAILED1 = []; dfaces = []
    nloct = []; nofacet = [] # nbre de grilles pour la face c; no de la face
    for c, i in enumerate(faceList):

        print("========== face %d / %d ==========="%(i,nbFaces))

        wires = occ.meshEdgesOfFace(hook, i, dedges)
        edges = wires[0] # only first for now

        SUCCESS = False
        edges = switch2UV(edges)
        # scale uv
        T = _scaleUV(edges)
        # force la fermeture des boucles
        #edges = Generator.close(edges, 1.e-6) # the weakness
        # TFI dans espace uv
        try:
            als = Generator.allTFI(edges)
            # unscale uv
            _unscaleUV(als, T)
            for c, a in enumerate(als):
                # evaluation sur la CAD
                o = occ.evalFace(hook, a, i)
                dfaces.append(o)
                nofacet.append(i)
                nloct.append(c+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, "edgesOfFace%03d.plt"%i)
            FAILED1.append(i)
            dfaces.append(None)
            nofacet.append(i)
            nloct.append(0)

    FAIL1 = len(FAILED1)
    print("STRUCTFAILURE = %d / %d"%(FAIL1, nbFaces))
    for f in FAILED1:
        print("STRUCTFAILED on face = %03d_edgeUV.plt"%f)

    return dfaces, nloct, nofacet

# project arrays on faces
def _projectOnFaces(hook, a, faceList=None):
    """Project arrays on CAD."""
    if isinstance(a[0], list):
        for i in a: occ.projectOnFaces(hook, i, faceList)
    else:
        occ.projectOnFaces(hook, a, faceList)
    return None

#=============================================================================
# CAD fixing
#=============================================================================

# read CAD and return CAD hook
def readCAD(fileName, format='fmt_step'):
    """Read CAD file and return a CAD hook."""
    h = occ.readCAD(fileName, format)
    return h

# write CAD to file
def writeCAD(hook, fileName, format='fmt_step'):
    """Write CAD file."""
    occ.writeCAD(hook, fileName, format)
    return None

# Return the number of edges in CAD hook
def getNbEdges(hook):
    """Return the number of edges in CAD hook."""
    return occ.getNbEdges(hook)

# Return the number of faces in CAD hook
def getNbFaces(hook):
    """Return the number of faces in CAD hook."""
    return occ.getNbFaces(hook)

# Return the file and format used to load CAD in hook
def getFileAndFormat(hook):
    """Return file and format of associated CAD file."""
    return occ.getFileAndFormat(hook)

# Return the area of specified faces
def getFaceArea(hook, faceList=None):
    """Return the area of given faces."""
    return occ.getFaceArea(hook, faceList)

# Translate
def _translate(hook, vector, faceList=None):
    """Translate all or given faces."""
    occ.translate(hook, vector, faceList)
    return None

# Rotate
def _rotate(hook, Xc, axis, angle, faceList=None):
    """Rotate all or given faces."""
    occ.rotate(hook, Xc, axis, angle, faceList)
    return None

# Scale
def _scale(hook, factor, X, faceList=None):
    """Scale all or given faces."""
    occ.scale(hook, factor, X, faceList)
    return None

# sew a set of faces
# faces: face list numbers
def _sewing(hook, faceList=None, tol=1.e-6):
    """Sew some faces (suppress redundant edges)."""
    occ.sewing(hook, faceList, tol)
    return None

# add fillet from edges with given radius
def _addFillet(hook, edges, radius, new2OldEdgeMap=[], new2OldFaceMap=[]):
    occ.addFillet(hook, edges, radius, new2OldEdgeMap, new2OldFaceMap)
    return None

# edgeMap and faceMap are new2old maps
def _removeFaces(hook, faceList, new2OldEdgeMap=[], new2OldFaceMap=[]):
    """Remove given faces."""
    occ.removeFaces(hook, faceList, new2OldEdgeMap, new2OldFaceMap)
    return None

# fill hole from edges
# edges: edge list numbers (must be ordered)
def _fillHole(hook, edges, faceList=None, continuity=0):
    """Fill hole defined by close loop of edges."""
    occ.fillHole(hook, edges, faceList, continuity)
    return None

# trim two set of surfaces
def _trimFaces(hook, faceList1, faceList2):
    """Trim a set of faces with another set of faces."""
    occ.trimFaces(hook, faceList1, faceList2)
    return None

# split all faces
def _splitFaces(hook, area):
    """Split all faces to be less than area."""
    occ.splitFaces(hook, area)
    return None

# merge faces
def _mergeFaces(hook, faceList=None):
    """Merge some faces."""
    occ.mergeFaces(hook, faceList)
    return None

# identify tag component
def identifyTags__(a):
    return occ.identifyTags(a)

# print OCAF document
def printOCAF(h):
    """Print OCAF document."""
    occ.printOCAF(h)

def getFaceNameInOCAF(h):
    """Return face names in OCAF."""
    return occ.getFaceNameInOCAF2(h)

def getEdgeNameInOCAF(h):
    """Return edge names in OCAF."""
    return occ.getEdgeNameInOCAF2(h)
