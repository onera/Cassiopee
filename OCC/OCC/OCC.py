"""OpenCascade definition module.
"""
__version__ = '3.5'
__author__ = "Sam Landier"

from . import occ

import Converter
import KCore

__all__ = ['convertCAD2Arrays', 'switch2UV', '_scaleUV', '_unscaleUV',
'allTFI', 'meshSTRUCT', 'meshSTRUCT__', 'meshTRI', 'meshTRI__', 
'meshTRIHO', 'meshQUAD', 'meshQUAD__', 'meshQUADHO', 'meshQUADHO__']

# algo=0: mailleur open cascade (chordal_error)
# algo=1: algorithme T3mesher (h, chordal_error, growth_ratio)
# algo=2: algorithme T3mesher (h, chordal_error, growth_ratio, merge_tol)
def convertCAD2Arrays(fileName, format='fmt_iges', 
                      h=0., chordal_err=0., growth_ratio=0., 
                      merge_tol=-1, algo=1):
    """Convert a CAD (IGES or STEP) file to arrays.
    Usage: a = convertCAD2Arrays(fileName, options)"""
    if algo == 0: # pure OCC
        if chordal_err == 0.: chordal_err = 1.
        a = occ.convertCAD2Arrays0(fileName, format, "None", "None", chordal_err)
        try: import Generator; a = Generator.close(a)
        except: pass
        return a
    elif algo == 1: # OCC+T3Mesher
    	return  occ.convertCAD2Arrays1(fileName, format, h, chordal_err, growth_ratio)
    else: # OCC+T3Mesher v2
    	return  occ.convertCAD2Arrays2(fileName, format, h, chordal_err, growth_ratio, merge_tol)

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

# Build a TFI for a set of edges
# IN: edges: list of arrays defining a loop
# OUT: list of surface meshes
def allTFI(edges):
    import Generator
    nedges = len(edges)
    if nedges == 4:
        return [Generator.TFI(edges)]
    elif nedges == 1:
        return Generator.TFIO(edges[0])
    elif nedges == 2:
        return Generator.TFIHalfO(edges[0], edges[1])
    elif nedges == 3:
        return Generator.TFITri(edges[0],edges[1],edges[2])
    else:
        # TFIStar
        #return Generator.TFIStar2(edges)
        return Generator.TFIStar(edges)
        
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
    import Generator, Converter
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N, -1.)
        #print("Face %d has %d edges."%(i+1,len(edges)))
        # edges dans espace uv
        edges = switch2UV(edges)
        # scale uv
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-6) # the weakness
        # TFI dans espace uv
        try:
            als = allTFI(edges)
            # unscale uv
            _unscaleUV(als, T)
            for a in als:
                # evaluation sur la CAD
                o = occ.evalFace(hook, a, i+1)
                out.append(o)
                if faceNo is not None: faceNo.append(i+1)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, "edges%d.plt"%i)
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
def meshTRI__(hook, N=11, hmax=-1., order=1, faceSubset=None, faceNo=None):
    """Return a TRI discretisation of CAD."""
    if hmax > 0: out = meshTRIH__(hook, hmax, faceSubset, faceNo)
    else: out = meshTRIN__(hook, N, order, faceSubset, faceNo)
    return out

# mesh with constant N
def meshTRIN__(hook, N=11, order=1, faceSubset=None, faceNo=None):
    import Generator, Converter, Transform
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # maille les edges de la face i avec N pt et parametres
        edges = occ.meshEdgesByFace(hook, i+1, N, -1.)
        # edges dans espace uv
        edges = switch2UV(edges)
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-4) # the weakness
        # Delaunay dans espace uv
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
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
            Converter.convertArrays2File(edges, 'edges%d.plt'%i)
    return out

# mesh with hmax
def meshTRIH__(hook, hmax=1., faceSubset=None, faceNo=None):
    import Generator, Converter, Transform
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in range(nbFaces):
        # edges de la face i mailles avec hmax et parametres
        edges = occ.meshEdgesByFace(hook, i+1, -1, hmax)
        _scaleUV(edges, vu='u', vv='v')
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-4)
        #edges = Transform.reorder(edges, (-1,))
        #Converter.convertArrays2File(edges, "edges.plt")
        #edges doit contenir les coords + uv normalises pour entrer dans trimesh
        try:
            a = occ.trimesh(hook, edges, i+1, hmax, 0.02)
            out.append(a)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edges%d.plt'%i)
    return out

def meshTRIHO(fileName, format="fmt_step", N=11):
    """Return a TRI HO discretisation of CAD."""
    import Converter
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
    import Converter
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
    import Generator, Converter
    nbFaces = occ.getNbFaces(hook)
    if faceSubset is None: flist = list(range(nbFaces))
    else: flist = faceSubset
    out = []
    for i in flist:
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N, -1.)
        #print("Face %d has %d edges."%(i+1,len(edges)))
        # edges dans espace uv
        edges = switch2UV(edges)
        # scale uv
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-6) # the weakness
        # TFI dans espace uv
        try:
            als = allTFI(edges)
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
            Converter.convertArrays2File(edges, "edges%d.plt"%i)
    return out