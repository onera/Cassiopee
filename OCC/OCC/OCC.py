"""OpenCascade definition module.
"""
__version__ = '3.2'
__author__ = "Sam Landier"

from . import occ

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
# OUT: liste d'arrays STRUCT ayant uv dans x,y
def switch2UV(edges):
    import Converter
    out = []
    for e in edges:
        ni = e[2]; nj = e[3]; nk = e[4]
        uv = Converter.array('x,y,z',ni,nj,nk)
        uv[1][0,:] = e[1][3,:]
        uv[1][1,:] = e[1][4,:]
        uv[1][2,:] = 0.
        out.append(uv)
    return out

# IN: edges: liste d'arrays
def _scaleUV(edges):
    import Converter
    umax = Converter.getMaxValue(edges, 'x')
    umin = Converter.getMinValue(edges, 'x')
    vmax = Converter.getMaxValue(edges, 'y')
    vmin = Converter.getMinValue(edges, 'y')
    du = max(umax-umin, 1.e-10); du = 1./du
    dv = max(vmax-vmin, 1.e-10); dv = 1./dv
    for e in edges:
        e[1][0,:] = (e[1][0,:]-umin)*du
        e[1][1,:] = (e[1][1,:]-vmin)*dv
    return (umin,umax,vmin,vmax)

# IN: edges: liste d'arrays
def _unscaleUV(edges, T):
    (umin,umax,vmin,vmax) = T
    du = umax-umin
    dv = vmax-vmin
    for e in edges:
        e[1][0,:] = e[1][0,:]*du+umin
        e[1][1,:] = e[1][1,:]*dv+vmin
    return None    

# IN: edges: list of arrays defining a loop
# OUT: list of surface meshes
def allTFI(edges):
    import Generator
    nedges = len(edges)
    corners = nedges # not true in general
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
        return Generator.TFIStar(edges)

# Mailleur de CAD structure
def meshSTRUCT(fileName, format="fmt_step", N=11):
    import Generator, Converter
    hook = occ.readCAD(fileName, format)
    nbFaces = occ.getNbFaces(hook)
    
    out = []
    for i in range(nbFaces):
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N)
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
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, "edges%d.plt"%i)
    return out
    
# Mailleur CAD non structure
def meshTRI(fileName, format="fmt_step", N=11):
    import Generator, Transform, Converter
    hook = occ.readCAD(fileName, format)
    nbFaces = occ.getNbFaces(hook)

    out = []
    for i in range(nbFaces):
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, N)
        # edges dans espace uv
        edges = switch2UV(edges)
        T = _scaleUV(edges)
        # force la fermeture de la boucle
        edges = Generator.close(edges, 1.e-4) # the weakness
        # Delaunay dans espace uv
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        try:
            a = Generator.T3mesher2D(edges)
            _unscaleUV([a], T)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
        except Exception as e:
            print(str(e))
            Converter.convertArrays2File(edges, 'edges%d.plt'%i)
    return out
