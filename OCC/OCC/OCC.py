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

# Mailleur de CAD structure
def meshSTRUCT(fileName, format="fmt_step"):
    import Generator, Converter
    hook = occ.readCAD(fileName, format)
    nbFaces = occ.getNbFaces(hook)
    
    out = []
    for i in range(nbFaces):
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, 10)
        print("Face %d has %d edges."%(i+1,len(edges)))
        # edges dans espace uv
        edges = switch2UV(edges)
        # TFI dans espace uv
        try: 
            a = Generator.TFI(edges)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
        except: Converter.convertArrays2File(edges, "edges%d.plt"%i)
    return out
    
# Mailleur CAD non structure
def meshTRI(fileName, format="fmt_step"):
    import Generator, Transform, Converter
    hook = occ.readCAD(fileName, format)
    nbFaces = occ.getNbFaces(hook)

    out = []
    for i in range(nbFaces):
        # edges de la face i
        edges = occ.meshEdgesByFace(hook, i+1, 10)
        # edges dans espace uv
        edges = switch2UV(edges)
        # Delaunay dans espace uv
        edges = Converter.convertArray2Tetra(edges)
        edges = Transform.join(edges)
        edges = Generator.close(edges, 1.e-6) # fix?
        try: 
            a = Generator.T3mesher2D(edges)
            # evaluation sur la CAD
            o = occ.evalFace(hook, a, i+1)
            out.append(o)
        except: Converter.convertArrays2File(edges, 'edges%d.plt'%i)
    return out
