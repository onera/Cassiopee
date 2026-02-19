# scanner
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.Vector as Vector
import OCC

def checkOrder(z, axis):
    P1 = C.getValue(z, 'GridCoordinates', 0)
    P2 = C.getValue(z, 'GridCoordinates', 1)
    P3 = C.getValue(z, 'GridCoordinates', 2)
    v1 = Vector.sub(P2, P1)
    v2 = Vector.sub(P3, P1)
    v3 = Vector.cross(v1, v2)
    s = Vector.dot(v3, axis)
    if s > 0: return True
    else: return False

def getSliceZ(t, zval, axis):
    iso = P.isoSurfMC(t, 'CoordinateZ', zval)
    iso = T.join(iso)
    isos = T.splitConnexity(iso)
    isos = C.convertBAR2Struct(isos)
    for c, iso in enumerate(isos):
        if not checkOrder(iso, axis):
            isos[c] = T.reorder(iso, (-1,2,3))
    return isos

#=============================================
# scan a surface mesh, return a CAD surface
# IN: axis: scan axis
# IN: nz: number of discretisation
#=============================================
def scan(t, axis, nz):

    hook = OCC.occ.createEmptyCAD("empty.step", "fmt_step")

    # rotate t to get axis in z
    T._rotate(t, axis, (0,0,1))

    bb = G.bbox(t)
    zmin = bb[2]+1.e-6; zmax = bb[5]-1.e-6
    for n in range(nz):
        zval = zmin + n*(zmax-zmin)/(nz-1.)
        iso = getSliceZ(t, zval, axis)[0]
        isof = C.getAllFields(iso, 'nodes', api=1)[0]
        OCC.occ.addSpline(hook, isof[1], 1)
    OCC.occ.loft(hook, [i for i in range(1,nz+1)], [])
    return hook
