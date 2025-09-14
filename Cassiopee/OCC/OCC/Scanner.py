# scanner
import Converter.PyTree as C
import Post.PyTree as P
import Transform.PyTree as T
import Generator.PyTree as G
import OCC

def getSlice(t, zval):
    iso = P.isoSurfMC(t, 'CoordinateZ', zval)
    iso = T.join(iso)
    iso = C.convertBAR2Struct(iso)
    # order must be kept consistent
    return iso

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
        iso = getSlice(t, zval)
        isof = C.getAllFields(iso, 'nodes', api=1)[0]
        OCC.occ.addSpline(hook, isof[1], 1)
    OCC.occ.loft(hook, [i for i in range(1,nz+1)], [])
    return hook
