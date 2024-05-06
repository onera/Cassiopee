# - projectOnFaces (array) -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

# project on all CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC.occ.projectOnFaces(hook, a, None)
C.convertArrays2File(a, 'out.plt')

# project on certain CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC.occ.projectOnFaces(hook, a, [1,2,3,4,5,6])
C.convertArrays2File(a, 'out.plt')
