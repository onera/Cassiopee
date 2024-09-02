# - projectOnFaces (array) -
import OCC
import Converter as C
import Generator as G
import KCore.test as test

hook = OCC.occ.readCAD("cube.step", "fmt_step")

# project on all CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC.occ.projectOnFaces(hook, a, None)
test.testA(a, 1)

# project on certain CAD faces
a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC.occ.projectOnFaces(hook, a, [1,2,3,4,5,6])
test.testA(a, 2)