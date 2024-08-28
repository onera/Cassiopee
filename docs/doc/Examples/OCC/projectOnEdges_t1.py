# - projectOnEdges (array) -
import OCC
import Generator as G
import KCore.test as test

hook = OCC.occ.readCAD("cube.step", "fmt_step")

# project on certain CAD edges
a = G.cart((-25,-25,1), (1,1,1), (25,1,1))
OCC.occ.projectOnEdges(hook, a, [1])
test.testA(a, 1)