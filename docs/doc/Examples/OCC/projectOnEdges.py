# - projectOnEdges (array) -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

# project on certain CAD edges
a = G.cart((-25,-25,1), (1,1,1), (25,1,1))
OCC.occ.projectOnEdges(hook, a, [1])
C.convertArrays2File(a, 'out.plt')
