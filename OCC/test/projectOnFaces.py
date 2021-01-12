# - projectOnFaces (array) -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

a = G.cart((-27,-17,-33), (1,1,1), (1,10,10))
OCC.occ.projectOnFaces(hook, a)
C.convertArrays2File(a, 'out.plt')
