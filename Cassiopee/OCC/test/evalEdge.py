# - evalEdge (array) -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

N=10; h = 1./(N-1)
a = G.cart((0,0,0), (h,1,1), (N,1,1))
a1 = OCC.occ.evalEdge(hook, a, 1)

C.convertArrays2File(a1, 'out.plt')