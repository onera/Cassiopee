# - evalFace (array) -
import OCC
import Converter as C
import Generator as G

hook = OCC.occ.readCAD("cube.step", "fmt_step")

N=10; h = 1./(N-1)
a = G.cart((0,0,0), (h,h,1), (N,N,1))
a1 = OCC.occ.evalFace(hook, a, 1)

C.convertArrays2File(a1, 'out.plt')
