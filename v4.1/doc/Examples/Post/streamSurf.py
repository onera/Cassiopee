# - streamSurf (array) -
import Converter as C
import Post as P
import Generator as G
import Geom as D
import math as M

ni = 30; nj = 40

# Node mesh
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,5))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,5))
b = D.line((0.1,5.,0.1), (0.1,5.,3.9), N=10)
b = C.convertArray2Tetra(b)
m = [m1,m2]
def F(x): return M.cos(x)

m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p = P.streamSurf(m, b,['rou','rov','row'])
C.convertArrays2File(m+[p], 'out.plt')
