# - streamLine (array) -
import Converter as C
import Post as P
import Generator as G
import math as M

ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,2))

def F(x): return M.cos(x)

m = [m1,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
x0=0.1; y0=5.; z0=0.5; p = P.streamLine(m, (x0,y0,z0),['rou','rov','row'])
C.convertArrays2File(m+[p], 'out.plt')
