# - streamRibbon (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import math as M

def F(x): return M.cos(x)
ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m2 = G.cart((5.5,0,0), (9./(ni-1),9./(nj-1),1), (ni,nj,2))
t = C.newPyTree(['Base','StreamR']); t[2][1][2] = [m1,m2]
t = C.initVars(t, 'vx', 1.)
t = C.initVars(t, 'vy', F, ['CoordinateX'])
t = C.initVars(t, 'vz', 0.)
x0=0.1; y0=5.; z0=0.5
p = P.streamRibbon(t, (x0,y0,z0),(0.,0.2,0.),['vx','vy','vz'])
C.convertPyTree2File(p, "out.cgns")
