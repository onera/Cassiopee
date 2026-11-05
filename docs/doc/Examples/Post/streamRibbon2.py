# - streamRibbon2 (array) -
import Converter as C
import Post as P
import Generator as G
import numpy

ni = 30; nj = 40; nk = 10;
def F(x): return numpy.cos(x)
def H(x): return numpy.sin(x)

m1 = G.cart((0,0,-1.), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m2 = G.cart((5.5,0,-1.), (9./(ni-1),9./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m1,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'], isVectorized=True)
m = C.initVars(m, 'row', H, ['x'], isVectorized=True)
x0=0.1; y0=5.; z0=0.; p = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'],dir=1,width=0.5)[0]
C.convertArrays2File(m+[p], 'out.plt')
