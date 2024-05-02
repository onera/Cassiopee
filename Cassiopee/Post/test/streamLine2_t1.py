# - streamLine2 (array) -
import Converter as C
import Post as P
import Generator as G
import math as M
import KCore.test as test

ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))

def F(x): return M.cos(x)

x0=0.1; y0=5.; z0=0.5

m = [m1]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 1)
#C.convertArrays2File(p2, 'outStruct.plt')

m2 = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = [m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 2)
#C.convertArrays2File(p2, 'outHexa.plt')

m3 = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = [m3,]#,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 3)
#C.convertArrays2File(p2, 'outTetra.plt')

m4 = G.cartPenta((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = [m4]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 4)
#C.convertArrays2File(p2, 'outPenta.plt')

m5 = G.cartPyra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = [m5]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F,['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 5)
#C.convertArrays2File(p2, 'outPyra.plt')

m6 = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = [m6]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', 0.)
p2 = P.streamLine2(m, (x0,y0,z0),['rou','rov','row'])
test.testA(p2, 6)
#C.convertArrays2File(p2, 'outNGon.plt')
