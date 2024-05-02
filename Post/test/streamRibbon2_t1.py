# - streamRibbon2 (array) -
import Converter as C
import Post as P
import Generator as G
import math as M
import KCore.test as test

def F(x): return M.cos(x)
def H(x): return M.sin(x)

ni = 30; nj = 40; nk = 10;
m1 = G.cart((0,0,-1.), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))


x0=0.1; y0=5.; z0=0.5

print(10*'#'+" Calcul sur maillage structure "+10*'#')
m = [m1]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', H, ['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 1)
#C.convertArrays2File(p2, 'outStruct.plt')

print(10*'#'+" Calcul sur maillage non structure hexa "+10*'#')
m2 = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', H, ['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 2)
#C.convertArrays2File(p2, 'outHexa.plt')

print(10*'#'+" Calcul sur maillage non structure tetra "+10*'#')
m3 = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m3,]#,m2]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', H, ['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 3)
#C.convertArrays2File(p2, 'outTetra.plt')

print(10*'#'+" Calcul sur maillage non structure penta "+10*'#')
m4 = G.cartPenta((0,0,0), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m4]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', H, ['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 4)
#C.convertArrays2File(p2, 'outPenta.plt')

print(10*'#'+" Calcul sur maillage non structure Pyra "+10*'#')
m5 = G.cartPyra((0,0,0), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m5]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F,['x'])
m = C.initVars(m, 'row', H,['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 5)
#C.convertArrays2File(p2, 'outPyra.plt')

print(10*'#'+" Calcul sur maillage non structure NGon "+10*'#')
m6 = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),4./(nk-1)), (ni,nj,nk))
m = [m6]
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', F, ['x'])
m = C.initVars(m, 'row', H, ['x'])
p2 = P.streamRibbon2(m, (x0,y0,z0),['rou','rov','row'], dir=1,width=0.5)
test.testA(p2, 6)
#C.convertArrays2File(p2, 'outNGon.plt')
