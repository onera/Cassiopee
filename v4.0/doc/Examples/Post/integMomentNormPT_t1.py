# - integMomentNorm (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

ni = 11; nj = 11
def f1(x,y): return 2*x + y
def f2(x,y): return 3*x*y + 4

# STRUCT 2D
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m = C.initVars(m,'vx', f1, ['CoordinateX','CoordinateY'])
m = C.node2Center(m, Internal.__GridCoordinates__)
m = C.initVars(m,'centers:vy', f2, ['CoordinateX','CoordinateY'])
res = P.integMomentNorm(m,(5.,5.,1.),'vx')+\
    P.integMomentNorm(m,(5.,5.,1.),'centers:vy')
test.testO(res,1)

#TRI
ni = 30; nj = 40
m2 = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
m2 = C.node2Center(m, Internal.__GridCoordinates__)
m2 = C.initVars(m2,'centers:vy', f2, ['CoordinateX','CoordinateY'])
res = P.integMomentNorm(m2,(5.,5.,1.),'vx')+\
    P.integMomentNorm(m2,(5.,5.,1.),'centers:vy')
test.testO(res,2)

# ARBRE
t = C.newPyTree(['Base',2,'Base2',2])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6); t[2][1][2].append(m)
t[2][2] = C.addState(t[2][2], 'Mach', 0.6); t[2][2][2].append(m2)
res = P.integMomentNorm(t,(5.,5.,1.),'vx')+\
    P.integMomentNorm(t,(5.,5.,1.),'centers:vy')
test.testO(res,3)
