# - integMoment (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

ni = 11; nj = 11
def f1(x,y): return 2*x + y
def f2(x,y): return 3*x*y + 4

# STRUCT 2D
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
C._initVars(m,'vx', f1, ['CoordinateX','CoordinateY'])
C._initVars(m,'vy', f2, ['CoordinateX','CoordinateY'])
C._initVars(m,'vz', 1.)
res = P.integMoment(m, center=(5.,5.,1.), vector=['vx','vy','vz'])
test.testO(res,1)

#TRI
ni = 30; nj = 40
m2 = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
C._initVars(m2,'vx', f1, ['CoordinateX','CoordinateY'])
C._initVars(m2,'vy', f2, ['CoordinateX','CoordinateY'])
C._initVars(m2,'vz', 1.)
res = P.integMoment(m2,center=(5.,5.,1.),vector=['vx','vy','vz'])
test.testO(res,2)

# ARBRE
t = C.newPyTree(['Base',2,'Base2',2])
t[2][1] = C.addState(t[2][1], 'Mach', 0.6); t[2][1][2].append(m)
t[2][2] = C.addState(t[2][2], 'Mach', 0.6); t[2][2][2].append(m2)
res = P.integMoment(t,center=(5.,5.,1.),vector=['vx','vy','vz'])
test.testO(res,3)
