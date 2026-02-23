# - computeDiv (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def Fx(x,y,z): return (x-10)*(x-10)+(y-5)*(y-5)
def Fy(x,y,z): return 4.*(y-3)
def Fz(x,y,z): return (x-10)*(y-2)+z*z

#-------------------------------
# 2D structure + raccords match
#-------------------------------
ni = 10; nj = 10; nk = 1
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))
a = C.initVars(a, 'fldX', Fx, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'fldY', Fy, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'fldZ', Fz, ['CoordinateX','CoordinateY','CoordinateZ'])
b = G.cart((9,0,0), (1,1,1), (ni,nj,nk))
b = C.initVars(b, 'fldX', Fx, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'fldY', Fy, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'fldZ', Fz, ['CoordinateX','CoordinateY','CoordinateZ'])

a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imax', b, 'imin',
                 trirac=[1,2])
b = C.addBC2Zone(b, 'match2', 'BCMatch', 'imin', a, 'imax',
                 trirac=[1,2])

t = C.newPyTree(['Base',2]); t[2][1][2] += [a,b]
t = P.computeDiv(t, 'fld')
test.testT(t,1)
