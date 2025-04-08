# - computeDiv (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as T

def Fx(x,y,z): return 2*x+x*y
def Fy(x,y,z): return 4.*y
def Fz(x,y,z): return x*y+z*z

# cas 3D structure
ni = 30; nj = 40; nk = 3
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'veloX', Fx, ['x','y','z'])
m = C.initVars(m, 'veloY', Fy, ['x','y','z'])
m = C.initVars(m, 'veloZ', Fz, ['x','y','z'])
p = P.computeDiv(m, ['veloX','veloY','veloZ']) # p is defined on centers
T.testA([p], 1)

# cas 2D structure
ni = 10; nj = 20; nk = 1
m2 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m2 = C.initVars(m2, 'veloX', Fx, ['x','y','z'])
m2 = C.initVars(m2, 'veloY', Fy, ['x','y','z'])
m2 = C.initVars(m2, 'veloZ', Fz, ['x','y','z'])
p2 = P.computeDiv(m2, ['veloX','veloY','veloZ']) # p is defined on centers
T.testA([p2], 2)

# cas 2D structure (xz)
ni = 10; nj = 1; nk = 20
m21 = G.cart((0,0,0), (10./(ni-1),1,10./(nk-1)), (ni,nj,nk))
m21 = C.initVars(m21, 'veloX', Fx, ['x','y','z'])
m21 = C.initVars(m21, 'veloY', Fy, ['x','y','z'])
m21 = C.initVars(m21, 'veloZ', Fz, ['x','y','z'])
p21 = P.computeDiv(m21, ['veloX','veloY','veloZ']) # p is defined on centers
T.testA([p21], 21)

# cas 2D structure (yz)
ni = 1; nj = 10; nk = 20
m22 = G.cart((0,0,0), (1,10./(nj-1),10./(nk-1)), (ni,nj,nk))
m22 = C.initVars(m22, 'veloX', Fx, ['x','y','z'])
m22 = C.initVars(m22, 'veloY', Fy, ['x','y','z'])
m22 = C.initVars(m22, 'veloZ', Fz, ['x','y','z'])
p22 = P.computeDiv(m22, ['veloX','veloY','veloZ']) # p is defined on centers
T.testA([p22], 22)

# test sur une liste
P = P.computeDiv([m,m2], ['veloX','veloY','veloZ']) # p is defined on centers
T.testA(P, 3)
