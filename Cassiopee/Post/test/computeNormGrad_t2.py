# - computeNormGrad 3D 2D non structure -
import Converter as C
import Post as P
import Generator as G
import KCore.test as T

def F(x,y): return 2*x+x*y

# cas 3D hexa
ni = 30; nj = 40; nk = 3
m = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeNormGrad(m, 'ro') # p is defined on centers
T.testA([p], 1)

# cas 3D tetra
ni = 30; nj = 40; nk = 3
m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeNormGrad(m, 'ro') # p is defined on centers
T.testA([p], 2)
#
# cas 3D penta
ni = 30; nj = 40; nk = 3
m = G.cartPenta((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeNormGrad(m, 'ro') # p is defined on centers
T.testA([p], 3)
#
# cas 2D QUAD
ni = 30; nj = 40; nk = 1
m = G.cartHexa((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeNormGrad(m, 'ro') # p is defined on centers
T.testA([p], 4)

# cas 2D tri
ni = 30; nj = 40; nk = 1
m = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'ro', F, ['x','y'])
p = P.computeNormGrad(m, 'ro') # p is defined on centers
T.testA([p], 5)
