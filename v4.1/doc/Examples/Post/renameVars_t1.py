# - renameVars (array) -
import Converter as C
import Post as P
import Generator as G
import KCore.test as test
ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'ro', 1.)
m = C.initVars(m, 'rou', 1.)
m = C.initVars(m, 'rov', 0.)
m = C.initVars(m, 'roE', 1.)
m2 = P.renameVars(m, ['ro','rou'],['Density','MomentumX'])
test.testA([m2],1)

# sur une liste de zones - variables differentes
m2 = G.cart((10,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
C._addVars(m2, ['ro','rou','rov','roE'])
A = [m,m2]
B = P.renameVars(A, ['ro','rou','rov'],['Density','MomentumX','MomentumY'])
test.testA(B,2)
