# - computeDiff (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

def F(x):
    if ( x > 5. ): return True
    else: return False
def celln(y):
    if ( y > 5. ): return True
    else: return False

#-------------------------------
# centres 2D structure
#-------------------------------
ni = 30; nj = 40; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'Density', F, ['CoordinateX'])
m = C.node2Center(m, 'Density')
m = P.computeDiff(m, 'centers:Density')
test.testT(m,1)

# Prise en compte du cellN
m = C.initVars(m, 'cellN', celln, ['CoordinateY'])
m = C.node2Center(m, 'cellN')
m = P.computeDiff(m, 'centers:Density')
test.testT(m,2)

#-------------------------------
# centres 3D structure
#-------------------------------
ni = 30; nj = 40; nk = 11
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1./(nk-1)), (ni,nj,nk))
m = C.initVars(m, 'Density', F, ['CoordinateX'])
m = C.node2Center(m,'Density')
m = P.computeDiff(m,'centers:Density')
test.testT(m,3)
#
# Prise en compte du cellN
m = C.initVars(m, 'cellN', celln, ['CoordinateY'])
m = C.node2Center(m, 'cellN')
m = P.computeDiff(m, 'centers:Density')
test.testT(m,4)
