# - computeDiff (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import Connector.PyTree as X
import KCore.test as test

def F(x):
    if x > 5.: return True
    else: return False
def F2(x):
    return x
def celln(y):
    if y > 5.: return True
    else: return False

#---------------------
# noeuds 2D structure
#---------------------
ni = 30; nj = 40; nk = 1
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m = C.initVars(m, 'Density', F, ['CoordinateX'])
m = P.computeDiff(m, 'Density')
test.testT(m,1)
#
# Prise en compte du cellN
#
m = C.initVars(m, 'cellN', celln, ['CoordinateY'])
m = P.computeDiff(m,'Density')
test.testT(m,2)

#
# Prise en compte des raccords pour les champs aux centres des cellules
#
m1 = G.cart((0.,0.,0.), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m2 = G.cart((10.,0.,0.), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
m1 = C.initVars(m1, 'Density', F2, ['CoordinateX'])
m2 = C.initVars(m2, 'Density', F2, ['CoordinateX'])
m1 = C.node2Center(m1,'Density')
m2 = C.node2Center(m2,'Density')
t = C.newPyTree(['Base',2]);
t[2][1][2] += [m1,m2]
t = C.rmVars(t,'Density')
t = X.connectMatch(t, dim=2)
t = P.computeDiff(t,'centers:Density')
test.testT(t,21)

#------------------------
# 3D structure aux noeuds
#------------------------
ni = 30; nj = 40; nk = 11
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
m = C.initVars(m, 'Density', F, ['CoordinateX'])
m = P.computeDiff(m,'Density')
test.testT(m,3)
#
# Prise en compte du cellN
#
m = C.initVars(m, 'cellN', celln, ['CoordinateY'])
m = P.computeDiff(m,'Density')
test.testT(m,4)

#
# Prise en compte des raccords pour les champs aux centres des cellules
#
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
m2 = G.cart((10.,0,0), (10./(ni-1),10./(nj-1),10./(nk-1)), (ni,nj,nk))
m1 = C.initVars(m1, 'Density', F2, ['CoordinateX'])
m2 = C.initVars(m2, 'Density', F2, ['CoordinateX'])
m1 = C.node2Center(m1,'Density')
m2 = C.node2Center(m2,'Density')
t = C.newPyTree(['Base'])
t[2][1][2] += [m1,m2]
t = C.rmVars(t,'Density')
t = X.connectMatch(t, dim=3)
t = P.computeDiff(t,'centers:Density')
test.testT(t,41)
