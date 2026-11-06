import Converter as C
import Connector as X
import Generator as G
import Geom as D
import KCore.test as test

# 2D QUAD
s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5.)
res = C.initVars(res, 'cellN', 1.)
ca = C.extractVars(res,['cellN'])
res = C.extractVars(res,['x','y','z'])
celln = X.blankCells([res], [ca], [s], blankingType=0, delta=0.)
res2 = C.addVars([res,celln[0]])
test.testA([res2],1)

# 3D HEXA
s = D.sphere((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5.)
res = C.initVars(res, 'cellN', 1.)
ca = C.extractVars(res,['cellN'])
res = C.extractVars(res,['x','y','z'])
celln = X.blankCells([res], [ca], [s], blankingType=0, delta=0.)
res = C.addVars([res,celln[0]])
test.testA([res],2)

# 2D TRI
s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5.)
res = C.convertArray2Tetra(res)
res = C.initVars(res, 'cellN', 1.)
ca = C.extractVars(res,['cellN'])
res = C.extractVars(res,['x','y','z'])
celln = X.blankCells([res], [ca], [s], blankingType=0, delta=0.)
res2 = C.addVars([res,celln[0]])
test.testA([res2],3)

# 3D TETRA
s = D.sphere((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s], [snear], dfar=5.)
res = C.convertArray2Tetra(res)
res = C.initVars(res, 'cellN', 1.)
ca = C.extractVars(res,['cellN'])
res = C.extractVars(res,['x','y','z'])
celln = X.blankCells([res], [ca], [s], blankingType=0, delta=0.)
res = C.addVars([res,celln[0]])
test.testA([res],4)
