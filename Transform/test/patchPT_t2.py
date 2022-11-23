# - patch (pyTree) - global nodes indices
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test
import numpy

def dens(x,y): return 3*x*y

# cas 2D
c1 = G.cart((0,0,0), (0.01,0.01,1), (201,101,2))
c1 = C.addBC2Zone(c1,'wall1','BCWall','imin')
c1 = C.addBC2Zone(c1,'match1','BCMatch','imax',c1,'imin',[1,2,3])
c1 = C.addBC2Zone(c1, 'overlap1', 'BCOverlap', 'jmax')
c1 = C.initVars(c1,'centers:celln',1.)
c1 = C.initVars(c1,'Density',dens,['CoordinateX','CoordinateY'])

c2 = G.cart((0,0,0), (0.01,0.01,1), (51,81,2))
t = C.newPyTree(['Base']); t[2][1][2].append(c2)
c2 = T.rotate(c2, (0,0,0),(0,0,1),0.2)
c2 = C.addBC2Zone(c2,'wall1','BCWall','imin')
c2 = C.addBC2Zone(c2,'overlap1','BCOverlap','imax')
c2 = C.initVars(c2, 'centers:celln',1.)
c2 = C.initVars(c2,'Density',dens,['CoordinateX','CoordinateY'])
# indice a partir duquel la zone est patchee
im1 = 51; jm1 = 81; km1 = 2   # dimensions de la zone patchee
im2 = 201; jm2 = 101; km2 = 2 # dimensions de la zone a patcher
# nodes: tableau des noeuds globaux correspondant au patch sur c1 en (1,1,1)
ip = 1; jp = 1; kp = 1
indp = (kp-1)*im2*jm2 + (jp-1)*im2 + ip-1
nodes = numpy.arange((im1*jm1*km1), dtype=Internal.E_NpyInt)
c = 0
for k in range(0,km1):
    for j in range(0,jm1):
        for i in range(0,im1):
            nodes[c] = indp+k*im2*jm2+j*im2+i+1; c += 1
a = T.patch(c1, c2, nodes=nodes)

t = C.newPyTree(['Base',3, a])
test.testT(a, 1)

# cas 3D
c1 = G.cart((0,0,0), (0.01,0.01,1), (201,101,20))
c1 = C.addBC2Zone(c1,'wall1','BCWall','imin')
c1 = C.addBC2Zone(c1,'match1','BCMatch','imax',c1,'imin',[1,2,3])
c1 = C.addBC2Zone(c1, 'overlap1', 'BCOverlap', 'jmax')
c1 = C.initVars(c1,'centers:celln',1.)
c1 = C.initVars(c1,'Density',dens,['CoordinateX','CoordinateY'])

c2 = G.cart((0,0,0), (0.01,0.01,1), (51,81,20))
c2 = T.rotate(c2, (0,0,0),(0,0,1),0.2)
c2 = C.addBC2Zone(c2,'wall1','BCWall','imin')
c2 = C.addBC2Zone(c2,'overlap1','BCOverlap','imax')
c2 = C.initVars(c2, 'centers:celln',1.)
c2 = C.initVars(c2, 'Density', dens, ['CoordinateX','CoordinateY'])
# indice a partir duquel la zone est patchee
im1 = 51; jm1 = 81; km1 = 20   # dimensions de la zone patchee
im2 = 201; jm2 = 101; km2 = 20 # dimensions de la zone a patcher
# nodes: tableau des noeuds globaux correspondant au patch sur c1 en (1,1,1)
ip = 1; jp = 1; kp = 1
indp = (kp-1)*im2*jm2 + (jp-1)*im2 + ip-1
nodes = numpy.arange((im1*jm1*km1),dtype=Internal.E_NpyInt)
c = 0
for k in range(0,km1):
    for j in range(0,jm1):
        for i in range(0,im1):
            nodes[c] = indp+k*im2*jm2+j*im2+i+1; c += 1
a = T.patch(c1, c2, nodes=nodes)
t = C.newPyTree(['Base', a])
test.testT(t, 2)
