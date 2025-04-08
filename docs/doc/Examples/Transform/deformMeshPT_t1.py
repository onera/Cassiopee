# - deformMesh(pyTree) -
# tests en structure
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

a1 = D.sphere6((0,0,0),1,N=18)
a1 = C.convertArray2Tetra(a1); a1 = T.join(a1)
point = (C.getValue(a1, 'CoordinateX', 0),C.getValue(a1, 'CoordinateY', 0),
         C.getValue(a1, 'CoordinateZ', 0))
a2 = T.deformPoint(a1, point, (0.1,0.05,0.2), 0.5, 2.)
delta = C.diffArrays(a2,a1)
deltax = C.getField('DCoordinateX',delta)
deltay = C.getField('DCoordinateY',delta)
deltaz = C.getField('DCoordinateZ',delta)
for noz in range(len(deltax)):
    deltax[noz][0] = 'dx'
    deltay[noz][0] = 'dy'
    deltaz[noz][0] = 'dz'
a1 = C.setFields(deltax,a1,'nodes')
a1 = C.setFields(deltay,a1,'nodes')
a1 = C.setFields(deltaz,a1,'nodes')

# 2D paroi
m = D.sphere6((0,0,0),1,N=18)
m = T.deformMesh(m, a1)
test.testT(m,1)

# 2D
m = D.sphere6((0,0,0),2,N=18)
m2 = T.deformMesh(m, a1)
test.testT(m2,2)

# 3D
m = D.sphere6((0,0,0),2,N=18)
dk = G.cart((0,0,0),(0.1,1,1),(3,1,1))
m = G.addNormalLayers(m, dk, niter=100)
m2 = T.deformMesh(m, a1)
test.testT(m2,3)
