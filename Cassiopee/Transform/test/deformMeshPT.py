# - deformMesh (pyTree) -
import Transform.PyTree as T
import Converter.PyTree as C
import Geom.PyTree as D

a1 = D.sphere6((0,0,0),1,20)
a1 = C.convertArray2Tetra(a1); a1 = T.join(a1)
point = C.getValue(a1, 'GridCoordinates', 0)
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

m = D.sphere6((0,0,0),2,20)
m = T.deformMesh(m, a1)
C.convertPyTree2File(m, "out.cgns")
