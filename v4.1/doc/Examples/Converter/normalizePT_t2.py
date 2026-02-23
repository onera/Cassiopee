# - normalize (pyTree) -
# aux centres
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

def f1(x1,x2,x3): return 3*x1*x2 + x3

# Structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.node2Center(a, 'GridCoordinates')
a = C.initVars(a, 'centers:sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'centers:sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'centers:sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.normalize(a, ['centers:sx','centers:sy','centers:sz'])
test.testT(a, 1)

# Non structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.convertArray2Hexa(a)
a = C.node2Center(a, 'GridCoordinates')
a = C.initVars(a, 'centers:sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'centers:sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'centers:sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.normalize(a,['centers:sx','centers:sy','centers:sz'])
test.testT(a,2)

# On lists
b = D.sphere((2.,0,0), 1., 50 )
b = C.node2Center(b, 'GridCoordinates')
b = C.initVars(b, 'centers:sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'centers:sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'centers:sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.normalize(b,['centers:sx','centers:sy','centers:sz'])
test.testT(a,3)

# arbre
t = C.newPyTree(['Base1',2,'Base2',2])
t[2][1][2].append(b); t[2][2][2].append(b)
t = C.normalize(t,['centers:sx','centers:sy','centers:sz'])
test.testT(t,4)
