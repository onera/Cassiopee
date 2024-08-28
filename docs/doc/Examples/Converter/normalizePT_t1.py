# - normalize (pyTree) -
# Aux noeuds
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

def f1(x1,x2,x3): return 3*x1*x2 + x3

# Structure
a = D.sphere((0,0,0), 1., 50)
a = C.initVars(a, 'sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.normalize(a, ['sx','sy','sz'])
test.testT(a,1)

# Non structure
a = D.sphere( (0,0,0), 1., 50 )
a = C.convertArray2Hexa(a)
a = C.initVars(a, 'sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.initVars(a, 'sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
a = C.normalize(a, ['sx','sy','sz'])
test.testT(a,2)

# On lists
b = D.sphere( (2.,0,0), 1., 50 )
b = C.initVars(b, 'sx', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'sy', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
b = C.initVars(b, 'sz', f1, ['CoordinateX','CoordinateY','CoordinateZ'])
out = C.normalize([a,b], ['sx','sy','sz'])
test.testT(out, 3)

# arbre
t = C.newPyTree(['Base1',2,'Base2',2])
t[2][1][2].append(b); t[2][2][2].append(b)
t = C.normalize(t, ['sx','sy','sz'])
test.testT(t,4)
