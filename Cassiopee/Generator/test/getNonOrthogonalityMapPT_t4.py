# - getNonOrthogonalityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

import Converter.Internal as Internal
import numpy

# check 45-degree non-ortho faces

# 2D unstructured tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,1), (2,2,1))
x = numpy.array([0,1+0.5,0-0.5,1], dtype=numpy.float64)
y = numpy.array([0,0+0.5,1-0.5,1], dtype=numpy.float64)
Internal.getNodeFromName(a, 'CoordinateX')[1] = x
Internal.getNodeFromName(a, 'CoordinateY')[1] = y
t = G.getNonOrthogonalityMap(a)
nonOrtho = Internal.getNodeFromName(t, 'nonOrthogonality')[1]
print('test %d: min/max = %f/%f'%(2, numpy.min(nonOrtho), numpy.max(nonOrtho)))
test.testT(t, 2)

# 2D unstructured quad
a = G.cartHexa((0.,0.,0.), (0.1,0.1,1), (3,2,1))
x = numpy.array([0,1,3,0,2,3], dtype=numpy.float64)
y = numpy.array([0,0,0,1,1,1], dtype=numpy.float64)
Internal.getNodeFromName(a, 'CoordinateX')[1] = x
Internal.getNodeFromName(a, 'CoordinateY')[1] = y
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getNonOrthogonalityMap(a)
nonOrtho = Internal.getNodeFromName(t, 'nonOrthogonality')[1]
print('test %d: min/max = %f/%f'%(3, numpy.min(nonOrtho), numpy.max(nonOrtho)))
test.testT(t, 3)

# 3D unstructured hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (3,2,2))
x = numpy.array([0,1,3,0,2,3,0,1,3,0,2,3], dtype=numpy.float64)
y = numpy.array([0,0,0,1,1,1,0,0,0,1,1,1], dtype=numpy.float64)
z = numpy.array([0,0,0,0,0,0,1,1,1,1,1,1], dtype=numpy.float64)
Internal.getNodeFromName(a, 'CoordinateX')[1] = x
Internal.getNodeFromName(a, 'CoordinateY')[1] = y
Internal.getNodeFromName(a, 'CoordinateZ')[1] = z
t = G.getNonOrthogonalityMap(a)
nonOrtho = Internal.getNodeFromName(t, 'nonOrthogonality')[1]
print('test %d: min/max = %f/%f'%(4, numpy.min(nonOrtho), numpy.max(nonOrtho)))
test.testT(t, 4)

# 3D unstructured penta
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (2,2,2))
x = numpy.array([0,1+0.5,0-0.5,1, 0,1+0.5,0-0.5,1], dtype=numpy.float64)
y = numpy.array([0,0+0.5,1-0.5,1, 0,0+0.5,1-0.5,1], dtype=numpy.float64)
z = numpy.array([0,0,0,0,1,1,1,1], dtype=numpy.float64)
Internal.getNodeFromName(a, 'CoordinateX')[1] = x
Internal.getNodeFromName(a, 'CoordinateY')[1] = y
Internal.getNodeFromName(a, 'CoordinateZ')[1] = z
t = G.getNonOrthogonalityMap(a)
nonOrtho = Internal.getNodeFromName(t, 'nonOrthogonality')[1]
print('test %d: min/max = %f/%f'%(5, numpy.min(nonOrtho), numpy.max(nonOrtho)))
test.testT(t, 5)
