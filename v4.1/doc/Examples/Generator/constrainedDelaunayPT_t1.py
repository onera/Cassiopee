# - constrainedDelaunay (pyTree) -
import Converter
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

#---------
# 1er test
#---------
# generation du contour
a = Converter.array('x,y,z', 4, 4, 'BAR')

a[1][0][1] = 1.
a[1][0][2] = 1.
a[1][1][2] = 1.
a[1][1][3] = 1.

for i in range(4):
    a[2][0][i] = i+1
    a[2][1][i] = i+2

a[2][1][3] = 1
z = C.convertArrays2ZoneNode('contour', [a])

def dens(x,y): return 3*x*y

z = C.initVars(z,'centers:cellN',1.); z = C.initVars(z, 'Density', dens, ['CoordinateX','CoordinateY'])

# constrained Delaunay
tri = G.constrainedDelaunay( z )

test.testT(tri, 1)
