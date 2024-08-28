# - checkDelaunay (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test
A = D.text1D('STEPHANIE', font='text1')
A = C.convertArray2Tetra(A); a = T.join(A)
a = C.addVars(a, 'Density'); a = C.addVars(a, 'centers:F')

# Triangulation respecting given contour
tri = G.constrainedDelaunay(a)
tri = C.addVars(tri, 'Density'); tri = C.addVars(tri, 'centers:F')
res = G.checkDelaunay(a, tri)
test.testT([res])
