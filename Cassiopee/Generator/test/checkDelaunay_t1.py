# - checkDelaunay (array) -
import Converter as C
import Generator as G
import Transform as T
import Geom as D
import KCore.test as test
A = D.text1D('STEPHANIE', font='text1')
A = C.convertArray2Tetra(A); a = T.join(A)
# Triangulation respecting given contour
tri = G.constrainedDelaunay(a)
res = G.checkDelaunay(a, tri)
test.testA([res])
