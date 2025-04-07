# - constrainedDelaunay (array) -
import Converter as C
import Generator as G
import Transform as T
import Geom as D

A = D.text1D('CASSIOPEE')
A = C.convertArray2Tetra(A); a = T.join(A)
# Triangulation respecting given contour
tri = G.constrainedDelaunay(a)
C.convertArrays2File([a,tri], "out.plt")
