# - signNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (3,3,2))
a = C.signNGonFaces(a)
test.testT(a, 1)