# - signNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartNGon((0,0,0), (1,1,1), (3,5,7))
a = C.makeParentElements(a)
test.testT(a, 1)
