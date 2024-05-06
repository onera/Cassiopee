# - compressCoords (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (5,3,4))
Compressor._compressCoords(a, tol=1.e-7, ctype=0)
Compressor._uncompressAll(a)
test.testT(a, 1)
