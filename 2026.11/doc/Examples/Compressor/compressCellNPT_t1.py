# - compressCellN (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,11,12))
C._initVars(a, '{centers:cellN}=1.')
Compressor._compressCellN(a)
Compressor._uncompressAll(a)
test.testT(a, 1)
