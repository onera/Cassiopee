# - compressFields (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateY}')
Compressor._compressFields(a, tol=1.e-6, ctype=0)
Compressor._uncompressAll(a)
test.testT(a, 1)

Compressor._compressFields(a, tol=1.e-6, ctype=0, varNames=['F', 'centers:G'])
Compressor._uncompressAll(a)
test.testT(a, 2)
