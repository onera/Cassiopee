# - compressCoords (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Compressor._compressCoords(a, tol=1.e-7)
C.convertPyTree2File(a, 'out.cgns')
