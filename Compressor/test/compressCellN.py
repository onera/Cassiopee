# - compressCellN (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(a, '{centers:cellN}=1.')
Compressor._compressCellN(a)
C.convertPyTree2File(a, 'out.cgns')
