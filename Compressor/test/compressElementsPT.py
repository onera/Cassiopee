# - compressElements (pyTree) -
import Compressor.PyTree as Compressor
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartHexa((0,0,0), (1,1,1), (25,23,24))
Compressor._compressElements(a)
#Compressor._uncompressAll(a)
C.convertPyTree2File(a, 'out.cgns')
