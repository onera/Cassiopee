# - compressNGonIndices (pyTree) -
import Compressor.compressor as Co
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
n = Internal.getNodeFromName(a,"ElementConnectivity")[1]
print(n)
print(f"Original data size : {n.size*4}")
comp = Co.compressNGonIndices(n)
print(comp)
print(f"Compressed data size : {comp[1].size}")
print(f"Ratio de compression : {comp[1].size/(n.size*4.)}")
n2 = Co.uncompressNGonIndices(comp)[0]
diff = n - n2
err = numpy.max(numpy.abs(diff))
print(f"Erreur max faite sur les indices : {err}")
