# - compressIndices (pyTree) -
import Compressor.compressor as Co
import Generator.PyTree as G
import Converter.Internal as Internal
import numpy as np

a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
n = Internal.getNodeFromName(a,"ElementConnectivity")[1]
print(f"Original data size : {n.size*4}")
comp = Co.compressIndices((8,n))
print(f"Compressed data size : {comp[2].size}")
print(f"Ratio de compression : {comp[2].size/(n.size*4.)}")
n2 = Co.uncompressIndices(comp)[0]
diff = n - n2
err = np.max(np.abs(diff))
print(f"Erreur max faite sur les indices : {err}")
