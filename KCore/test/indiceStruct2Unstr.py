# - indiceStruct2Unstr (array) -
import KCore
import Geom as D
import Converter as C
import Transform as T
import Generator as G
import numpy

# Maillage structure de 6 blocks
a = D.sphere6((0,0,0), 1., N=20)
# Maillage QUAD equivalent
b = C.convertArray2Hexa(a); b = T.join(b); b = G.close(b)

# Indices des vertex du maillage structure dont on cherche la correspondance
npts = a[0][1].shape[1]
indicesS = numpy.arange(0, npts, dtype=numpy.int32)

# Liste des indices des vertex correspondants dans le maillage QUAD b
indicesU = KCore.indiceStruct2Unstr(a[0], b, indicesS, 1.e-14)
print indicesU

# Liste de tous les indices correspondants
indicesU = KCore.indiceStruct2Unstr2(a, b, 1.e-14)
print indicesU
