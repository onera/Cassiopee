# - getVolumeMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as T

# method=1
a = G.cartNGon((0,0,0),(1,1,1),(3,5,7))
a = G.getVolumeMap(a, method=1)
T.testT(a,1)
