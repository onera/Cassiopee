# - expandLayer (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0.,0.,0.), 1., N=100)
o = G.octree([s], [0.1], dfar=1., balancing=1)
o2 = G.expandLayer(o, level=0)
C.convertArrays2File([o, o2], "out.plt")
