# - octree2Struct (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s],[snear], dfar=5., balancing=1)
res = G.octree2Struct(res, vmin=5, ext=2, optimized=1)
C.convertArrays2File([s]+res, "out.plt")
