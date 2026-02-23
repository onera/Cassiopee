# - conformOctree3 (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0,0,0), 1., N=100); snear = 1.
res = G.octree([s], [snear], dfar=5., balancing=1, ratio=3)
res = G.conformOctree3(res)
C.convertArrays2File([res], 'out.plt')

#s = D.sphere((0,0,0), 1., N=100); snear = 100.
#res = G.octree([s], [snear], dfar=5., balancing=1, ratio=3)
#res = G.conformOctree3(res)
#C.convertArrays2File([res], 'out.plt')
