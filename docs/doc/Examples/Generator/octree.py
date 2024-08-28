# - octree (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0,0,0), 1., N=100); snear = 0.01
res = G.octree([s], [snear], dfar=5., balancing=2)
C.convertArrays2File([res], 'out.plt')
