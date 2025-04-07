# - adaptOctree (array) -
import Generator as G
import Converter as C
import Geom as D

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=5., balancing=1)
indic = C.node2Center(o)
indic = C.initVars(indic, 'indicator', 1.)
res = G.adaptOctree(o, indic)
C.convertArrays2File([o,res], "out.plt")
