# - adaptOctree (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

s = D.circle((0,0,0), 1., N=100); snear = 0.1
o = G.octree([s], [snear], dfar=5.,balancing=1)
o = C.initVars(o, 'centers:indicator', 1.)
res = G.adaptOctree(o)
C.convertPyTree2File(res, 'out.cgns')
