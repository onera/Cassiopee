# - conformOctree3 (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

s = D.circle((0,0,0), 1., N=100); snear = 1.
res = G.octree([s], [snear], dfar=5., balancing=1, ratio=3)
res = G.conformOctree3(res)
C.convertPyTree2File(res, 'out.cgns')
