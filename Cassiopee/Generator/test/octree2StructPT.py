# - octree2Struct (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

s = D.circle((0,0,0), 1., N=100); snear = 0.1
res = G.octree([s],[snear], dfar=5., balancing=1)
res = G.octree2Struct(res, vmin=5, ext=2,merged=1)
C.convertPyTree2File(res, 'out.cgns')
