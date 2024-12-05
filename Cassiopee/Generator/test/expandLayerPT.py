# - expandLayer (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

s = D.circle((0.,0.,0.),1.,N=100)
o = G.octree([s], [0.1], dfar=1., balancing=1)
o2 = G.expandLayer(o, level=0)
C.convertPyTree2File(o2, 'out.cgns')
