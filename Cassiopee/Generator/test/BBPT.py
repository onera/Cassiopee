# - BB (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

s = D.circle((0,0,0), 1., N=100)
a = G.BB(s); a[0] = 'bbox'
C.convertPyTree2File([s,a], 'out.cgns')
