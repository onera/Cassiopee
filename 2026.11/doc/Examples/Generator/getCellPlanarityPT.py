# - getCellPlanarity (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

a = D.sphere((0,0,0), 1., 10)
a = G.getCellPlanarity(a)
C.convertPyTree2File(a, 'out.cgns')
