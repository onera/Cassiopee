# - applyBCOverlaps (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (30,30,10))
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmin')
t = C.newPyTree(['Base',a])
t = X.applyBCOverlaps(t)
C.convertPyTree2File(t, 'out.cgns')
