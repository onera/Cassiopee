# - setDegeneratedBC (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
a = G.cylinder((0,0,0), 0., 1., 360., 0., 1, (21,21,21))
t = C.newPyTree(['Base',a])
t = X.setDegeneratedBC(t)
C.convertPyTree2File(t, "out.cgns")
