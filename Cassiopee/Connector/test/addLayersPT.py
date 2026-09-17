# - addLayers (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

N = 101
h = 2./(N-1)
a = G.cart((-1.,-1.,-1.), (h,h,h), (N,N,N))
t = C.newPyTree(['Cart', a])
C._initVars(t, '{centers:cellN}=(1.-({centers:CoordinateX}*{centers:CoordinateX}+{centers:CoordinateY}*{centers:CoordinateY}+{centers:CoordinateZ}*{centers:CoordinateZ}<0.25))')

X._addLayers(t, depth=1, cellNName="centers:cellN")
C.convertPyTree2File(t, 'out.cgns')
