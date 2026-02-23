# - buildMaskFiles (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.elsAProfile as elsAProfile

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
t = C.newPyTree(['Cart',a])
C._initVars(t, 'centers:cellN=({centers:CoordinateX}>0.)')
t = X.cellN2OversetHoles(t)
C.convertPyTree2File(t, "in.cgns")
tp = elsAProfile.buildMaskFiles(t, keepOversetHoles=False, prefixBase=True)
C.convertPyTree2File(tp, 'out.cgns')
