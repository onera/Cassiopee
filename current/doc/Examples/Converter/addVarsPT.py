# - addVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,11))
a = C.addVars(a, 'rou')
a = C.addVars(a, 'centers:cellN')
a = C.addVars(a, ['Density', 'Hx', 'centers:Hy'])
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, 'out.cgns')
