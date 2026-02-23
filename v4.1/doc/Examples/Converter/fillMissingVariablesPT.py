# - fillMissingVariables (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,11)); a[0] = 'cart1'
b = G.cart((1,0,0), (2,1,1), (10,10,11)); b[0] = 'cart2'
C._addVars(a, 'rou'); C._addVars(a, 'rov')
C._addVars(a, 'centers:cellN')
C._addVars(a, ['Density', 'Hx', 'centers:Hy'])

t = C.newPyTree(['Base',a,b])
t = C.fillMissingVariables(t)
C.convertPyTree2File(t, 'out.cgns')
