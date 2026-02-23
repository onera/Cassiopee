# - addVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Sur une zone
a = G.cart((0,0,0), (1,1,1), (10,10,11))
a = C.addVars(a, 'rou')
a = C.addVars(a, 'centers:cellN')
a = C.addVars(a, ['Density', 'Hx', 'centers:Hy'])
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Sur un arbre
a = G.cart((0,0,0), (1,1,1), (10,10,11))
b = G.cart((10,0,0), (1,1,1), (10,10,11))
t = C.newPyTree(['Base',a,b])
t = C.addVars(t, 'rou')
test.testT(t, 2)

# Sur une liste de zones
A = C.addVars([a,b], 'rou')
t = C.newPyTree(['Base']); t[2][1][2] += A
test.testT(t, 3)
#
# sur un NGON
a = G.cartNGon((0,0,0), (1,1,1), (10,10,11))
a = C.addVars(a, 'rou')
#a = C.addVars(a, 'centers:cellN')
t = C.newPyTree(['Base',a])
test.testT(t, 4)
