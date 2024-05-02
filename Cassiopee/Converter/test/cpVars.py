# - cpVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10)); a[0] = 'cart1'
b = G.cart((0,0,0),(1,1,1),(10,10,10)); b[0] = 'cart2'
C._initVars(a, 'Density', 2.)
C._initVars(b, 'centers:H', 4.)
a = C.cpVars(a, 'Density', a, 'G') # copy in a
C._cpVars(a, 'Density', b, 'Density') # copy from a to b
a = C.cpVars(b, 'centers:H', a, 'centers:H') # copy from b to a
t = C.newPyTree(['Base',a,b])
C.convertPyTree2File(t, 'out.cgns')
