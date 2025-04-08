# - setValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Array structure
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
a = C.initVars(a, 'centers:Density', 1)
C.setValue(a, 'CoordinateX', (10,1,1), 0.25)
C.setValue(a, 'centers:Density', (10,1,1), 0.5)
C.setValue(a, 'GridCoordinates', (11,1,1), [0.3,0.2,0.1])
t = C.newPyTree(['Base',a])
test.testT(t, 1)

# Array non structure
Ni = 40; Nj = 50; Nk = 20
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
a = C.initVars(a, 'centers:Density', 1)
C.setValue(a, 'CoordinateX', 9, 0.1 )
C.setValue(a, 'centers:Density', 8, 0.5)
t = C.newPyTree(['Base',a])
test.testT(t, 2)
