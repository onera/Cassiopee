# - subzone (pyTree)-
import Converter.PyTree as C
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.test as test

# type=nodes, with different output dimensions

# 3D non structure
# TETRA -> TRI
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
t = T.subzone(t, [1, 2, 11, 12, 201, 202, 211], dimOut=2)
test.testT(t,1)

# TETRA -> NODE
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
C._initVars(a,'{Density}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base']); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 0)
t = T.subzone(t, [1, 2, 11, 12, 201, 202, 211], dimOut=0)
test.testT(t,2)

# 3D ME -> 2D ME
"""a = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 0)
t = T.subzone(t, [i for i in range(1, 401)], dimOut=2)
C.convertPyTree2File(t, "out1.cgns"); exit()
test.testT(t,10)
"""


# 2D non structure
# TRI -> BAR
a = G.cartTetra((0,0,0), (1,1,1), (10,20,1))
C._initVars(a,'{F}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 1)
t = T.subzone(t, [1, 2, 3, 4, 11, 12, 13, 14, 15], dimOut=1)
test.testT(t,3)

# TRI -> NODE
a = G.cartTetra((0,0,0), (1,1,1), (10,20,1))
C._initVars(a,'{F}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 0)
t = T.subzone(t, [1, 2, 3, 4, 11, 12, 13, 14, 15], dimOut=0)
test.testT(t,4)

# 2D ME -> BAR
"""a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b], None)
t = C.newPyTree(['Base',2]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 1)
C.convertPyTree2File(t, "out1.cgns"); exit()
#t = T.subzone(t, [1, 2, 3, 4, 11, 12, 13, 14, 15], dimOut=1)
test.testT(t,11)"""


# 1D non structure
# BAR -> NODE
a = G.cartTetra((0,0,0), (1,1,1), (10,1,1))
C._initVars(a,'{F}=3*{CoordinateX}*{CoordinateY}')
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
t[2][1] = C.addState(t[2][1], 'EquationDimension', 0)
t = T.subzone(t, [1, 2, 3, 4], dimOut=0)
test.testT(t,5)
