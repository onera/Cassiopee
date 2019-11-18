import Transform.PyTree as T
import Generator.PyTree as G
import Converter.PyTree as C
import Intersector.PyTree as XOR
import Converter.Internal as I
import time
import os, sys
import KCore.test as test

t1 = G.cart((0.,0.,0.), (0.1,0.1,0.1), (20,20,20))
#C.convertPyTree2File(t1, 'out.cgns')

t2 = T.translate(t1, (.43,.56,.27))
T._rotate(t2, (0.,0.,0.), (0.,10.,0.))

t3 = T.translate(t2, (.43,.56,.27))
T._rotate(t3, (0.,0.,0.), (0.,10.,0.))

z1s = I.getZones(t1)
z2s = I.getZones(t2)
z3s = I.getZones(t3)

t = C.newPyTree(['Base1', z1s[0], 'Base2', z2s[0], 'Base3', z3s[0]])

priorities = []
priorities.append((0,1))
priorities.append((1,2))
priorities.append((0,2))

t0 = time.time()
t = XOR.unify(t, priorities)
t1 = time.time()
print(' - UNIFY CPU time : ',t1-t0,'s')
test.testT(t,1)

priorities = []
priorities.append((1,0))
priorities.append((1,2))
priorities.append((2,0))

t0 = time.time()
t = XOR.unify(t, priorities)
t1 = time.time()
print(' - UNIFY CPU time : ',t1-t0,'s')
test.testT(t,2)

priorities = []
priorities.append((1,0))
priorities.append((2,1))
priorities.append((2,0))

t0 = time.time()
t = XOR.unify(t, priorities)
t1 = time.time()
print(' - UNIFY CPU time : ',t1-t0,'s')
test.testT(t,3)
