# - addVars (array) -
import Converter as C
import Generator as G
import KCore.test as test

# STRUCT
a = G.cart((0,0,0), (1,1,1), (10,10,11))
b = C.array('cell', a[2], a[3], a[4])
c = C.array('t,u', a[2], a[3], a[4])
f = C.addVars([a, b, c])
test.testA([f], 1)

# TETRA
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
ne = a[2][0].shape[0]
np = a[1][0].shape[0]
b = C.array('cell', np, ne, a[3])
c = C.array('t,u', np, ne, a[3])
f = C.addVars([a, b, c])
test.testA([f], 2)

# TETRA centers
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.node2Center(a)
ne = a[2][0].shape[0]
np = a[1][0].shape[0]
b = C.array('cell', np, ne, a[3])
c = C.array('t,u', np, ne, a[3])
f = C.addVars([a, b, c])
test.testA([f], 21)

# On lists
a = G.cart((0,0,0), (1,1,1), (10,10,11))
b = G.cart((0,0,0), (1,1,1), (5,5,2))
c = [a,b]
d = [C.array('cell', a[2], a[3], a[4]),C.array('cell', b[2], b[3], b[4])]
f = C.addVars([c,d])
test.testA(f, 3)

# NGON
a = G.cartNGon((0,0,0),(1,1,1),(11,11,11))
b = C.addVars(a, 'F')
test.testA([b], 4)

f = ['f',a[1][0:1,:],a[2],a[3]]
b = C.addVars([a,f])
test.testA([b], 5)
