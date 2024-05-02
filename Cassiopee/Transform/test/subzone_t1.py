# - subzone (array) -
import Converter as C
import Transform as T
import Generator as G
import KCore.test as test

# Structure + champ
a = G.cart((0,0,0), (1,1,1), (10,20,10))
a = C.initVars(a, 'celln', 1)
b = T.subzone(a, (3,3,3), (7,8,5))
test.testA([b], 1)

# Structure + champ / Indices negatifs
a = G.cart((0,0,0), (1,1,1), (10,20,10))
a = C.initVars(a, 'celln', 1)
b = T.subzone(a, (3,3,3), (10,19,8))
b = T.subzone(a, (3,3,3), (-1,-2,-3))
test.testA([b], 5)

# TETRA, type=nodes
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [10,19,20,220], type='nodes')
test.testA([b], 2)

# TETRA, type=elements
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [1], type='elements')
test.testA([b], 6)

# TETRA, type=faces
# le nbre de faces d'un tetra est 4
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [0*4+0+1], type='faces')
test.testA([b], 7)

# NGON, type=elements
a = G.cartNGon((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [0], type='elements')
test.testA([b], 3)

# NGON, type=faces
a = G.cartNGon((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [1], type='faces')
test.testA([b], 4)
