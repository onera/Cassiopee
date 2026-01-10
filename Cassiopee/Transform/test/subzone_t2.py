# - subzone (array) -
import Transform as T
import Generator as G
import KCore.test as test

# type=nodes, with different output dimensions

# TETRA -> TETRA (same as default, dimOut=-1)
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [10, 19, 20, 220], type='nodes', dimOut=3)
test.testA([b], 1)

# TETRA -> TRI
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 20, 220], type='nodes', dimOut=2)
test.testA([b], 2)

# TETRA -> NODE
a = G.cartTetra((0,0,0), (1,1,1), (10,20,10))
b = T.subzone(a, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 20, 220], type='nodes', dimOut=0)
test.testA([b], 3)

# TRI -> BAR
a = G.cartTetra((0,0,0), (1,1,1), (10,20,1))
b = T.subzone(a, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 20], type='nodes', dimOut=1)
test.testA([b], 4)

# QUAD -> BAR
a = G.cartHexa((0,0,0), (1,1,1), (10,20,1))
b = T.subzone(a, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 20], type='nodes', dimOut=1)
test.testA([b], 5)
