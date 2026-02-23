# - stack (array) -
import Generator as G
import Transform as T
import Converter as C
import KCore.test as test

# Struct -> Struct
a = G.cylinder( (0,0,0), 1, 1.3, 360, 0, 1., (50,10,1))
b = T.rotate(a, (0,0,0), (1,0,0), 5.)
b = T.translate(b, (0,0,0.5))
c = G.stack(a, b)
test.testA([c], 1)

# BAR -> QUAD
a = G.cartTetra((0,0,0), (1,1,1), (50,1,1))
b = G.cartTetra((0,1,0), (1,1,1), (50,1,1))
c = G.stack(a, b)
test.testA([c], 2)

# QUAD -> HEXA
a = G.cartHexa((0,0,0), (1,1,1), (5,5,1))
b = G.cartHexa((0,0,1), (1,1,1), (5,5,1))
c = G.stack(a, b)
test.testA([c], 3)

# TRI -> HEXA
a = G.cartTetra((0,0,0), (1,1,1), (3,3,1))
b = G.cartTetra((0,0,1), (1,1,1), (3,3,1))
c = G.stack(a, b)
test.testA([c], 4)
