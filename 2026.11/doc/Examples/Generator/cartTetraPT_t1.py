# - cartTetra (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# BAR
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,1,1)) # BAR - dir x
test.testT(a, 1)

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (1,10,1)) # BAR - dir y
test.testT(a, 11)

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (1,1,10)) # BAR - dir z
test.testT(a, 12)

a = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,1,1))
b = G.cartTetra((0.,0.,0.), (1.,1.,1.), (1,2,1))
c = G.cartTetra((0.,0.,0.), (1.,1.,1.), (1,1,2))
a = C.newPyTree(['Base',[a,b,c]]) # mono-cell in each directions
test.testT(a, 13)

# TRI
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1)) # TRI - dir x/y
test.testT(a, 2)

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,1,10)) # TRI - dir x/z
test.testT(a, 21)

a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (1,10,10)) # TRI - dir y/z
test.testT(a, 22)

a = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,2,1))
b = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,1,2))
c = G.cartTetra((0.,0.,0.), (1.,1.,1.), (1,2,2))
a = C.newPyTree(['Base',[a,b,c]]) # mono-cell in each directions
test.testT(a, 23)

# TETRA
a = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,2,2))
test.testT(a, 3)
test.writeCoverage(100)
