# - normL2 (pyTree) -
import Converter.PyTree as C
import KCore.test as test
import Generator.PyTree as G

# Compute L2 norm on structured mesh for a variable F(at nodes) and G(at centers), with celln
a = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a = C.initVars(a, "F", 1.)
a = C.initVars(a, "centers:G", 1.)
a = C.initVars(a, "centers:cellN", 1.)
norm = C.normL2(a, "F")
test.testO(norm, 1)
norm = C.normL2(a, "centers:G")
test.testO(norm, 11)

# Compute L2 norm on unstructured mesh for a variable F(at nodes) and G(at centers), with celln
a = G.cartTetra( (0,0,0), (0.1,0.1,0.1), (10,10,10) )
a = C.initVars(a, "F", 1.)
a = C.initVars(a, "centers:G", 1.)
a = C.initVars(a, "centers:cellN", 1.)
norm = C.normL2(a, "F")
test.testO(norm, 2)
norm = C.normL2(a, "centers:G")
test.testO(norm, 21)

# On lists
a1 = G.cart( (0,0,0), (1,1,1), (11,11,11) )
a2 = G.cart( (0,0,0), (1,1,1), (11,11,11) )
l = [a1,a2]
l[0] = C.initVars(l[0], "F", 1.)
l[1] = C.initVars(l[1], "F", 1.2)
norm = C.normL2(l, "F")
test.testO(norm, 3)
