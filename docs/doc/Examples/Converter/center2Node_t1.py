# - center2Node (array) -
import Converter as C
import Generator as G
import KCore.test as test

def F(x,y):
    return 2*x+y
def H(x,y):
    if (x+y > 5): return 0
    else: return 1

# Structure 3D
ni = 30; nj = 40; nk = 10
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'Density',  F, ['x','y'])
a = C.initVars(a, 'cellN',  H, ['x','y'])
an = C.center2Node(a)
test.testA([an], 1)

# Test sur des listes
ni = 20; nj = 10; nk = 5
b = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
b = C.initVars(b, 'ro',  F, ['x','y'])
b = C.initVars(b, 'cellN',  H, ['x','y'])
B = C.center2Node([a,b])
test.testA(B, 2)

# Structured 2D
ni = 20; nj = 10; nk = 1
b = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
b = C.initVars(b, 'ro',  F, ['x','y'])
b = C.center2Node(b)
test.testA(B, 3)

# Structured 1D
ni = 20; nj = 1; nk = 1
b = G.cart((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
b = C.initVars(b, 'ro',  F, ['x','y'])
b = C.center2Node(b)
test.testA(B, 4)
