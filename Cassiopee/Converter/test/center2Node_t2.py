# - center2Node (array) -
import Converter as C
import Generator as G
import KCore.test as test

def F(x,y):
    return 2*x+y

# Non structure TETRA
ni = 30; nj = 40; nk = 12
a = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'ro', F, ['x','y'])
ac = C.node2Center(a)
an = C.center2Node(ac)
test.testA([an], 1)

# Test sur une liste
ni = 10; nj = 15; nk = 2
b = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
b = C.initVars(b, 'ro', F, ['x','y'])
bc = C.node2Center(b)
B = C.center2Node([ac,bc])
test.testA(B, 2)

# Non structure TRI
ni = 30; nj = 40; nk = 1
a = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'ro', F, ['x','y'])
ac = C.node2Center(a)
an = C.center2Node(ac)
test.testA([an], 3)

# Non structure BAR
ni = 30; nj = 1; nk = 1
a = G.cartTetra((0,0,0), (10./(ni-1),1,1), (ni,nj,nk))
a = C.initVars(a, 'ro', F, ['x','y'])
ac = C.node2Center(a)
an = C.center2Node(ac)
test.testA([an], 4)

# Non structure NGON
ni = 30; nj = 40; nk = 12
a = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'ro', F, ['x','y'])
ac = C.node2Center(a)
an = C.center2Node(ac)
test.testA([an], 5)
