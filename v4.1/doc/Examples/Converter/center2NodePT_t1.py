# - center2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+y

def H(x,y):
    if (x+y > 5): return 0
    else: return 1

# center2Node: cree une nouvelle zone (structure)
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(a,'centers:F',1.)
C._initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX','CoordinateY'])
b = C.center2Node(a); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t, 1)

# center2Node: modifie une variable (structure)
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX','CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t, 2)

# center2Node: cree une nouvelle zone (non structure)
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(a, 'centers:F', 1.)
C._initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX','CoordinateY'])
b = C.center2Node(a,'centers:F'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t, 3)

# center2Node: modifie une variable (non structure)
ni = 30; nj = 40; nk = 2
a = G.cartTetra((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX','CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t, 4)

# center2Node: cree une nouvelle zone (non structure)
ni = 30; nj = 40; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
C._initVars(a, 'centers:F', 1.)
C._initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'cellN', H, ['CoordinateX','CoordinateY'])
b = C.center2Node(a,'centers:F'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t,5)

# center2Node: modifie une variable (non structure)
ni = 30; nj = 40; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a, 'GridCoordinates')
C._initVars(a, 'centers:Density', F, ['CoordinateX','CoordinateY'])
C._initVars(a, 'centers:cellN', H, ['CoordinateX','CoordinateY'])
a = C.rmVars(a,['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
b = C.center2Node(a, 'centers:cellN'); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
test.testT(t,6)
