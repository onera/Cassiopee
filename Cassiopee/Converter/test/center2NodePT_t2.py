# - center2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+y

def H(x,y):
    return (x+y <= 5).astype(float)

# center2Node: cree une nouvelle zone
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'], isVectorized=True)
a = C.initVars(a, 'cellN', H, ['CoordinateX','CoordinateY'], isVectorized=True)
t = C.newPyTree(['Base1',3,a])
t = C.center2Node(t)
test.testT(t, 1)

# center2Node: modifie une variable (structure)
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.node2Center(a,'GridCoordinates')
a = C.initVars(a, 'centers:Density', F, ['CoordinateX','CoordinateY'], isVectorized=True)
a = C.initVars(a, 'centers:cellN', H, ['CoordinateX','CoordinateY'], isVectorized=True)
a = C.rmVars(a,['centers:CoordinateX','centers:CoordinateY','centers:CoordinateZ'])
t = C.newPyTree(['Base1',3,a])
t = C.center2Node(t, 'cellN')
test.testT(t, 2)

# center2Node: modifie une variable (NGon)
ni = 11; nj = 11; nk = 2
a = G.cartNGon((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'centers:Density', F, ['centers:CoordinateX','centers:CoordinateY'], isVectorized=True)
a = C.center2Node(a, 'centers:Density')
a = C.rmVars(a, 'centers:Density')
test.testT(a, 3)
