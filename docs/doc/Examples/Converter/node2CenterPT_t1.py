# - node2Center (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

def F(x,y): return 2*x+y

# Sur une zone structuree + sur des champs
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'Helio', F, ['CoordinateX','CoordinateZ'])
a = C.initVars(a, 'centers:cellN', 1.)
a = C.node2Center(a, ['Density', 'Helio'])
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
test.testT(t, 1)

# Sur une zone structuree (creation)
b = G.cart((0,0,0), (1,1,1), (10,10,1))
b = C.addVars(b, ['Density', 'H', 'Helio'])
b = C.initVars(b, 'centers:cellN', 1.)
b = C.node2Center(b)
t = C.newPyTree(['Base',2]); t[2][1][2].append(b)
test.testT(t, 2)

# Sur un arbre
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addVars(a, ['Density', 'H', 'Helio'])
b = G.cart((10,0,0), (1,1,1), (10,10,10))
b = C.addVars(b, ['Density', 'H', 'Helio'])
t = C.newPyTree(['Base',3]); t[2][1][2] += [a,b]
t = C.initVars(t, 'centers:cellN', 1.)
t = C.node2Center(t)
test.testT(t, 3)

# Sur une liste de zones structurees
A = C.node2Center([a,b], ['Density'])
t = C.newPyTree(['Base',3]); t[2][1][2] += A
test.testT(t, 4)

# Sur une zone non structuree
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = C.addVars(a, ['Density', 'H', 'Helio'])
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'Helio', F, ['CoordinateX','CoordinateY'])
a = C.node2Center(a)
a = C.initVars(a, 'centers:cellN', 1.)
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
test.testT(t, 51)

# Sur une zone non structuree + sur des champs
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
a = C.addVars(a, ['Density', 'H', 'Helio'])
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'Helio', F, ['CoordinateX','CoordinateY'])
a = C.node2Center(a, ['Density', 'Helio'])
a = C.initVars(a, 'centers:cellN', 1.)
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
test.testT(t, 5)

# Sur une zone NGON
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.addVars(a, ['Density', 'H', 'Helio'])
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.initVars(a, 'Helio', F, ['CoordinateX','CoordinateY'])
a = C.node2Center(a, ['Density', 'Helio'])
a = C.initVars(a, 'centers:cellN', 1.)
t = C.newPyTree(['Base',3]); t[2][1][2].append(a)
test.testT(t,6)
