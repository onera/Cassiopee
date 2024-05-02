# - blankCellsTetra (array) - 'NODE IN'
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import Post.PyTree as P
import KCore.test as test
import Transform.PyTree as T

# Test 1
# Tet mask
m = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
m = P.exteriorFaces(m)
m = C.convertArray2Tetra(m)
m = T.reorder(m, (-1,))
#C.convertPyTree2File(m, 'm.plt')
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
t = C.newPyTree(['Cart',a])
# celln init
C._initVars(t, 'nodes:cellN', 1.)
#C.convertPyTree2File(t, 'b.plt')
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="node_in", tol=1.e-12)
#C.convertPyTree2File(t, 'out1t.cgns')
test.testT(t,1)

# Test 2
# Tet mask
m = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
m = P.exteriorFaces(m)
m = C.convertArray2Tetra(m)
m = T.reorder(m, (-1,))
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
#C.convertPyTree2File(a, 'bgm.plt')
t = C.newPyTree(['Cart',a])
# celln init
C._initVars(t, 'centers:cellN', 1.)
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="center_in", tol=1.e-12)
#C.convertPyTree2File(t, 'out2t.cgns')
test.testT(t,2)

# Test 3
# Tet mask
m = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
m = P.exteriorFaces(m)
m = C.convertArray2Tetra(m)
m = T.reorder(m, (-1,))
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
t = C.newPyTree(['Cart',a])
# celln init
C._initVars(t, 'centers:cellN', 1.)
#C.convertPyTree2File(t, 'b.plt')
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="cell_intersect", tol=1.e-12)
#C.convertPyTree2File(t, 'out3t.cgns')
test.testT(t,3)

# Test 4
# Tri mask
m = D.sphere((0,0,0), 20., 50)
m = C.convertArray2Tetra(m)
m= G.close(m) ## IMPORTANT : sinon des triangles degeneres sont engendres aux poles
#C.convertPyTree2File(m, 'before.mesh')
#C.convertPyTree2File(m, 'after.mesh')
#C.convertPyTree2File(m, 'm.plt')
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
t = C.newPyTree(['Cart',a])
# celln init
C._initVars(t, 'nodes:cellN', 1.)
#C.convertPyTree2File(t, 'b.plt')
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="node_in", tol=1.e-12)
#C.convertPyTree2File(t, 'out4t.cgns')
test.testT(t,4)

# Test 5
# Tet mask
m = D.sphere((0,0,0), 20., 50)
m = C.convertArray2Tetra(m)
m= G.close(m) ## IMPORTANT : sinon des triangles degeneres sont engendres aux poles
#C.convertPyTree2File(m, 'sph.plt')
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
t = C.newPyTree(['Cart',a])
# celln init
C._initVars(t, 'centers:cellN', 1.)
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="center_in", tol=1.e-12)
#C.convertPyTree2File(t, 'out5t.cgns')
test.testT(t,5)

# Test 6
# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))
t = C.newPyTree(['Cart',a])
#C.convertPyTree2File(t, 'bgm.plt')
# celln init
C._initVars(t, 'centers:cellN', 1.)
# Blanking
t = X.blankCellsTri(t, [[m]], [], blankingType="cell_intersect", tol=1.e-12)
#C.convertPyTree2File(t, 'out6t.cgns')
test.testT(t,6)
