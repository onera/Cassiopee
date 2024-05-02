# - blankCellsTetra (pyTree) - 'NODE IN'
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G

# Tet mask
mT4 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
mT4 = C.convertArray2Tetra(mT4)

# Mesh to blank
a = G.cart((-5.,-5.,-5.), (0.5,0.5,0.5), (100,100,100))

t = C.newPyTree(['Cart',a])
C._initVars(t, 'centers:cellN', 1.)

masks = [[mT4]]
# Matrice de masquage (arbre d'assemblage)
import numpy
BM = numpy.array([[1]])

t1 = X.blankCellsTetra(t, masks, BM, blankingType='node_in', tol=1.e-12)
C.convertPyTree2File(t1, 'out.cgns')

t2 = C.convertArray2Tetra(t)
t2 = X.blankCellsTetra(t2, masks, BM, blankingType='node_in', tol=1.e-12)
C.convertPyTree2File(t2, 'out2.cgns')

t3 = C.convertArray2NGon(t)
t3 = X.blankCellsTetra(t3, masks, BM, blankingType='node_in', tol=1.e-12)
C.convertPyTree2File(t3, 'out3.cgns')
