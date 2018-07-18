# - blankCells (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D

surf = D.sphere((0,0,0), 0.5, 20)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
t = C.newPyTree(['Cart',a])
C._initVars(t, 'centers:cellN', 1.)

bodies = [[surf]]
# Matrice de masquage (arbre d'assemblage)
import numpy
BM = numpy.array([ [ 1] ] )

t = X.blankCells(t, bodies, BM, blankingType='cell_intersect', delta=0.)
C.convertPyTree2File(t, 'out.cgns')
