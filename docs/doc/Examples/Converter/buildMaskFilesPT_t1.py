# - blankCells (pyTree) -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import Converter.elsAProfile as elsAProfile
import KCore.test as test
LOCAL = test.getLocal()

surf = D.sphere((0,0,0), 0.5, 20)

a = G.cart((-1.,-1.,-1.),(0.1,0.1,0.1), (20,20,20))
t = C.newPyTree(['Cart']); t[2][1][2].append(a)
t = C.initVars(t, 'centers:cellN', 1.)

bodies = [[surf]]
# Matrice de masquage (arbre d'assemblage)
import numpy
BM = numpy.array([ [ 1] ] )

t = X.blankCells(t, bodies, BM, blankingType='cell_intersect', delta=0.)
t = X.cellN2OversetHoles(t)
tp = elsAProfile.buildMaskFiles(t,keepOversetHoles=False,fileDir=LOCAL)
test.testT(tp, 1)
