# - blanCells (pyTree) -
# - cas 2D -
import Converter.PyTree as C
import Connector.PyTree as X
import Generator.PyTree as G
import Geom.PyTree as D
import KCore.test as test

surf = D.circle((0,0,0), 0.5, 20)

a = G.cart((-1.,-1.,0.),(0.1,0.1,1), (20,20,1))
C._addBC2Zone(a, 'ov', 'BCOverlap', 'jmin')
t = C.newPyTree(['Cart',2,a])
t[2][1] = C.addState(t[2][1], 'EquationDimension', 2)
C._fillEmptyBCWith(t, 'wall', 'BCWall')
C._addVars(t, 'Density')
C._initVars(t, 'centers:cellN', 1.)
bodies = [[surf]]
criterions = ['cell_intersect', 'cell_intersect_opt', 'center_in']
# Matrice de masquage (arbre d'assemblage)
import numpy
BM = numpy.array([[1]])
c = 1
for delta in [0.,0.1]:
    for type in criterions:
        if (type == 'cell_intersect_opt' and delta > 0): c += 1
        else:
            t2 = X.blankCells(t, bodies, BM,
                              blankingType=type, delta=delta,
                              dim=2)
            test.testT(t2,c)
            c += 1
