# driver: parametric revolve surface
import Roms.Driver as D
import Converter.PyTree as C

# Create parameter
epaisseur = D.Scalar('epaisseur')
epaisseur.range = [0.5, 2., 0.1]

# create Points
P1 = D.Point('P1', (0.1,0,0))
P2 = D.Point('P2', (1.,0.,1))
P2.x.range = [0.5, 2., 0.1]
P3 = D.Point('P3', (0.1,0.,2))

D.Eq(P2.x, epaisseur)
D.Eq(P1.y, epaisseur*0.2)

# Create profile
spline1 = D.Spline1('spline1', [P1,P2,P3])

# Create sketch 1
sketch1 = D.Sketch('sketch1', [spline1])

# surface
surface1 = D.Revolve('surface1', sketch1, center=(0,0,0), axis=(0,0,1), angle=90., h=[0.05,0.05,0.1])

# test
D.DRIVER.solve()

D.DRIVER.instantiate({'epaisseur': 1.2})
#surface1.writeCAD('out.step')
#D.DRIVER._dXdmu(surface1, m)

m = surface1.MeshAsReference()
#m = surface1.Mesh()
#C.convertPyTree2File(m, 'out.cgns')

import CPlot.PyTree as CPlot, time
point = D.DRIVER.walkDOE()
while point is not None:
    D.DRIVER.instantiate(point)
    #m = surface1.Mesh()
    m = surface1.Dmesh()
    CPlot.display(m)
    point = D.DRIVER.walkDOE()
    time.sleep(0.5)