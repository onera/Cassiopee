# driver: parametric revolve surface
import Roms.Driver as D
import Converter.PyTree as C

T1 = D.Part("Part1")

# Create parameter
epaisseur = T1.Scalar('epaisseur')
epaisseur.range = [0.5, 2., 0.1]

# create Points
P1 = T1.Point('P1', (0.1,0,0))
P2 = T1.Point('P2', (1.,0.,1))
P2.x.range = [0.5, 2., 0.1]
P3 = T1.Point('P3', (0.1,0.,2))

T1.Eq(P2.x, epaisseur)
T1.Eq(P1.y, epaisseur*0.2)

# Create profile
spline1 = T1.Spline1('spline1', [P1,P2,P3])

# Create sketch 1
sketch1 = T1.Sketch('sketch1', [spline1])

# surface
surface1 = T1.Revolve('surface1', sketch1, center=(0,0,0), axis=(0,0,1), angle=90., h=[0.05,0.05,0.1])

# test
T1.solve()

T1.instantiate({'epaisseur': 1.2})
#surface1.writeCAD('out.step')
#D.DRIVER._dXdmu(surface1, m)

m = surface1.MeshAsReference()
#m = surface1.Mesh()
#C.convertPyTree2File(m, 'out.cgns')

import CPlot.PyTree as CPlot, time
point = T1.walkDOE()
while point is not None:
    T1.instantiate(point)
    #m = surface1.Mesh()
    m = surface1.Dmesh()
    CPlot.display(m)
    point = T1.walkDOE()
    time.sleep(0.5)