# driver: parametric loft surface
import OCC.Driver as D
import Geom
import Converter
import Generator

# Create parameter
epaisseur = D.Scalar(12., name='epaisseur')
epaisseur.range = [10, 15]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create grid
grid1 = D.Grid(bbox[0:3], bbox[3:], N=(3,3,1))
D.Eq(epaisseur.s, grid1.P[1][2][0].y.s)

# Create profile
spline1 = D.Spline3(grid1, mesh=naca, name='spline1')

# Create sketch 1
sketch1 = D.Sketch([spline1], name='sketch1')

# Create sketch 2
sketch2 = D.Sketch([spline1], name='sketch2')
sketch2.position.z.v = 1.
#sketch2.rotAngle.v = 1. # cause trouble
sketch2.update()

# surface
surface1 = D.Loft([sketch1, sketch2], name='surface1')

# test
D.DRIVER.solve2()

D.DRIVER.instantiate({'epaisseur': 0.8})

surface1.writeCAD('out.step')

#mesh = sketch1.mesh(0.01, 0.01, 0.01)
#mesh += sketch2.mesh(0.01, 0.01, 0.01)
mesh = surface1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(surface1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'epaisseur': 0.3+i/50.})
    mesh = surface1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)