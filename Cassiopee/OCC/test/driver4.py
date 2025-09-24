# test driver
import OCC.Driver as D
import Geom
import Generator
import Converter

# Create parameter
epaisseur = D.Scalar(12., name='epaisseur')
epaisseur.range = [10,15]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create grid
grid1 = D.Grid(bbox[0:3], bbox[3:], N=(3,3,1))
grid1.P[1][2][0].y.range = [0, 5]
D.Eq(epaisseur.s, grid1.P[1][2][0].x.s)

# Create profile
spline1 = D.Spline3( grid1, mesh=naca, name='spline1' )

# Create sketch
sketch1 = D.Sketch([spline1], name='sketch1')

# test
D.DRIVER.solve2()

#grid1.P[1][2][0].y.print()

D.DRIVER.instantiate({'P.1.2.0.y': 0.8})

sketch1.writeCAD('out.step')

mesh = sketch1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(sketch1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'P.1.2.0.y': 0.3+i/50.})
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)