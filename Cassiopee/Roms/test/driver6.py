# driver: parametric loft surface
import Roms.Driver as D
import Geom
import Converter
import Generator

# Create parameter
epaisseur = D.Scalar('epaisseur')
epaisseur.range = [10, 15, 1.]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create grid
grid1 = D.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
D.Eq(epaisseur, grid1.P[1][2][0].y)

# Create profile
spline1 = D.Spline3('spline1', grid1, mesh=naca)

# Create sketch 1
sketch1 = D.Sketch('sketch1', [spline1])

# Create sketch 2
sketch2 = D.Sketch('sketch2', [spline1])
sketch2.position.z.v = 1.
sketch2.update()

# surface
surface1 = D.Loft('surface1', [sketch1, sketch2])

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