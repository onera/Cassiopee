# driver: parametric loft surface
import Roms.Driver as D
import Geom
import Converter
import Generator

T1 = D.Part("Part1")

# Create parameter
epaisseur = T1.Scalar('epaisseur')
epaisseur.range = [10, 15, 1.]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create grid
grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(3,3,1))
T1.Eq(epaisseur, grid1.P[1][2][0].y)

# Create profile
spline1 = T1.Spline3('spline1', grid1, mesh=naca)

# Create sketch 1
sketch1 = T1.Sketch('sketch1', [spline1])

# Create sketch 2
sketch2 = T1.Sketch('sketch2', [spline1])
sketch2.position.z.v = 1.

# surface
surface1 = T1.Loft('surface1', [sketch1, sketch2], h=[0.01,0.01,0.01])

# test
T1.solve()

T1.instantiate({'epaisseur': 0.8})

surface1.writeCAD('out.step')

#mesh = sketch1.mesh()
#mesh += sketch2.mesh()
mesh = surface1.mesh()
#D.DRIVER._dXdmu(surface1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

import CPlot, time
for i in range(50):
    T1.instantiate({'epaisseur': 0.3+i/50.})
    mesh = surface1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)
