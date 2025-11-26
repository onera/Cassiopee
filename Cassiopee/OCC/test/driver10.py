# driver: parametric boolean
import OCC.Driver as D
import Geom
import Converter
import Generator

# Create parameter
length = D.Scalar('length', 12.)
length.range = [1, 15, 1]

# Create sketch 1
circle1 = D.Circle('circle1', (0,0,0), 1.)
sketch1 = D.Sketch('sketch1', [circle1])

# Create sketch 2
circle2 = D.Circle('circle2', (0,0,5), 1.)
sketch2 = D.Sketch('sketch2', [circle2])

# surface1
surface1 = D.Loft('surface1', [sketch1, sketch2])

# Create sketch 3
circle3 = D.Circle('circle3', (0,0,0), 0.5)
sketch3 = D.Sketch('sketch3', [circle3])

# Create sketch 4
circle4 = D.Circle('circle4', (0,0,1), 0.5)
sketch4 = D.Sketch('sketch4', [circle4])
D.Eq(circle4.P[0].z.s, length.s)

# surface2
surface2 = D.Loft('surface2', [sketch3, sketch4])
surface2.rotAxis.x.v = 0
surface2.rotAxis.y.v = 1
surface2.rotAxis.z.v = 0
surface2.rotAngle.v = 90.
surface2.position.z.v = 2.

# surface finale
#surface = D.Merge('surface', listSurfaces=[surface1,surface2])
surface = D.Union('surface', listSurfaces1=[surface1], listSurfaces2=[surface2])

# test
D.DRIVER.solve2()
D.DRIVER.instantiate({'length': 10})

surface.writeCAD('out.step')

import CPlot, time
pt = D.DRIVER.walkDOE()
while pt is not None:
    D.DRIVER.instantiate(pt)
    mesh = surface.mesh(0.1, 0.1, 0.1)
    CPlot.display(mesh)
    pt = D.DRIVER.walkDOE()
    time.sleep(0.5)
