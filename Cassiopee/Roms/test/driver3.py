# driver: parametric sketch from lines and spline
import Roms.Driver as D

# Create parameter
hauteur = D.Scalar('hauteur')
hauteur.range = [0., 1., 0.1]

# Create points
P1 = D.Point('P1', (0,0,0))

P2 = D.Point('P2', (1,0,0))

P3 = D.Point('P3', (1.5,1,0))
D.Eq(P3.y, hauteur)

P4 = D.Point('P4', (2.5,1,0))
D.Eq(P4.y, hauteur)

P5 = D.Point('P5', (3,0,0))

P6 = D.Point('P6', (4,0,0))

P7 = D.Point('P7', (4,2,0))

P8 = D.Point('P8', (0,2,0))

# Create lines
spline1 = D.Spline1('spline1', [P1,P2,P3,P4,P5,P6])

line1 = D.Line('line1', P6, P7)
line2 = D.Line('line2', P7, P8)
line3 = D.Line('line3', P8, P1)

# Create sketch
sketch1 = D.Sketch('sketch1', [spline1, line1, line2, line3], h=[0.01,0.01,0.01])

# solve
D.DRIVER.solve()

D.DRIVER.instantiate({'hauteur': 0.5})
sketch1.writeCAD('out.step')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'hauteur': 0.3+i/100.})
    mesh = sketch1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)
