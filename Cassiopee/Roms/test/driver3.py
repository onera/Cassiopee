# driver: parametric sketch from lines and spline
import Roms.Driver as D

T1 = D.Part("Part1")

# Create parameter
hauteur = T1.Scalar('hauteur')
hauteur.range = [0., 1., 0.1]

# Create points
P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (1,0,0))

P3 = T1.Point('P3', (1.5,1,0))
T1.Eq(P3.y, hauteur)

P4 = T1.Point('P4', (2.5,1,0))
D.Eq(P4.y, hauteur)

P5 = T1.Point('P5', (3,0,0))

P6 = T1.Point('P6', (4,0,0))

P7 = T1.Point('P7', (4,2,0))

P8 = T1.Point('P8', (0,2,0))

# Create lines
spline1 = T1.Spline1('spline1', [P1,P2,P3,P4,P5,P6])

line1 = T1.Line('line1', P6, P7)
line2 = T1.Line('line2', P7, P8)
line3 = T1.Line('line3', P8, P1)

# Create sketch
sketch1 = T1.Sketch('sketch1', [spline1, line1, line2, line3], h=[0.01,0.01,0.01])

# solve
T1.solve()

T1.instantiate({'hauteur': 0.5})
sketch1.writeCAD('out.step')

import CPlot, time
for i in range(50):
    T1.instantiate({'hauteur': 0.3+i/100.})
    mesh = sketch1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)
