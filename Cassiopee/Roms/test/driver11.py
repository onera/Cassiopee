# driver: parametric sketch from lines and spline
# + volume meshing
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
sketch1 = D.Sketch('sketch1', [spline1, line1, line2, line3], h=[0.1,0.1,0.1])

# Create volume
vol1 = D.Volume2D('vol1', [sketch1])

# solve
D.DRIVER.solve()

D.DRIVER.instantiate({'hauteur': 0.5})
m = vol1.MeshAsReference()
#D.DRIVER.instantiate({'hauteur': 0.8})
#m = vol1.Dmesh()

import CPlot.PyTree as CPlot
import time
pt = D.DRIVER.walkDOE()
while pt is not None:
    #m = vol1.Mesh()
    #m = sketch1.Mesh()
    m = vol1.Dmesh()
    CPlot.display(m)
    pt = D.DRIVER.walkDOE()
    time.sleep(0.5)