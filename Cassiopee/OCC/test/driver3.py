# driver: parametric sketch from lines and spline
import OCC.Driver as D

# Create parameter
hauteur = D.Scalar(1., name='hauteur')
hauteur.range = [0,1]

# Create points
P1 = D.Point( (0,0,0), name='P1' )

P2 = D.Point( (1,0,0), name='P2' )

P3 = D.Point( (1.5,1,0), name='P3' )
D.Eq(P3.y.s, hauteur.s)

P4 = D.Point( (2.5,1,0), name='P4' )
D.Eq(P4.y.s, hauteur.s)

P5 = D.Point( (3,0,0), name='P5' )

P6 = D.Point( (4,0,0), name='P6' )

P7 = D.Point( (4,2,0), name='P7' )

P8 = D.Point( (0,2,0), name='P8' )

# Create lines
spline1 = D.Spline1( [P1,P2,P3,P4,P5,P6], name='spline1' )

line1 = D.Line( P6, P7, name='line1' )
line2 = D.Line( P7, P8, name='line2' )
line3 = D.Line( P8, P1, name='line3' )

# Create sketch
sketch1 = D.Sketch([spline1, line1, line2, line3], name='sketch1')

# solve
D.DRIVER.solve2()
D.DRIVER.instantiate({'hauteur': 0.5})

sketch1.writeCAD('out.step')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'hauteur': 0.3+i/100.})
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)