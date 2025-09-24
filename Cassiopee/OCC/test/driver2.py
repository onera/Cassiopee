# test driver
import OCC.Driver as D

# Create parameters
hauteur = D.Scalar(1., name='hauteur')
hauteur.range = [0,1]

largeur = D.Scalar(1., name='largeur')
largeur.range = [0,2]

# equation
D.Eq(largeur.s, 2*hauteur.s)

# Create point
P1 = D.Point( (0,0,0), name='P1' )
P1.x.range = [-1,1]
P1.y.range = [-1,1]
D.Eq(P1.x.s, 0.)
D.Eq(P1.y.s, 0.)

P2 = D.Point( (1,0,0), name='P2' )
P2.x.range = [-5,5]

D.Eq(P2.x.s, P1.x.s + largeur.s)

P3 = D.Point( (1,1,0), name='P3' )
P3.x.range = [-5,5]
P3.y.range = [-5,5]

P4 = D.Point( (0,1,0), name='P4' )
P4.y.range = [-5,5]

D.Eq(P3.y.s, P4.y.s)
D.Eq(P3.y.s, P1.y.s + hauteur.s)

# Create lines
line1 = D.Line( P1, P2, name='line1' )
#line1.eq = 'length(line1) = 1.'
#D.Eq('length(line1) = 1.')

line2 = D.Line( P2, P3, name='line2' )
line3 = D.Line( P3, P4, name='line3' )
line4 = D.Line( P4, P1, name='line4' )

# Create sketch
sketch1 = D.Sketch([line1, line2, line3, line4], name='sketch1')
#sketch1.eq = 'area(sketch) = 1.'
#sketch1.eq = 'length(line3) = 2.'
#sketch1.eq = 'angle(line1, line2) = 45.'
#sketch1.update()
#sketch1.writeCAD('out.step')
#sketch1.print()

# solve
solution, freevars = D.DRIVER.solve2()
D.DRIVER.instantiate({'P3.x': 3, 'P4.y': 1.})

sketch1.writeCAD('out.step')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'P3.x': 3+i/50., 'P4.y': 1.})
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)
