# driver: parametric sketch from lines
import Roms.Driver as D

# Create parameters
hauteur = D.Scalar('hauteur')
hauteur.range = [0., 1., 0.1]

largeur = D.Scalar('largeur')
largeur.range = [0., 2., 0.1]

# equation
D.Eq(largeur, 2*hauteur)

# Create point
P1 = D.Point('P1', (0,0,0))
D.Eq(P1.x, 0.)
D.Eq(P1.y, 0.)

P2 = D.Point('P2', (1,0,0))
D.Eq(P2.x, P1.x + largeur)

P3 = D.Point('P3', (1,1,0))
P3.x.range = [-5., 5.]
P4 = D.Point('P4', (0,1,0))
D.Eq(P3.y, P4.y)
D.Eq(P3.y, P1.y + hauteur)

# Create lines
line1 = D.Line('line1', P1, P2)
line2 = D.Line('line2', P2, P3)
line3 = D.Line('line3', P3, P4)
line4 = D.Line('line4', P4, P1)

# Create sketch
sketch1 = D.Sketch('sketch1', [line1, line2, line3, line4], h=[0.01,0.01,0.01])
#sketch1.eq = 'area(sketch) = 1.'
#sketch1.eq = 'length(line3) = 2.'
#sketch1.eq = 'angle(line1, line2) = 45.'
#sketch1.update()
#sketch1.writeCAD('out.step')
#sketch1.print()

# solve
solution, freevars = D.DRIVER.solve()

D.DRIVER.instantiate({'P3.x': 3, 'hauteur': 1.})

sketch1.writeCAD('out.step')

import CPlot, time
for i in range(50):
    D.DRIVER.instantiate({'P3.x': 3+i/50., 'hauteur': 1.})
    mesh = sketch1.mesh()
    CPlot.display(mesh)
    time.sleep(0.5)
