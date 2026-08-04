# example of parametric sketch from lines
import Roms.Driver as D

T1 = D.Part('part1')

# Create parameters
hauteur = T1.Scalar('hauteur')
hauteur.range = [0., 1., 0.1]

largeur = T1.Scalar('largeur')
largeur.range = [0., 2., 0.1]

# equation
# equations can use standard math functions for instance: D.sympy.cos, D.sympy.exp, ...
T1.Eq(largeur, 2*hauteur)

# Create points
P1 = T1.Point('P1', (0,0,0))
T1.Eq(P1.x, 0.)
T1.Eq(P1.y, 0.)

P2 = T1.Point('P2', (1,0,0))
T1.Eq(P2.x, P1.x + largeur)

P3 = T1.Point('P3', (1,1,0))
P3.x.range = [-5., 5.]
P4 = T1.Point('P4', (0,1,0))
T1.Eq(P3.y, P4.y)
T1.Eq(P3.y, P1.y + hauteur)

# Create lines
line1 = T1.Line('line1', P1, P2)
line2 = T1.Line('line2', P2, P3)
line3 = T1.Line('line3', P3, P4)
line4 = T1.Line('line4', P4, P1)

# Create sketch
# a sketch gathers entities that should form a single contiguous wire (open or closed)
sketch1 = T1.Sketch('sketch1', [line1, line2, line3, line4], h=[0.01,0.01,0.01])

# solve
freevars = T1.solve()
print("freevars=", freevars)

# instantiate. You must set every variables in freevars
T1.instantiate({'P3.x': 3, 'hauteur': 1.})

#sketch1.writeCAD('out.step')
T1.writeCAD('out.step')

# display in image
import CPlot
mesh = sketch1.mesh()
CPlot.display(mesh, export="image.png", offscreen=2)
CPlot.finalizeExport(2)
import os; os._exit(0)
