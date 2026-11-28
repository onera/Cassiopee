# strongly non linear equations

import Roms.Driver as D
import numpy, sympy

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., numpy.pi, 0.1]

v2 = T1.Scalar('v2', 0.)
v2.range = [0., numpy.pi, 0.1]

T1.Eq(sympy.sqrt(v1), v2*v2)

P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2')
T1.Eq(P2.x, sympy.cos(v1))
T1.Eq(P2.y, sympy.sin(v1))

P3 = T1.Point('P3')
T1.Eq(P3.x, sympy.cos(v2))
T1.Eq(P3.y, sympy.sin(v2))

line1 = T1.Line('line1', P1, P2)
line2 = T1.Line('line2', P1, P3)
line3 = T1.Line('line3', P2, P3)

sketch1 = T1.Sketch('sketch1', [line1,line3,line2])

T1.solve()
#T1.freeParams = [v1]
T1.instantiate({'v1': 0.})

D.DRIVER.plot()
