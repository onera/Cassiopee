# - rotate (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C
import math
import KCore.test as test
vectnames = [['centers:VelocityX','centers:VelocityY','centers:VelocityZ']]

ni = 91; nj = 21
a = G.cylinder((0,0,0),0.5,1.,0.,90.,1.,(ni,nj,1))
a = C.initVars(a, 'centers:VelocityX', 0.)
a = C.initVars(a, 'centers:VelocityY', 0.)
a = C.initVars(a, 'centers:VelocityZ', 1.)
vx = C.getField('centers:VelocityX', a)[0]
vy = C.getField('centers:VelocityY', a)[0]
nic = vx[2]; njc = vx[3]
for j in range(njc):
    for i in range(nic):
        i0 = i * math.pi/2./90.
        vx[1][0,i+j*nic] = math.cos(i0)
        vy[1][0,i+j*nic] = math.sin(i0)
C.setFields([vx], a, loc='centers')
C.setFields([vy], a, loc='centers')
a = C.center2Node(a, 'centers:VelocityX')
# Rotate with axis+angle
b = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 90., vectors=vectnames)
test.testT(b, 1)
# Rotate with axis transformations
c = T.rotate(a, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
             ((1,1,0), (1,-1,0), (0,0,1)), vectors=vectnames)
test.testT(c, 2)

# Rotate with three angles
d = T.rotate(a, (0.,0.,0.), (0.,0.,90.), vectors=vectnames)
test.testT(d, 3)
