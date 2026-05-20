# parametric 2D nozzle
import Roms.Driver as D

#==============
# parameters
#==============

# throat height
Ht = D.Scalar('Ht', 0.5)
Ht.range = [0.1, 2., 0.1]

# inlet height
Hi = D.Scalar('Hi', 1.)
Hi.range = [0.1, 2., 0.1]

# exit height
He = D.Scalar('He', 1.)
He.range = [0.1, 2., 0.1]

# convergent length
Lc = D.Scalar('Lc', 2.)
Lc.range = [0.1,3]

# divergent length
Ld = D.Scalar('Ld', 2.)
Ld.range = [0.1,3]

# convergent curvature
Rc = D.Scalar('Rc', 1.)

# divergent curvature
Rd = D.Scalar('Rd', 1.)

#=======================
# create points
#=======================

# Create points P1 (throat)
P1 = D.Point('P1')
P1.x.range = [-1,1]
D.Eq(P1.y, Ht)

# Create points P2 (inlet)
P2 = D.Point('P2')
D.Eq(P2.x, -Lc)
D.Eq(P2.y, Hi)

# Create points P3 (exit)
P3 = D.Point('P3')
D.Eq(P3.x, Ld)
D.Eq(P3.y, He)

#=================
# create entities
#=================
# basic with lines
line1 = D.Line('line1', P2, P1)
line2 = D.Line('line2', P1, P3)

# create symmetric part
P1S = D.Point('P1S')
D.Eq(P1S.x, P1.x)
D.Eq(P1S.y, -P1.y)
P2S = D.Point('P2S')
D.Eq(P2S.x, P2.x)
D.Eq(P2S.y, -P2.y)
P3S = D.Point('P3S')
D.Eq(P3S.x, P3.x)
D.Eq(P3S.y, -P3.y)
line3 = D.Line('line3', P2S, P1S)
line4 = D.Line('line4', P1S, P3S)

line5 = D.Line('line5', P2, P2S)
line6 = D.Line('line6', P3, P3S)

sketch1 = D.Sketch('sketch1', [line1,line2,line5,line4,line3,line6])

#=======================
# solve and instantiate
#=======================
solution, freevars = D.DRIVER.solve()

D.DRIVER.instantiate({'P1.x':0.,
                      'Hi': 1.,
                      'Ht': 0.2,
                      'He': 1.,
                      'Lc': 2.,
                      'Ld': 2.})
sketch1.writeCAD('out.step')

sketch1.print()

D.DRIVER.plot(sketch1)
