# parametric 3D nozzle
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
Lc.range = [0.1, 3.]

# divergent length
Ld = D.Scalar('Ld', 2.)
Ld.range = [0.1, 3.]

# convergent curvature
Rc = D.Scalar('Rc', 1.)

# divergent curvature
Rd = D.Scalar('Rd', 1.)

#=======================
# create points
#=======================

# Create points P1 (throat)
P1 = D.Point('P1')
P1.x.range = [-2., 2.]
D.Eq(P1.y, Ht)

# Create points P2 (inlet)
P2 = D.Point('P2')
D.Eq(P2.x, -Lc)
D.Eq(P2.y, Hi)

P2a = D.Point('P2a')
D.Eq(P2a.x, -0.1*Lc)
D.Eq(P2a.y, Ht)

P2b = D.Point('P2b')
D.Eq(P2b.x, -0.5*Lc)
D.Eq(P2b.y, Ht+0.5*(Ht+Hi))

# Create points P3 (exit)
P3 = D.Point('P3')
D.Eq(P3.x, Ld)
D.Eq(P3.y, He)

P3a = D.Point('P3a')
D.Eq(P3a.x, 0.1*Ld)
D.Eq(P3a.y, Ht)

P3b = D.Point('P3b')
D.Eq(P3b.x, 0.5*Ld)
D.Eq(P3b.y, Ht+0.5*(Ht+He))

#=================
# create entities
#=================
# basic with lines
spline1 = D.Spline1('spline1', [P2, P2b, P2a, P1])
spline2 = D.Spline1('spline2', [P1, P3a, P3b, P3])

sketch1 = D.Sketch('sketch1', [spline1,spline2])
surface1 = D.Revolve('surface1', sketch1, center=(0,0,0), axis=(1,0,0), angle=360., h=(0.1,0.1,0.01))

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
surface1.writeCAD('out.step')

surface1.print()

D.DRIVER.plot(surface1)
