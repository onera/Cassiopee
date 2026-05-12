# parametric missile
import Roms.Driver as D

#===================
# Create parameters
#===================
radius1 = D.Scalar('radius1', 10.)
radius1.range = [5., 15., 1.]

radius2 = D.Scalar('radius2', 5.)
radius2.range = [0.1, 8., 1.]

# length back
length1 = D.Scalar('length1', 100.)
D.Eq(length1, 100.)

# length middle
length2 = D.Scalar('length2', 20.)
D.Eq(length2, 20.)

# length forward
length3 = D.Scalar('length3', 10.)
D.Eq(length3, 10.)

#=================
# Create points
#=================
P1 = D.Point('P1')
D.Eq(P1.y, radius1)

P2 = D.Point('P2')
D.Eq(P2.y, radius1)
D.Eq(P2.x, -length1)

P3 = D.Point('P3')
D.Eq(P3.y, radius2)
D.Eq(P3.x, -length1-length2)

P4 = D.Point('P4')
D.Eq(P4.x, -length1-length2-length3)
D.Eq(P4.y, P3.y/20.)

P5 = D.Point('P5')
D.Eq(P5.x, -length1-length2-length3)

line1 = D.Line('line1', P1, P2)
line2 = D.Line('line2', P2, P3)
spline1 = D.Spline1('spline1', [P3,P4,P5])

# sketch1
sketch1 = D.Sketch('sketch1', [line1, line2, spline1])

# surface1 - body
surface1 = D.Revolve('surface1', sketch1, center=(0,0,0), axis=(1,0,0), angle=360.)

# body closure
circle1 = D.Circle('circle1', (0,0,0), radius1)
sketch2 = D.Sketch('sketch2', [circle1])
surface2 = D.Fill('surface2', sketch2)

surface3 = D.Merge('surface3', [surface1, surface2])

#=======================
# solve and instantiate
#=======================
solution, freevars = D.DRIVER.solve()

D.DRIVER.instantiate({'radius1': 10., 'radius2': 4.})

surface3.writeCAD('out.step')
