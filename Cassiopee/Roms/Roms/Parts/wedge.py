# parametric wedge
import Roms.Driver as D

#==============
# parameters
#==============

# forward length
Lf = D.Scalar('Lf', 2.)
Lf.range = [0.1, 3.]

# backward length
Lb = D.Scalar('Lb', 0.1)
Lb.range = [0.1, 1.]

# top height
Ht = D.Scalar('Ht', 0.1)
Ht.range = [0.1, 1.]

# bottom height
Hb = D.Scalar('Hb', 0.1)
Hb.range = [0.1, 1.]

# side length
Ls = D.Scalar('Ls', 1.)
Ls.range = [0.1, 2.]

#=======================
# create points
#=======================

# forward point
P0 = D.Point('P0')
D.Eq(P0.x, -Lf)

P1 = D.Point('P1')
D.Eq(P1.y, -Ls)

P2 = D.Point('P2')
D.Eq(P2.y, +Ls)

P3 = D.Point('P3')
D.Eq(P3.z, Ht)

P4 = D.Point('P4')
D.Eq(P4.z, -Hb)

P5 = D.Point('P5')
D.Eq(P5.x, +Lb)

#=================
# create entities
#=================

surface1 = D.FillLinear('surface1', [P0,P3,P1])
surface2 = D.FillLinear('surface2', [P0,P3,P2])
surface3 = D.FillLinear('surface3', [P0,P4,P1])
surface4 = D.FillLinear('surface4', [P0,P4,P2])

surface5 = D.FillLinear('surface5', [P5,P3,P1])
surface6 = D.FillLinear('surface6', [P5,P3,P2])
surface7 = D.FillLinear('surface7', [P5,P4,P1])
surface8 = D.FillLinear('surface8', [P5,P4,P2])

surface = D.Merge('surface', [surface1,surface2,surface3,surface4,surface5,surface6,surface7, surface8], h=(0.1,0.1,0.01))

#=======================
# solve and instantiate
#=======================
solution, freevars = D.DRIVER.solve()

D.DRIVER.instantiate({'Lf': 2,
                      'Lb': 0.1,
                      'Ht': 0.1,
                      'Hb': 0.1,
                      'Ls': 1.})

surface.writeCAD('out.step')
D.DRIVER.plot(surface)
