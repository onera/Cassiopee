import Roms.Driver as D

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

P1 = T1.Point('P1', (0,0,0))
T1.Eq(P1.x, v1)

P2 = T1.Point('P2', (0,0,0))
T1.Eq(P2.y, v1)

P3 = P1.Shift('P3', vector=(1,0,0))

P4 = P3.Symm('P4', plane='yz', center=(0,0,0))

P5 = P1.Rotate('P5', center=(0,0,0), alpha=90.)

T1.solve()
T1.instantiate({'v1': 1.})

P1.print()
P2.print()
P3.print()
P4.print()
P5.print()

#D.DRIVER.plot()
