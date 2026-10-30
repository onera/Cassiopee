import Roms.Driver as D

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (0,0,0))
T1.Eq(P2.x, v1)

P3 = T1.Point('P3', (0,0,0))
T1.Eq(P3.x, v1)
T1.Eq(P3.y, v1)

superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
sketch1 = T1.Sketch('sketch1', [superellipse1])
sketch1.set(position=(0,0,0), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=10.)

sketch2 = sketch1.Copy('sketch2')
sketch2.set(position=(1,0,0))

T1.solve()
T1.instantiate({'v1': 1.})

D.DRIVER.plot()
