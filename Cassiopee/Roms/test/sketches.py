# test sketches
import Roms.Driver as D

TEST = 'airfoil'

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (0,0,0))
T1.Eq(P2.x, v1)

P3 = T1.Point('P3', (0,0,0))
T1.Eq(P3.x, v1)
T1.Eq(P3.y, v1)

if TEST == 'line':
    line1 = T1.Line('line1', P1, P3)
    sketch1 = T1.Sketch('sketch1', [line1])
    
elif TEST == 'superellipse':
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    sketch1 = T1.Sketch('sketch1', [superellipse1])

elif TEST == 'spline1':
    spline1 = T1.Spline1('spline1', [P1, P2, P3])
    sketch1 = T1.Sketch('sketch1', [spline1])

elif TEST == 'airfoil':
    airfoil1 = T1.Airfoil('airfoil1', 'ONERA/ONERAOA213AIRFOIL')
    sketch1 = T1.Sketch('sketch1', [airfoil1], h=0.01)

sketch1.set(position=(0,0,0), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=10.)

sketch2 = sketch1.Copy('sketch2')
sketch2.set(position=(1,0,0))

sketch3 = sketch1.Copy('sketch3')
sketch3.set(scaleCenter=(0,0,0), scaleFactor=0.5)

T1.solve()
T1.instantiate({'v1': 1.})

D.DRIVER.plot()
