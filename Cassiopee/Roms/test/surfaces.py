# test Driver surfaces
import Roms.Driver as D

TEST = 'loft4'

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (0,0,0))
T1.Eq(P2.x, v1)

P3 = T1.Point('P3', (0,0,0))
T1.Eq(P3.x, v1)
T1.Eq(P3.y, v1)

P4 = T1.Point('P4', (0,0,0))
T1.Eq(P4.y, v1)

if TEST == 'loft1':
    # loft avec une entite
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    sketch1 = T1.Sketch('sketch1', [superellipse1])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,2), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,4), rotAxis=(0,0,1), rotAngle=50.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'loft2':
    # loft avec 4 entites
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,2), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,4), rotAxis=(0,0,1), rotAngle=50.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'loft3':
    # loft avec des guides
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,2), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,4), rotAxis=(0,0,1), rotAngle=50.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'loft4':
    # loft avec un trou
    circle1 = T1.Circle('circle1', (0,0,0), v1)
    circle2 = T1.Circle('circle2', (0,0,0), 0.1)
    sketch1 = T1.Sketch('sketch1', [circle1,circle2])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,2), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,4), rotAxis=(0,0,1), rotAngle=50.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3], close=False)

T1.solve()
T1.instantiate({'v1': 1.})

D.DRIVER.plot()
