# test Driver surfaces
import Roms.Driver as D

TEST = 'sweep1'

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (0,0,0))
T1.Eq(P2.x, v1)

P3 = T1.Point('P3', (0,0,0))
T1.Eq(P3.x, v1)
T1.Eq(P3.y, v1*0.5)

P4 = T1.Point('P4', (0,0,0))
T1.Eq(P4.y, v1*0.5)

if TEST == 'loft1':
    # loft continue avec un seul edge
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    sketch1 = T1.Sketch('sketch1', [superellipse1])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=40.)

    sketch4 = sketch1.Copy('sketch4')
    sketch4.set(position=(0,0,6), scaleCenter=(0,0,0), scaleFactor=0.5)

    sketch5 = sketch1.Copy('sketch5')
    sketch5.set(position=(0,0,8))

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3,sketch4,sketch5])

elif TEST == 'loft2':
    # loft continue avec un sketch de 4 edges
    #line1 = T1.Line('line1', P1, P2)
    #line2 = T1.Line('line2', P2, P3)
    #line3 = T1.Line('line3', P3, P4)
    #line4 = T1.Line('line4', P4, P1)
    line1 = T1.Spline1('line1', [P1, P2]) # gordon needs spline and not line
    line2 = T1.Spline1('line2', [P2, P3])
    line3 = T1.Spline1('line3', [P3, P4])
    line4 = T1.Spline1('line4', [P4, P1])

    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,2), rotAxis=(0,0,1), rotAngle=20.)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), rotCenter=(0,0,4), rotAxis=(0,0,1), rotAngle=50.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'loft3':
    # loft continue deux sketchs avec 1 edge
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    sketch1 = T1.Sketch('sketch1', [superellipse1])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=20.)

    surface1 = T1.Loft('surface1', [sketch1,sketch2])

elif TEST == 'loft4':
    # loft continue deux sketchs avec 4 edges
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), rotCenter=(0,0,0), rotAxis=(0,0,1), rotAngle=v1)

    surface1 = T1.Loft('surface1', [sketch1,sketch2])

elif TEST == 'loft5':
    # loft continue deux sketchs avec 2 edges
    spline1 = T1.Spline1('spline1', [P1,P2,P3,P4])
    spline2 = spline1.Symm(center=(0,0,0), plane='yz')
    sketch1 = T1.Sketch('sketch1', [spline1,spline2])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), scaleCenter=(0,0,0), scaleFactor=0.6)

    surface1 = T1.Loft('surface1', [sketch1,sketch2])

elif TEST == 'loft6':
    # loft continue deux sketchs avec 1 circle
    circle1 = T1.Circle('circle1', (0,0,0), v1)
    sketch1 = T1.Sketch('sketch1', [circle1])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), scaleCenter=(0,0,0), scaleFactor=0.6)

    surface1 = T1.Loft('surface1', [sketch1,sketch2])

elif TEST == 'loft7':
    # loft linear 3 sketchs avec 1 circle
    circle1 = T1.Circle('circle1', (0,0,0), v1)
    sketch1 = T1.Sketch('sketch1', [circle1])

    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(position=(0,0,2), scaleCenter=(0,0,0), scaleFactor=0.6)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), scaleCenter=(0,0,0), scaleFactor=0.8)

    surface1 = T1.LoftSet('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'loft8':
    # loft linear 3 sketchs non isomorphe
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    circle1 = T1.Circle('circle1', (0,0,0), v1)
    sketch2 = T1.Sketch('sketch2', [circle1])
    sketch2.set(position=(0,0,2), scaleCenter=(0,0,0), scaleFactor=0.6)

    sketch3 = sketch1.Copy('sketch3')
    sketch3.set(position=(0,0,4), scaleCenter=(0,0,0), scaleFactor=0.8)

    surface1 = T1.LoftSet('surface1', [sketch1,sketch2,sketch3])

elif TEST == 'extrude1':
    # extrude following vector
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    surface1 = T1.Extrude('surface1', sketch1, vector=(0,0,1))

elif TEST == 'sweep1':
    # sweep following path
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])

    spline1 = T1.Spline1('spline1', [(0,0,0), (0,0,1), (0,2,2)])
    sketch2 = T1.Sketch('sketch2', [spline1])

    surface1 = T1.Sweep('surface1', sketch1, sketch2)

T1.solve()
T1.instantiate({'v1': 1.})

D.DRIVER.plot()
