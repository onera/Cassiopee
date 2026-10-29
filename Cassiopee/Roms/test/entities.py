import Roms.Driver as D

TEST = 'polyline'

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
    # line
    line1 = T1.Line('line1', P1, P3)
    line2 = line1.Shift('line2', vector=(0,5,0))
    line3 = line1.Symm('lines3', plane='yz', center=(0,0,0)) 
    
    line1.print()
    line2.print()
    
    # sketch for viz
    sketch1 = T1.Sketch('sketch1', [line1, line2, line3])

elif TEST == 'polyline':
    # polyline
    polyline1 = T1.Polyline('polyline1', [P1, P2, P3])
    polyline2 = polyline1.Shift('polyline2', vector=(0,3,0))
    # sketch for viz
    sketch1 = T1.Sketch('sketch1', [polyline1,polyline2])

elif TEST == 'spline1':
    spline1 = T1.Spline1('spline1', [P1, P2, P3])
    spline2 = spline1.Shift('spline2', vector=(0,3,0))
    spline3 = spline1.Symm('spline3', plane='yz', center=(0,0,0))
    sketch1 = T1.Sketch('sketch1', [spline1,spline2,spline3])

elif TEST == 'circle':
    circle1 = T1.Circle('circle1', P2, v1)
    circle2 = circle1.Shift('circle2', vector=(0,3,0))
    circle3 = circle1.Symm('circle3', plane='yz', center=(0,0,0))
    sketch1 = T1.Sketch('sketch1', [circle1,circle2,circle3])

elif TEST == 'arc':
    arc1 = T1.Arc('arc1', P1, P2, P3)
    arc2 = arc1.Shift('arc2', vector=(0,3,0))
    arc3 = arc1.Symm('arc3', plane='yz', center=(0,0,0))
    sketch1 = T1.Sketch('sketch1', [arc1,arc2,arc3])

elif TEST == 'superellipse':
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    superellipse2 = superellipse1.Shift('superellipse2', vector=(0,3,0))
    superellipse3 = superellipse1.Symm('superellipse3', plane='yz', center=(0,0,0))
    sketch1 = T1.Sketch('sketch1', [superellipse1,superellipse2,superellipse3])

T1.solve()
T1.instantiate({'v1': 1.})



D.DRIVER.plot()
