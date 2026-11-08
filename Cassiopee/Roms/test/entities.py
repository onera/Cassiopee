# test Driver entities
import Roms.Driver as D

TEST = 'line'

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
    line4 = line1.Rotate('line4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [line1, line2, line3, line4])

elif TEST == 'polyline':
    # polyline
    polyline1 = T1.Polyline('polyline1', [P1, P2, P3])
    polyline2 = polyline1.Shift('polyline2', vector=(0,3,0))
    polyline3 = polyline1.Symm('polyline3', plane='yz', center=(0,0,0))
    polyline4 = polyline1.Rotate('polyline4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [polyline1,polyline2,polyline3,polyline4])

elif TEST == 'spline1':
    spline1 = T1.Spline1('spline1', [P1, P2, P3])
    spline2 = spline1.Shift('spline2', vector=(0,3,0))
    spline3 = spline1.Symm('spline3', plane='yz', center=(0,0,0))
    spline4 = spline1.Rotate('spline4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [spline1,spline2,spline3,spline4])

elif TEST == 'spline2':
    spline1 = T1.Spline2('spline1', [P1, P2, P3])
    spline2 = spline1.Shift('spline2', vector=(0,3,0))
    spline3 = spline1.Symm('spline3', plane='yz', center=(0,0,0))
    spline4 = spline1.Rotate('spline4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [spline1,spline2,spline3,spline4])

elif TEST == 'spline3':
    import Geom, Generator
    naca = Geom.naca(12, N=51)
    bbox = Generator.bbox(naca)
    grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(3,2,1))
    T1.Eq(grid1.P[1][0][0].y, -0.1*v1)
    T1.Eq(grid1.P[1][1][0].y, +0.1*v1)
    spline1 = T1.Spline3('spline1', grid1, mesh=naca)
    spline2 = spline1.Shift('spline2', (0,2,0))
    spline3 = spline1.Symm('spline1', plane='yz', center=(0,0,0))
    spline4 = spline1.Rotate('spline4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [spline1,spline2,spline3,spline4], h=0.01)

elif TEST == 'circle':
    circle1 = T1.Circle('circle1', P2, v1)
    circle2 = circle1.Shift('circle2', vector=(0,3,0))
    circle3 = circle1.Symm('circle3', plane='yz', center=(0,0,0))
    circle4 = circle1.Rotate('circle4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [circle1,circle2,circle3,circle4])

elif TEST == 'arc':
    arc1 = T1.Arc('arc1', P1, P2, P3)
    arc2 = arc1.Shift('arc2', vector=(0,3,0))
    arc3 = arc1.Symm('arc3', plane='yz', center=(0,0,0))
    arc4 = arc1.Rotate('arc4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [arc1,arc2,arc3,arc4])

elif TEST == 'superellipse':
    superellipse1 = T1.SuperEllipse('superellipse1', P2, v1, 1.)
    superellipse2 = superellipse1.Shift('superellipse2', vector=(0,3,0))
    superellipse3 = superellipse1.Symm('superellipse3', plane='yz', center=(0,0,0))
    superellipse4 = superellipse1.Rotate('superellipse4', center=(0,0,0), alpha=-90.)
    sketch1 = T1.Sketch('sketch1', [superellipse1,superellipse2,superellipse3,superellipse4])

elif TEST == 'naca':
    naca1 = T1.Naca('naca1', M=v1, P=v1, e=12.)
    # naca canot be shifted or rotated
    sketch1 = T1.Sketch('sketch1', [naca1], h=0.01)

elif TEST == 'airfoil':
    airfoil1 = T1.Airfoil('airfoil1', 'ONERA/ONERAOA213AIRFOIL')
    # airfoil canot be shifted or rotated
    sketch1 = T1.Sketch('sketch1', [airfoil1], h=0.01)

T1.solve()
T1.instantiate({'v1': 1.})

D.DRIVER.plot()
