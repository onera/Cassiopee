# parametric fin
import Roms.Driver as D
import Geom
import Generator

# height of fin
height = D.Scalar('height', 3.)
height.range = [1.0, 10.0, 1.]

# length of fin at base
baseLength = D.Scalar('baseLength', 3.)
D.Eq(baseLength, 3.)

# thickness of fin at base
baseThickness = D.Scalar('baseThickness', 0.5)
D.Eq(baseThickness, 0.5)

# length of thin at top
topLength = D.Scalar('topLength', 2.)
D.Eq(topLength, 1.)

# thickness of fin at top
topThickness = D.Scalar('topThickness', 0.3)
D.Eq(topThickness, 0.3)

# shift of top
#finShift = D.Scalar('finShift', 0.1)

# Base
profile = Geom.naca(12, N=51)
#profile = Geom.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
bbox = Generator.bbox(profile)
grid1 = D.Grid('grid1', bbox[0:3], bbox[3:], N=(2,2,1))
D.Eq(grid1.P[1][0][0].x, baseLength)
D.Eq(grid1.P[1][1][0].x, baseLength)
D.Eq(grid1.P[0][1][0].y, 0.5*baseThickness)
D.Eq(grid1.P[1][1][0].y, 0.5*baseThickness)
D.Eq(grid1.P[0][0][0].y, -0.5*baseThickness)
D.Eq(grid1.P[1][0][0].y, -0.5*baseThickness)
spline1 = D.Spline3('spline1', grid1, mesh=profile)
sketch1 = D.Sketch('sketch1', [spline1])

# Top
profile = Geom.naca(12, N=51)
#naca = Geom.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
bbox = Generator.bbox(profile)
grid2 = D.Grid('grid2', bbox[0:3], bbox[3:], N=(2,2,1))
D.Eq(grid2.P[1][0][0].x, topLength)
D.Eq(grid2.P[1][1][0].x, topLength)
D.Eq(grid2.P[0][1][0].y, 0.5*topThickness)
D.Eq(grid2.P[1][1][0].y, 0.5*topThickness)
D.Eq(grid2.P[0][0][0].y, -0.5*topThickness)
D.Eq(grid2.P[1][0][0].y, -0.5*topThickness)
spline2 = D.Spline3('spline2', grid2, mesh=profile)
sketch2 = D.Sketch('sketch2', [spline2])
D.Eq(sketch2.position.x, baseLength-topLength)
D.Eq(sketch2.position.z, height)

P1 = D.Point('P1')
P2 = D.Point('P2')
D.Eq(P2.z, height)
D.Eq(P2.x, baseLength-topLength)

P3 = D.Point('P3')
D.Eq(P3.x, baseLength)
P4 = D.Point('P4')
D.Eq(P4.x, baseLength)
D.Eq(P4.z, height)

line1 = D.Line('line1', P1, P2)
line2 = D.Line('line2', P3, P4)

sketch3 = D.Sketch('sketch3', [line2])

# body Surface
surface1 = D.Loft('surface1', listSketches=[sketch1, sketch2], listGuides=[])
#surface1 = D.MergeEdges('surface1', listSketches=[sketch1, sketch2, sketch3])

# solve
solution, freevars = D.DRIVER.solve()
D.DRIVER.instantiate({'height': 1.0})

#sketch2.writeCAD('out.step')
surface1.writeCAD('out.step')
#sketch2.writeCAD('out.step')
