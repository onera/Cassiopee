# parametric fin
import Roms.Driver as D
import Geom
import Generator

def createPart(name):

    T1 = D.Part(name)

    # height of fin
    height = T1.Scalar('height', 3.)
    height.range = [1.0, 10.0, 1.]

    # length of fin at base
    baseLength = T1.Scalar('baseLength', 3.)
    D.Eq(baseLength, 3.)

    # thickness of fin at base
    baseThickness = T1.Scalar('baseThickness', 0.5)
    D.Eq(baseThickness, 0.5)

    # length of thin at top
    topLength = T1.Scalar('topLength', 2.)
    D.Eq(topLength, 1.)

    # thickness of fin at top
    topThickness = T1.Scalar('topThickness', 0.3)
    D.Eq(topThickness, 0.3)

    # shift of top
    #finShift = T1.Scalar('finShift', 0.1)

    # Base
    profile = Geom.naca(12, N=51)
    #profile = Geom.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
    bbox = Generator.bbox(profile)
    grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(2,2,1))
    D.Eq(grid1.P[1][0][0].x, baseLength)
    D.Eq(grid1.P[1][1][0].x, baseLength)
    D.Eq(grid1.P[0][1][0].y, 0.5*baseThickness)
    D.Eq(grid1.P[1][1][0].y, 0.5*baseThickness)
    D.Eq(grid1.P[0][0][0].y, -0.5*baseThickness)
    D.Eq(grid1.P[1][0][0].y, -0.5*baseThickness)
    spline1 = T1.Spline3('spline1', grid1, mesh=profile)
    sketch1 = T1.Sketch('sketch1', [spline1])

    # Top
    profile = Geom.naca(12, N=51)
    #naca = Geom.profile("SIKORSKY/SIKORSKYSSC-A07AIRFOIL")
    bbox = Generator.bbox(profile)
    grid2 = T1.Grid('grid2', bbox[0:3], bbox[3:], N=(2,2,1))
    D.Eq(grid2.P[1][0][0].x, topLength)
    D.Eq(grid2.P[1][1][0].x, topLength)
    D.Eq(grid2.P[0][1][0].y, 0.5*topThickness)
    D.Eq(grid2.P[1][1][0].y, 0.5*topThickness)
    D.Eq(grid2.P[0][0][0].y, -0.5*topThickness)
    D.Eq(grid2.P[1][0][0].y, -0.5*topThickness)
    spline2 = T1.Spline3('spline2', grid2, mesh=profile)
    sketch2 = T1.Sketch('sketch2', [spline2])
    D.Eq(sketch2.position.x, baseLength-topLength)
    D.Eq(sketch2.position.z, height)

    P1 = T1.Point('P1')
    P2 = T1.Point('P2')
    D.Eq(P2.z, height)
    D.Eq(P2.x, baseLength-topLength)

    P3 = T1.Point('P3')
    D.Eq(P3.x, baseLength)
    P4 = T1.Point('P4')
    D.Eq(P4.x, baseLength)
    D.Eq(P4.z, height)

    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P3, P4)

    sketch3 = T1.Sketch('sketch3', [line2])

    # body Surface
    surface1 = T1.Loft('surface1', listSketches=[sketch1, sketch2], listGuides=[])
    #surface1 = T1.MergeEdges('surface1', listSketches=[sketch1, sketch2, sketch3])

    # solve
    T1.solve()

    #T1.instantiate({'height': 1.0})

    return T1