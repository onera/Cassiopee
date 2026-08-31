# parametric missile
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    #===================
    # Create parameters
    #===================
    radius1 = T1.Scalar('radius1', 10.)
    radius1.range = [5., 15., 1.]

    radius2 = T1.Scalar('radius2', 5.)
    radius2.range = [0.1, 8., 1.]

    # length back
    length1 = T1.Scalar('length1', 100.)
    T1.Eq(length1, 100.)

    # length middle
    length2 = T1.Scalar('length2', 20.)
    T1.Eq(length2, 20.)

    # length forward
    length3 = T1.Scalar('length3', 10.)
    T1.Eq(length3, 10.)

    #=================
    # Create points
    #=================
    P1 = T1.Point('P1')
    T1.Eq(P1.y, radius1)

    P2 = T1.Point('P2')
    T1.Eq(P2.y, radius1)
    T1.Eq(P2.x, -length1)

    P3 = T1.Point('P3')
    T1.Eq(P3.y, radius2)
    T1.Eq(P3.x, -length1-length2)

    P4 = T1.Point('P4')
    T1.Eq(P4.x, -length1-length2-length3)
    T1.Eq(P4.y, P3.y/20.)

    P5 = T1.Point('P5')
    T1.Eq(P5.x, -length1-length2-length3)

    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    spline1 = T1.Spline1('spline1', [P3,P4,P5])

    # sketch1
    sketch1 = T1.Sketch('sketch1', [line1, line2, spline1])

    # surface1 - body
    surface1 = T1.Revolve('surface1', sketch1, center=(0,0,0), axis=(1,0,0), angle=360.)

    # body closure
    circle1 = T1.Circle('circle1', (0,0,0), radius1)
    sketch2 = T1.Sketch('sketch2', [circle1])
    surface2 = T1.Fill('surface2', sketch2)

    surface3 = T1.Merge('surface3', [surface1, surface2])

    #=======================
    # solve and instantiate
    #=======================
    T1.solve()

    #T1.instantiate({'radius1': 10., 'radius2': 4.})

    return T1