# dimaon parametric wing
import Roms.Driver as D
import Geom
import Generator
import Transform

def createPart(name):

    T1 = D.Part(name)

    #======================
    # Create parameters
    #======================
    chord1 = T1.Scalar('chord1', 1.)
    chord1.range = [0.1, 1.5, 0.1]
    thickness1 = T1.Scalar('thickness1', 0.1)
    thickness1.range = [0.02, 0.2, 0.01]
    angle1 = T1.Scalar('angle1', 1.)
    angle1.range = [0., 5., 1.]

    chord2 = T1.Scalar('chord2', 1.)
    chord2.range = [0.1, 1.5, 0.1]
    thickness2 = T1.Scalar('thickness2', 0.1)
    thickness2.range = [0.02, 0.2, 0.01]
    angle2 = T1.Scalar('angle2', 1.)
    angle2.range = [0., 5., 1.]

    spanLength = T1.Scalar('spanLength', 10.)
    spanLength.range = [0.1, 5., 0.1]
    tipShift = T1.Scalar('tipShift', 0.)
    tipShift.range = [0.,1.,0.1]

    #=================
    # create entities
    #=================

    # Sketch1
    P1 = T1.Point('P1', (0,0,0))
    P2 = T1.Point('P2', (0,0,0))
    T1.Eq(P2.y, chord1)
    P3 = T1.Point('P3', (0,0,0))
    T1.Eq(P3.y, 0.5*chord1)
    T1.Eq(P3.z, 0.5*thickness1)
    P4 = T1.Point('P4', (0,0,0))
    T1.Eq(P4.y, 0.5*chord1)
    T1.Eq(P4.z, -0.5*thickness1)

    line1 = T1.Line('line1', P1, P4)
    line2 = T1.Line('line2', P4, P2)
    line3 = T1.Line('line3', P2, P3)
    line4 = T1.Line('line4', P3, P1)
    sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4])
    T1.Eq(sketch1.rotAngle, angle1)

    # Sketch2
    P5 = T1.Point('P5', (0,0,0))
    P6 = T1.Point('P6', (0,0,0))
    T1.Eq(P6.y, chord2)

    P7 = T1.Point('P7', (0,0,0))
    T1.Eq(P7.y, 0.5*chord2)
    T1.Eq(P7.z, 0.5*thickness2)

    P8 = T1.Point('P8', (0,0,0))
    T1.Eq(P8.y, 0.5*chord2)
    T1.Eq(P8.z, -0.5*thickness2)

    line5 = T1.Line('line5', P5, P8)
    line6 = T1.Line('line6', P8, P6)
    line7 = T1.Line('line7', P6, P7)
    line8 = T1.Line('line8', P7, P5)
    sketch2 = T1.Sketch('sketch2', [line5,line6,line7,line8])
    T1.Eq(sketch2.rotAngle, angle2)
    T1.Eq(sketch2.position.x, spanLength)
    T1.Eq(sketch2.position.y, tipShift)

    surface1 = T1.Loft('surface1', listSketches=[sketch1, sketch2], close=True)

    # symetrisation
    symetric = True
    if symetric:
        
        # Sketch2S
        P5S = T1.Point('P5S', (0,0,0))
        P6S = T1.Point('P6S', (0,0,0))
        T1.Eq(P6S.y, chord2)

        P7S = T1.Point('P7S', (0,0,0))
        T1.Eq(P7S.y, 0.5*chord2)
        T1.Eq(P7S.z, 0.5*thickness2)

        P8S = T1.Point('P8S', (0,0,0))
        T1.Eq(P8S.y, 0.5*chord2)
        T1.Eq(P8S.z, -0.5*thickness2)

        line5S = T1.Line('line5S', P5S, P8S)
        line6S = T1.Line('line6S', P8S, P6S)
        line7S = T1.Line('line7S', P6S, P7S)
        line8S = T1.Line('line8S', P7S, P5S)
        sketch2S = T1.Sketch('sketch2S', [line5S,line6S,line7S,line8S])
        T1.Eq(sketch2S.rotAngle, angle2)
        T1.Eq(sketch2S.position.x, -spanLength)
        T1.Eq(sketch2S.position.y, tipShift)

        surface2 = T1.Loft('surface2', listSketches=[sketch1, sketch2S], close=True)

        surface = T1.Merge('surface', [surface1, surface2])
        #surface = T1.Union('surface', [surface1], [surface2])

    #=======================
    # solve and instantiate
    #======================

    # solve
    T1.solve()

    # instantiate
    #T1.instantiate({'tipPosx': 0., 'tipPosz': 2.})

    return T1
