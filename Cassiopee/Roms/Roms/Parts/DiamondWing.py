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
    thickness1.range = [0.02, 0.2, 0.01]
    angle2 = T1.Scalar('angle2', 1.)
    angle2.range = [0., 5., 1.]

    spanLength = T1.Scalar('spanLength', 10.)
    tipShift = T1.Scalar('tipShift', 0.)

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

    line1 = T1.Line('line1', P1, P3)
    line2 = T1.Line('line2', P3, P2)
    line3 = T1.Line('line3', P2, P4)
    line4 = T1.Line('line4', P4, P1)
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
    
    line5 = T1.Line('line5', P5, P7)
    line6 = T1.Line('line6', P7, P6)
    line7 = T1.Line('line7', P6, P8)
    line8 = T1.Line('line8', P8, P5)
    sketch2 = T1.Sketch('sketch2', [line5,line6,line7,line8])
    T1.Eq(sketch2.rotAngle, angle2)
    T1.Eq(sketch2.position.x, spanLength)
    T1.Eq(sketch2.position.y, tipShift)

    surface = T1.Loft('surface', listSketches=[sketch1, sketch2], close=False)

    #=======================
    # solve and instantiate
    #======================

    # solve
    T1.solve()

    # instantiate
    #T1.instantiate({'tipPosx': 0., 'tipPosz': 2.})

    return T1
