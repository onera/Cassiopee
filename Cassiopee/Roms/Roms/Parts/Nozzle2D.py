# parametric 2D nozzle
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    #==============
    # parameters
    #==============

    # throat height
    Ht = T1.Scalar('Ht', 0.5)
    Ht.range = [0.1, 2., 0.1]

    # inlet height
    Hi = T1.Scalar('Hi', 1.)
    Hi.range = [0.1, 2., 0.1]

    # exit height
    He = T1.Scalar('He', 1.)
    He.range = [0.1, 2., 0.1]

    # convergent length
    Lc = T1.Scalar('Lc', 2.)
    Lc.range = [0.1,3]

    # divergent length
    Ld = T1.Scalar('Ld', 2.)
    Ld.range = [0.1,3]

    #=======================
    # create points
    #=======================

    # Create points P1 (throat)
    P1 = T1.Point('P1')
    P1.x.range = [-1,1]
    T1.Eq(P1.y, Ht)

    # Create points P2 (inlet)
    P2 = T1.Point('P2')
    T1.Eq(P2.x, -Lc)
    T1.Eq(P2.y, Hi)

    # Create points P3 (exit)
    P3 = T1.Point('P3')
    T1.Eq(P3.x, Ld)
    T1.Eq(P3.y, He)

    #=================
    # create entities
    #=================
    # basic with lines
    line1 = T1.Line('line1', P2, P1)
    line2 = T1.Line('line2', P1, P3)

    # create symmetric part
    P1S = T1.Point('P1S')
    T1.Eq(P1S.x, P1.x)
    T1.Eq(P1S.y, -P1.y)
    P2S = T1.Point('P2S')
    T1.Eq(P2S.x, P2.x)
    T1.Eq(P2S.y, -P2.y)
    P3S = T1.Point('P3S')
    T1.Eq(P3S.x, P3.x)
    T1.Eq(P3S.y, -P3.y)
    line3 = T1.Line('line3', P2S, P1S)
    line4 = T1.Line('line4', P1S, P3S)

    line5 = T1.Line('line5', P2, P2S)
    line6 = T1.Line('line6', P3, P3S)

    sketch1 = T1.Sketch('sketch1', [line1,line2,line5,line4,line3,line6])

    #=======================
    # solve and instantiate
    #=======================
    T1.solve()

    #T1.instantiate({'P1.x':0.,
    #                'Hi': 1.,
    #                'Ht': 0.2,
    #                'He': 1.,
    #                'Lc': 2.,
    #                'Ld': 2.})

    return T1