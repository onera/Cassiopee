# Wedge part
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    #==============
    # parameters
    #==============

    # forward length
    Lf = T1.Scalar('Lf', 2.)
    Lf.range = [0.1, 3.]

    # backward length
    Lb = T1.Scalar('Lb', 0.1)
    Lb.range = [0.1, 1.]

    # top height
    Ht = T1.Scalar('Ht', 0.1)
    Ht.range = [0.1, 1.]

    # bottom height
    Hb = T1.Scalar('Hb', 0.1)
    Hb.range = [0.1, 1.]

    # side length
    Ls = T1.Scalar('Ls', 1.)
    Ls.range = [0.1, 2.]

    #=======================
    # create points
    #=======================

    # forward point
    P0 = T1.Point('P0')
    T1.Eq(P0.x, -Lf)

    P1 = T1.Point('P1')
    T1.Eq(P1.y, -Ls)

    P2 = T1.Point('P2')
    T1.Eq(P2.y, +Ls)

    P3 = T1.Point('P3')
    T1.Eq(P3.z, Ht)

    P4 = T1.Point('P4')
    T1.Eq(P4.z, -Hb)

    P5 = T1.Point('P5')
    T1.Eq(P5.x, +Lb)

    #=================
    # create entities
    #=================

    surface1 = T1.FillLinear('surface1', [P0,P3,P1])
    surface2 = T1.FillLinear('surface2', [P0,P3,P2])
    surface3 = T1.FillLinear('surface3', [P0,P4,P1])
    surface4 = T1.FillLinear('surface4', [P0,P4,P2])

    surface5 = T1.FillLinear('surface5', [P5,P3,P1])
    surface6 = T1.FillLinear('surface6', [P5,P3,P2])
    surface7 = T1.FillLinear('surface7', [P5,P4,P1])
    surface8 = T1.FillLinear('surface8', [P5,P4,P2])

    surface = T1.Merge('surface', [surface1,surface2,surface3,surface4,surface5,surface6,surface7, surface8], h=(0.1,0.1,0.01))

    #=======================
    # solve and instantiate
    #=======================
    T1.solve()

    # example of instance
    #T1.instantiate({'Lf': 2,
    #                'Lb': 0.1,
    #                'Ht': 0.1,
    #                'Hb': 0.1,
    #                'Ls': 1.})

    return T1
