# parametric 3D nozzle
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
    Lc.range = [0.1, 3.]

    # divergent length
    Ld = T1.Scalar('Ld', 2.)
    Ld.range = [0.1, 3.]

    # convergent curvature
    Rc = T1.Scalar('Rc', 1.)

    # divergent curvature
    Rd = T1.Scalar('Rd', 1.)

    #=======================
    # create points
    #=======================

    # Create points P1 (throat)
    P1 = T1.Point('P1')
    P1.x.range = [-2., 2.]
    T1.Eq(P1.y, Ht)

    # Create points P2 (inlet)
    P2 = T1.Point('P2')
    T1.Eq(P2.x, -Lc)
    T1.Eq(P2.y, Hi)

    P2a = T1.Point('P2a')
    T1.Eq(P2a.x, -0.1*Lc)
    T1.Eq(P2a.y, Ht)

    P2b = T1.Point('P2b')
    T1.Eq(P2b.x, -0.5*Lc)
    T1.Eq(P2b.y, Ht+0.5*(Ht+Hi))

    # Create points P3 (exit)
    P3 = T1.Point('P3')
    T1.Eq(P3.x, Ld)
    T1.Eq(P3.y, He)

    P3a = T1.Point('P3a')
    T1.Eq(P3a.x, 0.1*Ld)
    T1.Eq(P3a.y, Ht)

    P3b = T1.Point('P3b')
    T1.Eq(P3b.x, 0.5*Ld)
    T1.Eq(P3b.y, Ht+0.5*(Ht+He))

    #=================
    # create entities
    #=================
    # basic with lines
    spline1 = T1.Spline1('spline1', [P2, P2b, P2a, P1])
    spline2 = T1.Spline1('spline2', [P1, P3a, P3b, P3])

    sketch1 = T1.Sketch('sketch1', [spline1,spline2])
    surface1 = T1.Revolve('surface1', sketch1, center=(0,0,0), axis=(1,0,0), angle=360., h=(0.1,0.1,0.01))

    #=======================
    # solve and instantiate
    #=======================
    T1.solve()

    #T1.instantiate({'P1.x':0.,
    #                 'Hi': 1.,
    #                 'Ht': 0.2,
    #                 'He': 1.,
    #                 'Lc': 2.,
    #                 'Ld': 2.})

    return T1