# Cap part from a sketch
import Roms.Driver as D
import numpy

def createPart(name):

    T1 = D.Part(name)

    # starting sketch
    v1 = T1.Scalar('v1', 1.)
    v1.range = [0.1, 10., 1]

    P1 = T1.Point('P1', (1, 0, 0))
    T1.Eq(P1.x, v1)
    P2 = T1.Point('P2', (numpy.sqrt(2.)/2.,numpy.sqrt(2.)/2.,0))
    P3 = T1.Point('P3', (0, 1, 0))

    #P4 = T1.Point('P4', (0, 1, 0))
    P5 = T1.Point('P5', (0, numpy.sqrt(2.)/2.,numpy.sqrt(2.)/2.))
    P6 = T1.Point('P6', (0, 0, 1))

    #P7 = T1.Point('P7', (0, 0, 1))
    P8 = T1.Point('P8', (numpy.sqrt(2.)/2.,0, numpy.sqrt(2.)/2.))
    #P9 = T1.Point('P9', (1, 0, 0))

    spline1 = T1.Spline2('spline1', [P1, P2, P3])
    spline2 = T1.Spline2('spline2', [P3, P5, P6])
    spline3 = T1.Spline2('spline3', [P6, P8, P1])
    sketch1 = T1.Sketch('sketch1', [spline1,spline2,spline3])

    spline4 = spline1.Symm('spline4', plane='xz', center=(0,0,0))
    spline5 = spline2.Symm('spline5', plane='xz', center=(0,0,0))
    spline6 = spline3.Symm('spline6', plane='xz', center=(0,0,0))
    sketch2 = T1.Sketch('sketch2', [spline4,spline5,spline6])

    spline7 = spline1.Symm('spline7', plane='xy', center=(0,0,0))
    spline8 = spline2.Symm('spline8', plane='xy', center=(0,0,0))
    spline9 = spline3.Symm('spline9', plane='xy', center=(0,0,0))
    sketch3 = T1.Sketch('sketch3', [spline7,spline8,spline9])

    surface1 = T1.Fill('surface1', sketch1, continuity=1)
    surface2 = T1.Fill('surface2', sketch2, continuity=1)
    surface3 = T1.Fill('surface3', sketch3, continuity=1)

    # solve
    T1.solve()

    return T1

if __name__ == "__main__":

    T1 = createPart('cap')
    # instantiate
    T1.instantiate({'v1': 1.})
    D.DRIVER.plot()
