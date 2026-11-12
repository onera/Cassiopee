# Cap part from a sketch
import Roms.Driver as D
import numpy

def createPart(name):

    T1 = D.Part(name)

    # starting sketch
    thickness = T1.Scalar('thickness', 1.)
    thickness.range = [0.1, 10., 1]

    P1 = T1.Point('P1', (1, 0, 0))
    P2 = T1.Point('P2', (numpy.sqrt(2.)/2.,numpy.sqrt(2.)/2.,0))
    P3 = T1.Point('P3', (0, 1, 0))

    P4 = T1.Point('P4', (0, 1, 0))
    P5 = T1.Point('P5', (0, numpy.sqrt(2.)/2.,numpy.sqrt(2.)/2.))
    P6 = T1.Point('P6', (0, 0, 1))

    P7 = T1.Point('P7', (0, 0, 1))
    P8 = T1.Point('P8', (numpy.sqrt(2.)/2.,0, numpy.sqrt(2.)/2.))
    P9 = T1.Point('P9', (1, 0, 0))

    arc1 = T1.Arc('arc1', P1, P2, P3)
    arc2 = T1.Arc('arc2', P4, P5, P6)
    arc3 = T1.Arc('arc3', P7, P8, P9)

    sketch1 = T1.Sketch('sketch1', [arc1,arc2,arc3])

    surface1 = T1.Fill('surface1', sketch1, continuity=2)

    # solve
    T1.solve()

    return T1

if __name__ == "__main__":

    T1 = createPart('cap')
    # instantiate
    T1.instantiate({'thickness': 1.})
    D.DRIVER.plot()
