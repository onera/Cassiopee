# Square part
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    # Create parameters
    height = T1.Scalar('height')
    height.range = [0., 2., 0.1]

    width = T1.Scalar('width')
    width.range = [0., 2., 0.1]

    # Create points
    P1 = T1.Point('P1', (0,0,0))
    T1.Eq(P1.x, 0.)
    T1.Eq(P1.y, 0.)

    P2 = T1.Point('P2', (1,0,0))
    T1.Eq(P2.x, P1.x + width)

    P3 = T1.Point('P3', (1,1,0))
    T1.Eq(P3.x, P1.x + width)

    P4 = T1.Point('P4', (0,1,0))
    T1.Eq(P3.y, P4.y)
    T1.Eq(P3.y, P1.y + height)

    # Create lines
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)

    # Create sketch
    sketch1 = T1.Sketch('sketch1',
                        [line1, line2, line3, line4],
                        h=[0.01,0.01,0.01])

    # solve result
    T1.solve()

    # example of instantiation
    #T1.instantiate({'radius':1.})

    return T1