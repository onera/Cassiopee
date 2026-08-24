import Roms.Driver as D

# this is a part
def createPart(name):

    T1 = D.Part(name)

    # Create parameters
    radius = T1.Scalar('radius')
    radius.range = [0., 2., 0.1]

    # Create points
    P1 = T1.Point('P1', (0,0,0))
    T1.Eq(P1.x, 0.)
    T1.Eq(P1.y, 0.)

    # Create edges
    circle1 = T1.Circle('circle1', P1, radius)

    # Create sketch
    sketch1 = T1.Sketch('sketch1',
                        [circle1],
                        h=[0.01,0.01,0.01])

    # solve result
    #T1.solution = {P1.x: 0.0, P1.y: 0.0, P2.x: width, P3.x: width, P3.y: height, P4.y: height}
    #T1.params = [P4.y, P3.y, P3.x, P2.x, P1.y, P1.x, width, height]
    #T1.freeParams = [width, height]
    T1.solve()

    return T1
