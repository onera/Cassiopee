# Circle part
import Roms.Driver as D

def createPart(name):

    T1 = D.Part(name)

    # Create parameters
    radius = T1.Scalar('radius')
    radius.range = [0., 2., 0.1]

    # Create points
    P1 = T1.Point('P1', (0,0,0))

    # Create edges
    circle1 = T1.Circle('circle1', P1, radius)

    # Create sketch
    sketch1 = T1.Sketch('sketch1',
                        [circle1],
                        h=[0.01,0.01,0.01])

    # solve result
    T1.solve()

    return T1

if __name__ == "__main__":

    T1 = createPart('circle1')
    # instantiate
    T1.instantiate({'radius':1.})
    D.DRIVER.plot()
