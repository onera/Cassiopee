# 7A rotor blade with parametric twist
import Roms.Driver as D
import numpy

def createPart(name):

    T1 = D.Part(name)

    radius = 2.1         # blade raidus (m)
    solidity = 0.084     # rotor solidity
    nblades = 4          # number of blades
    chord = numpy.pi * radius * solidity / nblades

    angle = T1.Scalar('angle', 0.)
    angle.range = [0., 90., 1]

    # first sketch = OA213
    airfoil1 = T1.Airfoil('airfoil1', 'ONERA/ONERAOA213AIRFOIL')
    sketch1 = T1.Sketch('sketch1', [airfoil1], h=0.01)
    sketch1.set(rotCenter=(0,0,0), rotAxis=(0,1,0), rotAngle=90.)

    #sketch2 = sketch1.Copy('sketch1')
    #sketch2.shift(


    # solve
    T1.solve()

    return T1

if __name__ == "__main__":

    T1 = createPart('7ARotor1')
    # instantiate
    T1.instantiate({'angle': 0.})
    D.DRIVER.plot()
