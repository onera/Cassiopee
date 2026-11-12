# 7A rotor blade with parametric twist
import Roms.Driver as D
import numpy

# manque trailing edge, pb de croisement

def createPart(name):

    T1 = D.Part(name)

    radius = 2.1         # blade radius (m)
    solidity = 0.084     # rotor solidity
    nblades = 4          # number of blades
    chord = numpy.pi * radius * solidity / nblades

    angle = T1.Scalar('angle', 0.)
    angle.range = [-30., 0., 1]

    angle2 = T1.Scalar('angle2', 0.)
    T1.Eq(angle2, 0.8*angle)

    #angle3 = 0.5*angle
    angle3 = T1.Scalar('angle3', 0.)
    T1.Eq(angle3, 0.6*angle)

    # sketch = OA213, 0. incidence
    airfoil1 = T1.Airfoil('airfoil1', 'ONERA/ONERAOA213AIRFOIL')
    sketch1 = T1.Sketch('sketch1', [airfoil1], h=0.001)
    sketch1.set(scaleCenter=(0,0,0), scaleFactor=chord, position=(chord/4.,0,0.2*radius))
    sketch1.set(rotCenter=(chord/4.,0,0.2*radius), rotAxis=(0,0,1), rotAngle=0.)

    # OA213, incidence angle3
    sketch2b = sketch1.Copy('sketch2b')
    sketch2b.set(scaleCenter=(0,0,0), scaleFactor=chord, position=(chord/4.,0,0.5*radius))
    sketch2b.set(rotCenter=(chord/4.,0,0.5*radius), rotAxis=(0,0,1), rotAngle=angle3)

    # OA213, incidence angle
    sketch2 = sketch1.Copy('sketch2')
    sketch2.set(scaleCenter=(0,0,0), scaleFactor=chord, position=(chord/4.,0,0.75*radius))
    sketch2.set(rotCenter=(chord/4.,0,0.75*radius), rotAxis=(0,0,1), rotAngle=angle)

    # OA209, incidence angle2
    airfoil2 = T1.Airfoil('airfoil2', 'ONERA/ONERAOA209AIRFOIL')
    sketch3 = T1.Sketch('sketch3', [airfoil2], h=0.001)
    sketch3.set(scaleCenter=(0,0,0), scaleFactor=chord, position=(chord/4.,0,0.9*radius))
    sketch3.set(rotCenter=(chord/4.,0,0.9*radius), rotAxis=(0,0,1), rotAngle=angle2)

    # OA209, incidence angle
    sketch4 = sketch3.Copy('sketch4')
    sketch4.set(scaleCenter=(0,0,0), scaleFactor=chord, position=(chord/4.,0,1.*radius))
    sketch4.set(rotCenter=(chord/4.,0,1.*radius), rotAxis=(0,0,1), rotAngle=angle2)

    surface1 = T1.Loft('surface1', [sketch1,sketch2b,sketch2,sketch3,sketch4])

    # solve
    T1.solve()

    return T1

if __name__ == "__main__":

    T1 = createPart('7ARotor1')
    # instantiate
    T1.instantiate({'angle': 0.})
    D.DRIVER.plot()
