# - computeThrustAndTorque (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import math
import KCore.test as test

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

ret = PR.computeThrustAndTorque(teff, 419., PInf, center=(0,0,0), relativeShaft=0.)
test.testO(ret[0], 1)
