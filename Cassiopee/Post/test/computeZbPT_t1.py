# - computeZb (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import KCore.test as test
import math

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

zb = PR.computeZb(teff, 419., RoInf, AINF, Mtip, AR, SIGMA, relativeShaft=-12.12)
test.testO(zb, 1)
