# - extractSlices (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import KCore.test as test
import math
test.TOLERANCE=1.e-8

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

accu1 = {}; accu2 = {}
psi = 419.; radii = [1.2,1.3,2.]
slices, CnM2, CmM2 = PR.extractSlices(teff, 'Blade7A_00', psi, radii,
                                      RoInf, PInf, AINF, Mtip, AR, CHORD, MU,
                                      accumulatorCnM2=accu1, accumulatorCmM2=accu2,
                                      localFrame=True,
                                      relativeShaft=-12.12)
test.testT(slices, 1)

exp = PR.exportAccumulatorMap(accu1, vars=['CnM2x','CnM2y','CnM2z'])
test.testT(exp, 2)
exp = PR.exportAccumulatorMap(accu1, vars=['CmM2x','CmM2y','CmM2z'])
test.testT(exp, 3)
