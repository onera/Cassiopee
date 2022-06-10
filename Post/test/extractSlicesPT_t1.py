# - extractSlices (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import KCore.test as test
import math

MTIP = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
ROINF = 1.225; PINF = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

accu1 = {}; accu2 = {} 
psi = 419.; radius = [1.2,1.3,2.]
slices = PR.extractSlices(teff, 'Blade7A_00', psi, radius,
                          ROINF, PINF, AINF, MTIP, AR, CHORD, MU,
                          accumulatorCnM2=accu1, accumulatorCmM2=accu2, 
                          localFrame=True,
                          relativeShaft=-12.12)
test.testT(slices, 1)

exp = PR.exportAccumulatorMap(accu1, vars=['CnM2x','CnM2y','CnM2z'])
test.testT(exp, 2)
exp = PR.exportAccumulatorMap(accu1, vars=['CmM2x','CmM2y','CmM2z'])
test.testT(exp, 3)
