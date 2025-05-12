# - extractSlices (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import math

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

accu = {}
psi = 419.; radius = [1.2,1.3,2.]
slices, CnM2, CmM2 = PR.extractSlices(teff, 'Blade7A_00', psi, radius,
                                      RoInf, PInf, AINF, Mtip, AR, CHORD, MU,
                                      accumulatorCnM2=accu, localFrame=True,
                                      relativeShaft=-12.12)

# export for Cp, Cf
C.convertPyTree2File(slices, 'slice.cgns')

# export for CnM2 maps
exp = PR.exportAccumulatorMap(accu, vars=['CnM2x','CnM2y','CnM2z'])
C.convertPyTree2File(exp, 'map.cgns')
