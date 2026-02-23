# - computeZb (pyTree) -
import Converter.PyTree as C
import Post.Rotor as PR
import math

Mtip = 0.6462; MU = 0.4
CHORD = 0.14; AR = 2.1
RoInf = 1.225; PInf = 101325.; AINF = 340.1
SIGMA = 4*CHORD / (math.pi*AR)
teff = C.convertFile2PyTree('stress_419.cgns')

zb = PR.computeZb(teff, 419., RoInf, AINF, Mtip, AR, SIGMA, relativeShaft=-12.12); print(zb)
#> [-1.101453126880057, -0.11259440192933268, 10.18004327358801]

accu = {}
zb = PR.computeZb(teff, 419., RoInf, AINF, Mtip, AR, SIGMA, relativeShaft=-12.12, accumulatorZb=accu)
ret = PR.exportAccumulatorPerRadius(accu, vars=['Xb','Yb','Zb'])
C.convertPyTree2File(ret, 'Zb.cgns')
