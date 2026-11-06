# - polyC1Mesher (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

tb = C.convertFile2PyTree('curve1.svg', nptsCurve=100, nptsLine=400)
z = T.homothety(tb[2][1][2][0], (0.,0.,0.), 0.01)
z = T.reorder(z, (-1,2,3))
h = 0.1; hf = 0.001; density = 100; splitCrit = 10.
res = G.polyC1Mesher(z, h, hf, density, splitCrit)
zones = res[0]; h = res[1]; density = res[2]
t = C.newPyTree(['Base']); t[2][1][2] += zones
C.convertPyTree2File(t, 'out.cgns')
