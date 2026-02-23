# - polyLineMesher (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T

tb = C.convertFile2PyTree('fusee.plt')
tb = G.close(tb,1e-2); tb = T.reorder(tb,(-1,2,3))
h = 0.02; hf = 0.0001; density = 500

res = G.polyLineMesher(tb[2][1][2][0], h, hf, density)
zones = res[0]; h = res[1]; density = res[2]
t = C.newPyTree(['PolyC1']); t[2][1][2] += zones
C.convertPyTree2File(t, 'out.cgns')
