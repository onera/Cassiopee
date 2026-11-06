# - fittingPlaster (pyTree) -
import Generator.PyTree  as G
import Converter.PyTree  as C
import Geom.PyTree  as D
import KCore.test as test

a = D.circle( (0,0,0), 1, N=50 )
a = C.convertArray2Tetra(a)
a = G.close(a)
b = G.fittingPlaster(a, bumpFactor=0.5)
t = C.newPyTree(['Base',2]); t[2][1][2].append(b)
test.testT(t)
