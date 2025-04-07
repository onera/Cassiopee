# - fittingPlaster (pyTree) -
import Generator.PyTree  as G
import Converter.PyTree  as C
import Geom.PyTree  as D

a = D.circle( (0,0,0), 1, N=50 )
a = C.convertArray2Tetra(a)
a = G.close(a)
b = G.fittingPlaster(a, bumpFactor=0.5)
C.convertPyTree2File(b, 'out.cgns')
