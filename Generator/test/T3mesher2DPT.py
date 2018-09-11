# - T3mesher2D (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D

a = D.circle((0,0,0), 1, N=50)
a = C.convertArray2Tetra(a);
a = G.close(a)

b = G.T3mesher2D(a, triangulateOnly=0, grading=1.2, metric_interp_type=0) # linear metric interpolation
t = C.newPyTree(['Base',2,b])
C.convertPyTree2File(t, 'outL.cgns')

b = G.T3mesher2D(a, triangulateOnly=0, grading=1.2, metric_interp_type=1) # geometric metric interpolation
t = C.newPyTree(['Base',2,b])
C.convertPyTree2File(t, 'outG.cgns')
