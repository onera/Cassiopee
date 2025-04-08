# - smooth (pyTree) -
import Transform.PyTree as T
import Geom.PyTree as D
import Converter.PyTree as C

a = D.sphere6((0,0,0), 1, N=20)
b = T.smooth(a, eps=0.5, niter=20)
C.convertPyTree2File(b, "out.cgns")
