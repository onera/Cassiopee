# - deformNormals (array) -
import Converter as C
import Generator as G
import Geom as D
import Transform as T

a = D.sphere((0,0,0), 1., 50)
a = C.convertArray2Hexa(a)
a = G.close(a)
b = C.initVars(a, '{alpha}=0.5*{x}')
b = C.extractVars(b, ['alpha'])
b = T.deformNormals(a, b, niter=2)
C.convertArrays2File([b], 'out.plt')
