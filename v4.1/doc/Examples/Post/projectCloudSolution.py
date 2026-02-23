# -projectCloudSolution (array)
import Converter as C
import Geom as D
import Post as P
import Transform as T
import Generator as G

a = D.sphere((0,0,0),1.,N=20)
a = C.convertArray2Tetra(a); a = G.close(a)
b = D.sphere6((0,0,0),1.,N=15)
b = C.convertArray2Tetra(b); b = T.join(b)
pts = C.convertArray2Node(b)
pts = C.initVars(pts,'{F}={x}*{y}')
a = C.initVars(a,'F', 0.)
a = P.projectCloudSolution(pts,a, loc='nodes')
C.convertArrays2File(a,"out.plt")
