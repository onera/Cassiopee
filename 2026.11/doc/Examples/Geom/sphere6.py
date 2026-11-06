# - sphere6 (array) -
import Geom as D
import Converter as C

a = D.sphere6((0,0,0), 1., N=20)
b = D.sphere6((3,3,0), 1.2, N=20, ntype='QUAD')
C.convertArrays2File(a+[b], "out.plt")
