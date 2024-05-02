# - refine (array) -
import Post as P
import Converter as C
import Generator as G

# Using indic (linear)
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
indic = C.array('indic', a[2].shape[1], 1, 1)
indic = C.initVars(indic, 'indic', 0)
C.setValue(indic, 50, [1])
C.setValue(indic, 49, [1])
a = P.refine(a, indic)
C.convertArrays2File(a, 'out.plt')

# Using butterfly
a = G.cartTetra((0,0,0), (2,1,1), (3,3,3))
a = P.exteriorFaces(a)
#a = C.initVars(a, "z = 0.1*{x}*{x}+0.2*{y}")
for i in range(6):
    a = P.refine(a, w=1./64.)
C.convertArrays2File(a, 'out.plt')
