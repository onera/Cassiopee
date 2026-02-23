# - refine (array) -
import Geom as D
import Transform as T
import Converter as C
import KCore.test as test

# broken line (STRUCT)
a = D.line((0,0,0), (1,0,0), N=10)
b = D.line((1,0,0), (2,1,0), N=30)
a = T.join([a,b])
a = D.refine(a, N=30)
test.testA([a], 1)
#C.convertArrays2File(a, 'out.plt')

# fourche (BAR)
a = D.line((0,0,0), (1,0,0), N=10)
b = D.line((1,0,0), (2,1,0), N=50)
c = D.line((1,0,0), (2,-1,0), N=100)
A = [a,b,c]
A = C.convertArray2Hexa(A)
a = T.join(A)
a = D.refine(a, N=50)
test.testA([a], 2)
#C.convertArrays2File(a, 'out.plt')

# Circle (STRUCT)
a = D.circle((0,0,0), 1, N=120)
a = D.refine(a, N=120)
test.testA([a], 3)
#C.convertArrays2File(a, 'out.plt')
