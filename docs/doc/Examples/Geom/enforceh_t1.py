# - enforceh (array) -
import Geom as D
import Transform as T
import Converter as C
import KCore.test as test

# broken line (STRUCT)
a = D.line((0,0,0), (1,0,0), N=10)
D.setF(a, 0, 1.); D.setF(a, -1, 0.5)
b = D.line((1,0,0), (2,1,0), N=50)
D.setF(b, 0, 0.5); D.setF(b, -1, 1.)
a = T.join([a,b])
a = D.enforceh(a, h=0.05)
test.testA([a], 1)
#C.convertArrays2File(a, 'out.plt')

# fourche (BAR)
a = D.line((0,0,0), (1,0,0), N=10)
D.setH(a, 0, 0.1); D.setH(a, -1, 0.01)
b = D.line((1,0,0), (2,1,0), N=50)
D.setH(b, 0, 0.01); D.setH(b, -1, 0.1)
c = D.line((1,0,0), (2,-1,0), N=100)
D.setH(c, 0, 0.01); D.setH(c, -1, 0.1)
A = [a,b,c]
A = C.convertArray2Hexa(A)
a = T.join(A)
a = D.enforceh(a, N=100)
#C.convertArrays2File(a, 'out.plt')
test.testA([a],2)

# Circle (STRUCT)
a = D.circle((0,0,0), 1, N=120)
D.setH(a, 0, 0.1)
D.setH(a, 60, 0.01)
D.setH(a, -1, 0.1)
a = D.enforceh(a, N=120)
#C.convertArrays2File(a, 'out.plt')
test.testA([a],3)
