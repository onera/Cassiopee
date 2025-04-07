# - reorderAll (array) -
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

def f(x,y):
    return 3*x*y

ni = 30; nj = 40
m1 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,1))
#-------------------------------------------------
# test 1 : blocs recouvrants pareillement orientes
#-------------------------------------------------
m2 = T.rotate(m1, (0.2,0.2,0.), (0.,0.,1.), 15.)
A = [m1,m2]
A = C.initVars(A,'F',f,['x','y'])
B = T.reorderAll(A,1)
test.testA(B,1)

#-------------------------------------------------
# test 2 : blocs recouvrants orientes differemment
#-------------------------------------------------
m2 = T.reorder(m2,(-1,2,3))
A = [m1,m2]
A = C.initVars(A,'F',f,['x','y'])
B = T.reorderAll(A,1)
test.testA(B,2)

#---------------------------------
# test 3: blocs sans intersection
#---------------------------------
m2 = T.translate(m1, (0.,12,0.))
A = [m1,m2]
A = C.initVars(A,'F',f,['x','y'])
B = T.reorderAll(A)
test.testA(B,3)
