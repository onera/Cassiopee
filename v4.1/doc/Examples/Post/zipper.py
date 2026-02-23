# - zipper (array) -
import Converter as C
import Post as P
import Generator as G
import Transform as T

m1 = G.cylinder((0,0,0), 1, 5, 0., 360., 10., (50,50,1))
m1 = C.initVars(m1, 'cellN', 1.)

# Set cellN = 2 (interpolated points) to boundary
s = T.subzone(m1, (1,m1[3],1),(m1[2],m1[3],m1[4]))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,m1[3],1))
s = T.subzone(m1, (1,1,1),(m1[2],1,m1[4]))
s = C.initVars(s, 'cellN', 2)
m1 = T.patch(m1, s, (1,1,1))

ni = 30; nj = 40
m2 = G.cart((0,0,0), (10./(ni-1),10./(nj-1),-1), (ni,nj,1))
m2 = C.initVars(m2, 'cellN', 1.)

array = P.zipper([m1,m2],[])
C.convertArrays2File([array], 'out.plt')
