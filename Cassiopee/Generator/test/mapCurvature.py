# - mapCurvature (array) -
import Generator as G
import Converter as C
import Geom as D

ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,4.])
C.setValue(a, (1,3,1), [1.,3.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,2.])
b = D.bezier(a, density=10.)

b = G.mapCurvature(b, N=100, power=0.5, dir=1)
C.convertArrays2File([b], 'out.plt')
