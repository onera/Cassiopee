# - star (array) -
import Converter as C
import Modeler.Models as Models
import Transform as T

a = Models.star(R1=1., R2=2., shift=0.5, N=10, h=0.1)
a = T.rotate(a, (0,0,0), (0,1,0), 90.)

b = Models.star(R1=1.5, R2=2., shift=0., N=20, h=0.2)
b = T.rotate(b, (0,0,0), (0,1,0), 90.)
b = T.translate(b, (0,5,0))

c = Models.star(R1=0.5, R2=2., shift=0.5, N=40, h=0.2)
c = T.rotate(c, (0,0,0), (0,1,0), 90.)
c = T.translate(c, (0,10,0))

C.convertArrays2File([a,b,c], 'out.plt')
