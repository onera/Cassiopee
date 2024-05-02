# - cylinder (array) -
import Converter as C
import Transform as T
import Modeler.Models as Models

a = Models.cylinder(R1=1.,R2=1,N=20,h=10)

b = Models.cylinder(R1=0.1,R2=1,N=20,h=10)
b = T.translate(b, (0,4,0))

c = Models.cylinder(R1=0.1,R2=1,N=20,h=10,Rc=1.)
c = T.translate(c, (0,8,0))

C.convertArrays2File([a,b,c], 'out.plt')
