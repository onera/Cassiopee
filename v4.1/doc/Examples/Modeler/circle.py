# - circle (array) -
import Converter as C
import Modeler.Models as Models
import Transform as T

a = Models.circle1(R=1., Rd=0.8, Nd=6, fracD=0.5, N=180)

b = Models.circle2(R=1., Rd=0.8, Nd=6, fracD=0.2, N=180)
b = T.translate(b, (2.5,0,0))

c = Models.circle2(R=1., Rd=1.2, Nd=9, fracD=0.2, N=180)
c = T.translate(c, (5,0,0))

d = Models.circle2(R=0.1, Rd=1., Nd=3, fracD=0.3, N=180)
d = T.translate(d, (0,2.5,0))

e = Models.circle2(R=1., Rd=0.3, Nd=3, fracD=0.7, N=180)
e = T.translate(e, (2.5,2.5,0))

C.convertArrays2File([a,b,c,d,e], 'out.plt')
