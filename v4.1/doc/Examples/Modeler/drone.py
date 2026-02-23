# - drone (array) -
import Converter as C
import Geom as D
import Modeler.Models as Models
import Transform as T

# bras de drone
dx = 0.01; dy = 0.02; dz = 0.1
P0=(-dx/2,-dy/2,0); P1=(dx/2,dy/2,0)

section = Models.box2D(P0, P1, r=0.3, fill=False)
PA = (0,0,0); PB = (0,0,dz)
line = D.line(PA, PB)
arm = D.lineDrive(section, line)

cyl = Models.cylinder(R1=dx, R2=dx, N=120, h=dy)
cyl = T.rotate(cyl, (0,0,0), (1,0,0), 90.)
cyl = T.translate(cyl, (0,dy/2,0))
C.convertArrays2File([arm, cyl], 'out.plt')


# body de drone
dx = 0.1; dy = 0.2; dz = 0.05
P0=(-dx/2,-dy/2,-dz/2); P1=(dx/2,dy/2,dz/2)
body = Models.box(P0,P1,chamfer=0.01)
#C.convertArrays2File([body], 'out.plt')
