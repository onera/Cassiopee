import Models
import Converter.PyTree as C

def box(Pmin, Pmax, chamfer=-1.):
    a = Models.box(Pmin, Pmax, chamfer)
    return C.convertArrays2ZoneNode('box', [a])

def box2D(Pmin, Pmax, r=0.):
    a = Models.box2D(Pmin, Pmax, r)
    return C.convertArrays2ZoneNode('box', [a])

def ellipse2D(Pmin, Pmax):
    a = Models.ellipse2D(Pmin, Pmax)
    return C.convertArrays2ZoneNode('ellipse', [a])

def star(R1=1., R2=2., shift=0.5, N=10, h=1.):
    a = Models.star(R1, R2, shift, N, h)
    return C.convertArrays2ZoneNode('star', [a])

def circle1(R=1., Rd=0.8, Nd=10, fracD=0.5, N=180):
    a = Models.circle1(R, Rd, Nd, fracD, N)
    return C.convertArrays2ZoneNode('circle', [a])

def circle2(R=1., Rd=0.8, Nd=10, fracD=0.5, N=180):
    a = Models.circle2(R, Rd, Nd, fracD, N)
    return C.convertArrays2ZoneNode('circle', [a])

def column(R=0.2, N=10, h=1.):
    a = Models.column(R, N, h)
    return C.convertArrays2ZoneNode('column', [a])

def column2(R1=0.2, R2=0.2, N=10, h=1.):
    a = Models.column2(R1, R2, N, h)
    return C.convertArrays2ZoneNode('column', [a])

def column3(R1=0.2, R2=0.2, N=10, h=1.):
    a = Models.column3(R1, R2, N, h)
    return C.convertArrays2ZoneNode('column', [a])

def cylinder(R1=1.,R2=1.,N=10,h=1.,Rc=-1):
    a = Models.cylinder(R1, R2, N, h, Rc)
    return C.convertArrays2ZoneNode('cylinder', [a])
