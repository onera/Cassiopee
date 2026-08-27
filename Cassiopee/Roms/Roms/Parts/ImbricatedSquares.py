# imbricated squares
import Roms.Driver as D
import sympy

def createSquare(T1, n, P1, L, alpha, lines):
    P2 = T1.Point(f'P{n}#2', (0,0,0))
    T1.Eq(P1.x, P1.x+L*sympy.cos(alpha))
    T1.Eq(P1.y, P1.y+L*sympy.sin(alpha))
    P3 = T1.Point(f'P{n}#3', (0,0,0))
    T1.Eq(P3.x, P2.x+L*sympy.sin(alpha))
    T1.Eq(P3.y, P2.y+L*sympy.cos(alpha))
    P4 = T1.Point(f'P{n}#4', (0,0,0))
    T1.Eq(P4.x, P3.x+L*sympy.cos(alpha))
    T1.Eq(P4.y, P3.y+L*sympy.sin(alpha))
    line1 = T1.Line('line1', P1, P2)
    line2 = T1.Line('line2', P2, P3)
    line3 = T1.Line('line3', P3, P4)
    line4 = T1.Line('line4', P4, P1)
    lines += [line1,line2,line3,line4]

def createPart(name):

    T1 = D.Part(name)

    # Create parameters
    width = T1.Scalar('width')
    width.range = [0., 2., 0.1]

    angle = T1.Scalar('angle')
    angle.range = [0., 30., 0.1]

    ndepth = T1.Scalar('depth')
    ndepth.range = [0., 30., 1.]

    lines = []
    n = 0
    for n in range(1):
        P1 = T1.Point(f'P{n}#1', (0,0,0))
        T1.Eq(P1.x, 0.)
        T1.Eq(P1.y, 0.)
        alpha = angle.v
        L = width.v
        createSquare(T1, n, P1, L, alpha, lines)

    sketch1 = T1.Sketch('sketch1',
                        lines,
                        h=[0.01,0.01,0.01])

    # solve result
    T1.solve()

    #T1.instantiate({'P0#2.y':1., 'P0#2.x':1, 'depth':1, 'angle':20., 'width':1})

    return T1