# driver 9
import Roms.Driver as D
import Converter.Mpi as Cmpi

T1 = D.Part("Part1")

# parameters
a = T1.Scalar('a')
a.range = [0., 5., 0.5]
b = T1.Scalar('b')
b.range = [0., 5., 1.]
c = T1.Scalar('c')
c.range = [0., 5., 1.2]
d = T1.Scalar('d')
d.range = [0., 5., 1.2]

eq1 = T1.Eq(a, b)
eq2 = T1.Eq(d, a+b)
c1 = T1.Lt(b, c)
c2 = T1.Ne(a, 0.)

T1.solve()
T1.instantiate({'a':1., 'c':1.})

# replace here by walkDOE1 or walkDOE2 for parallel
point = T1.walkDOE()
while point is not None:
    # <work here with pt>
    point = T1.walkDOE()
    print('Executing', Cmpi.rank, point, flush=True)
