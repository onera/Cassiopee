# driver 9
import Roms.Driver as D
import Converter.Mpi as Cmpi

a = D.Scalar('a')
a.range = [0., 5., 0.5]
b = D.Scalar('b')
b.range = [0., 5., 1.]
c = D.Scalar('c')
c.range = [0., 5., 1.2]
d = D.Scalar('d')
d.range = [0., 5., 1.2]

eq1 = D.Eq(a, b)
eq2 = D.Eq(d, a+b)
ineq1 = D.Lt(b, c)

D.DRIVER.solve()

D.DRIVER.instantiate({'a':1., 'c':1.})

# replace here by walkDOE1 or walkDOE2 for parallel
pt = D.DRIVER.walkDOE()
while pt is not None:
    # <work here with pt>
    pt = D.DRIVER.walkDOE()
    print('Executing', Cmpi.rank, pt, flush=True)
