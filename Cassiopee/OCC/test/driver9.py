# driver 9
import OCC.Driver as D

a = D.Scalar('a', 1.)
a.range = [0,5]
b = D.Scalar('b', 1.)
b.range = [0,5]
c = D.Scalar('c', 1.)
c.range = [0,5]

eq1 = D.Eq(a.s, b.s)
ineq1 = D.Lt(b.s, c.s)

D.DRIVER.solve2()

D.DRIVER.instantiate({'a':1, 'c':1})