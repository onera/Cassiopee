# - expression -
import Generator as G
import Converter as C
import Converter.expression as expr
import numpy as np
from math import *

def test_function(formule, func) :
	print("testing formule ", formule)
	a = expr.ast(formule)
	N = 30
	coords = G.cart((0.01, 0.01, 0.01), (1, 1, 1), (N, 1, 1), api=1)
	C._addVars(coords, 'result')
	a.run(coords)
	x = coords[1][0,:]
	y = coords[1][1,:]
	z = coords[1][2,:]
	r = coords[1][3,:]
	for i in range(N):
		v = func(x[i],y[i],z[i])
		if abs(r[i] - v ) > 1.E-14 * abs(v):
			print("Erreur calcul : resultat trouve avec expression : ", r[i], " et avec formule : ", v)


a = expr.ast("{x}**2+{y}**2+{z}**2")

xs = 3.
ys = 2.
zs = -1.
print("{} avec x = {}, y = {} et z = {} => {}".format(a, xs, ys, zs,
                                                      a(x=xs, y=ys, z=zs)))

coords = C.array("x,y,z", 3, 1, 1)
coords[1][:, :] = np.array([[1., 2., 3.], [2., 3., 4.], [-1., -3., -5.]])

xv = np.array([1.0, 2.0, 3.0])
yv = np.array([2., 3., 4.])
zv = np.array([-1., -3., -5.])
print("{} avec x = {}, y = {} et z = {} => {}".format(a, xv, yv, zv,
                                                      a(x=xv, y=yv, z=zv)))

t = a(coords)
print("{} avec {} => {}".format(a, coords, a(coords)))

a.eval(coords, "x", coords)
print(coords)
a.eval(coords, "crdL2", coords)
print(coords)
res = C.array("norm", 3, 1, 1)
a.eval(res, "norm", coords)
print(res)

a = expr.ast("{norm} = {x}**2+{y}**2+{z}**2")
coords = C.array("x,y,z,norm", 3, 1, 1)
coords[1][:, :] = np.array([[1., 2., 3.], [2., 3., 4.], [-1., -3., -5.],
                            [0., 0., 0.]])
a.run(coords)
print(coords)

crds = G.cart((0, 0, 0), (1, 1, 1), (10, 10, 1), api=1)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cart((0, 0, 0), (1, 1, 1), (10, 10, 1), api=2)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cartHexa((0, 0, 0), (1, 1, 1), (10, 10, 1), api=1)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cartHexa((0, 0, 0), (1, 1, 1), (10, 10, 1), api=2)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cartNGon((0, 0, 0), (1, 1, 1), (10, 10, 1), api=1)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cartNGon((0, 0, 0), (1, 1, 1), (10, 10, 1), api=2)
C._addVars(crds, 'norm')
a.run(crds)
print(crds)

crds = G.cart((0, 0, 0), (1, 1, 1), (3, 3, 1), api=1)
C._addVars(crds, 'norm')
da = expr.derivate(a)
C._addVars(crds, 'd_norm')
shp = crds[1][0].shape
da.run(crds, d_x=np.ones(shp), d_y=np.zeros(shp), d_z=np.zeros(shp))
print(da)
print(crds)

test_function("{result} = 3 * {x} + sin({y})", lambda x,y,z : 3 * x + sin(y))
test_function("{result} = cos(pi*{x})", lambda x,y,z : cos(pi*x))
test_function("{result} = sin(pi*{x})", lambda x,y,z : sin(pi*x))
test_function("{result} = exp(pi*{x})", lambda x,y,z : exp(pi*x))
test_function("{result} = log(pi*{x})", lambda x,y,z : log(pi*x))
test_function("{result} = cosh(pi*{x})", lambda x,y,z : cosh(pi*x))
test_function("{result} = sinh(pi*{x})", lambda x,y,z : sinh(pi*x))
test_function("{result} = abs(pi-{x})", lambda x,y,z : abs(pi-x))

test_function("{result} = cos(2*pi*{x})*sin(0.03+1.97*{y})", lambda x,y,z : cos(2*pi*x)*sin(0.03+1.97*y))
#test_function("{result} = (cos(0.03+1.97*{x})+log(tan(1/2*(0.03+1.97*{x})+0.00000001)))", lambda x,y,z : cos(0.03+1.97*x)+log(tan(1/2*(0.03+1.97*x)+0.00000001)))
test_function("{result} = sin(2*pi*{x})*sin(0.03+1.97*{y})", lambda x,y,z : sin(2*pi*x)*sin(0.03+1.97*y) )

test_function("{result} = (-pi+4*pi*{x})-sin(-pi+4*pi*{x})*cosh(-2+4*{y})", lambda x,y,z : (-pi+4*pi*x)-sin(-pi+4*pi*x)*cosh(-2+4*y ) )
test_function("{result} = 4*sin(1/2*(-pi+4*pi*{x}))*sinh((-2+4*{y})/2)", lambda x,y,z : 4*sin(1./2.*(-pi+4*pi*x))*sinh((-2+4*y)/2) )
test_function("{result} = 1-cos(-pi+4*pi*{x})*cosh(-2+4*{y})", lambda x,y,z : 1-cos(-pi+4*pi*x)*cosh(-2+4*y) )

test_function("{result} = minimum(cos({x}), sin({y}))", lambda x, y,z : min(cos(x),sin(y)) )
test_function('{result}= %20.16g*log( ({x}/%20.16g)*(%20.16g/{y})**%20.16g ) '%(0.5,10.,1.,1.4),
              lambda x,y,z : 0.5*log( (x/10.)*(1./y)**1.4 ))
test_function('{result} = 0.3*{x}/{y} + 0.3*{z}/{y}', lambda x,y,z: 0.3*x/y + 0.3*z/y)
test_function('{result} = 0.5*({y}-0.3)**2*{y}', lambda x,y,z: 0.5*(y-0.3)**2*y)
