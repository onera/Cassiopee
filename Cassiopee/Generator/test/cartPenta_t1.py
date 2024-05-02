# - cartPenta (array) -
import Generator as G
import KCore.test as test

a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
test.testA([a],1)
