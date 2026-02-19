# - perlinNoise (pyTree) -
import Post.PyTree as P
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# struct
a = G.cart((0,0,0), (1,1,1), (10,10,1))
P._perlinNoise(a)
test.testT(a, 1)

# hexa
a = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
P._perlinNoise(a)
test.testT(a, 2)
