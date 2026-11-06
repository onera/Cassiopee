# - perlinNoise (pyTree) -
import Post.PyTree as P
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))
P._perlinNoise(a)
C.convertPyTree2File(a, 'out.cgns')
