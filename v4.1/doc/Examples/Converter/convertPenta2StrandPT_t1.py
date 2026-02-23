# - convertPenta2Strand (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cartPenta((0,0,0), (1,1,1), (3,3,3))

C._convertPenta2Strand(a)

test.testT(a, 1)
