import KCore

# Test Fld
#KCore.tester()

# Test Generator
#import Generator.PyTree as G
#import Converter.PyTree as C
#a = G.cart((0,0,0), (1,1,1), (10,10,10))
#p = KCore.tester(a, 'GridCoordinates/CoordinateX')
#print(p)

import Generator as G
import Converter as C
import numpy

ntry = 10000000
ntry = 1
ns = 500
# Array1 struct
#a = G.cart((0,0,0), (1,1,1), (ns,ns,ns))
#for i in xrange(ntry): KCore.tester(a)
# Array1 BE
#a = G.cartHexa((0,0,0), (1,1,1), (ns,ns,ns))
#for i in xrange(ntry): KCore.tester(a)
# Array1 NGON
#a = G.cartNGon((0,0,0), (1,1,1), (ns,ns,ns))
#for i in xrange(ntry): KCore.tester(a)

#C.convertArrays2File(a, 'out.plt')

# Array2 struct
#a = ['x,y,z', [numpy.zeros(ns*ns*ns),numpy.zeros(ns*ns*ns),numpy.zeros(ns*ns*ns)], ns,ns,ns]
#for i in xrange(ntry): KCore.tester(a)
a = KCore.tester()
print(a)
# a est directement transformable en zones
