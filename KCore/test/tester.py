import KCore


a = KCore.tester()
print a
#import Converter as C
#b = C.converter.convertHO2LO(a, 0)
#print b
#C.convertArrays2File(b, 'out.plt')
#import sys; sys.exit()

import Converter.PyTree as C
import Converter.Internal as Internal
import Converter
a = C.convertArrays2ZoneNode('triangle', [a])
C.convertPyTree2File(a, 'out.cgns')

#fc = C.getFields(Internal.__GridCoordinates__, a, api=2)[0]
#print fc
#fcp = Converter.converter.convertHO2LO(fc, 0)
#print fcp
#C.setFields([fcp], a, 'nodes', True)

#print a
#C.convertPyTree2File(a, 'out2.cgns')

b = C.convertHO2LO(a, mode=0)

print Internal.getZoneDim(b)
C.convertPyTree2File(b, 'out2.cgns')

import sys; sys.exit()

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
#for i in range(ntry): KCore.tester(a)
# Array1 BE
#a = G.cartHexa((0,0,0), (1,1,1), (ns,ns,ns))
#for i in range(ntry): KCore.tester(a)
# Array1 NGON
#a = G.cartNGon((0,0,0), (1,1,1), (ns,ns,ns))
#for i in range(ntry): KCore.tester(a)

#C.convertArrays2File(a, 'out.plt')

# Array2 struct
#a = ['x,y,z', [numpy.zeros(ns*ns*ns),numpy.zeros(ns*ns*ns),numpy.zeros(ns*ns*ns)], ns,ns,ns]
#for i in range(ntry): KCore.tester(a)
a = KCore.tester()
print(a)
# a est directement transformable en zones
