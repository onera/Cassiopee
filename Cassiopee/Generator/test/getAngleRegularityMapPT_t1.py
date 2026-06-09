# - getGridSkewnessMap (pyTree) -
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# Test 1D structure
a = G.cart((0,0,0), (1,1,1), (10,1,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,0,0), (0.,1.,0.), 0.5, 0.5)
t = G.getGridSkewnessMap(a)
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 1)
t = G.getGridSkewnessMap(t, normalized=True)
Internal._rmNodesByName2(t, 'regularityAngle')
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 11)

# Test 2D structure
a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,0), (1.,1.,0.), 0.5, 0.5)
t = G.getGridSkewnessMap(a)
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 2)
t = G.getGridSkewnessMap(t, normalized=True)
Internal._rmNodesByName2(t, 'regularityAngle')
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 21)

# Test 3D structure
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,9), (1.,1.,1.), 0.5, 0.5)
t = G.getGridSkewnessMap(a)
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 3)
t = G.getGridSkewnessMap(t, normalized=True)
Internal._rmNodesByName2(t, 'regularityAngle')
Internal._renameNode(t, 'gridSkewness', 'regularityAngle') # backward compatibility with test case references
test.testT(t, 31)