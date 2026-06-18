# - getCellSkewnessMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import KCore.test as test

# Test 3D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
a = C.addBC2Zone(a, 'wall1','BCWall','jmin')
a = C.addBC2Zone(a, 'wall2','BCWall','kmin')
a = C.addBC2Zone(a, 'match1','BCMatch','imin',a,'imax',[1,2,3])
a = C.addBC2Zone(a, 'match2','BCMatch','imax',a,'imin',[1,2,3])
a = C.fillEmptyBCWith(a,'overlap','BCOverlap')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
a = G.getCellSkewnessMap(a)
Internal._renameNode(a, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
t = C.newPyTree(['Base']); t[2][1][2].append(a)
test.testT(t, 1)
t = G.getCellSkewnessMap(t, normalized=True)
Internal._rmNodesByName(t, 'orthogonality')
Internal._renameNode(t, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
test.testT(t, 11)

# Test 2D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.addBC2Zone(a, 'wall1','BCWall','imin')
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getCellSkewnessMap(a)
Internal._renameNode(t, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
test.testT(t, 2)
t = G.getCellSkewnessMap(t, normalized=True)
Internal._rmNodesByName2(t, 'orthogonality')
Internal._renameNode(t, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
test.testT(t, 21)

# Test 1D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,10,1))
a = C.initVars(a,'Density',1.); a = C.initVars(a,'centers:cellN',1.)
t = G.getCellSkewnessMap(a)
Internal._renameNode(t, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
test.testT(t, 3)
t = G.getCellSkewnessMap(t, normalized=True)
Internal._rmNodesByName2(t, 'orthogonality')
Internal._renameNode(t, 'cellSkewness', 'orthogonality') # backward compatibility with test case references
test.testT(t, 31)
