# - addOutput (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
bcs = Internal.getNodesFromName(a, 'wall')
dict = {}
dict['var'] = 'convflux_rou convflux_rov convflux_row'
dict['loc'] = 1
dict['fluxcoef'] = 1.
dict['period'] = 10
for b in bcs: elsAProfile._addOutput(b,dict,name='')
test.testT(a,1)
#
a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._addBC2Zone(a, 'wall', 'BCWall', 'imin')
bcs = Internal.getNodesFromName(a, 'wall*')
for b in bcs: elsAProfile._addOutput(b,dict,name='#100')
test.testT(a,2)
#
bcs = Internal.getNodesFromName(a, 'wall*')
dict={}
dict['var'] = 'convflux_rou convflux_rov convflux_row'
dict['loc'] = 0
dict['fluxcoef'] = 1.
dict['period'] = 100
for b in bcs: elsAProfile._addOutput(b,dict,name='#100',update=True)
test.testT(a,3)
#
bcs = Internal.getNodesFromName(a, 'wall*')
dict={}
dict['var'] = 'convflux_rou convflux_rov convflux_row'
dict['loc'] = 1
dict['fluxcoef'] = 1.
dict['period'] = 10
for b in bcs: elsAProfile._addOutput(b,dict,name='#100',update=False)
test.testT(a,4)
