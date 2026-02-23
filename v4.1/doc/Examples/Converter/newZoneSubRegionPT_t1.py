# - newZoneSubRegion (pyTree) -
import Converter.Internal as Internal
import KCore.test as test

# Create a node linked to a BC
n = Internal.newZoneSubRegion(name='SubRegionBC', bcName='BC')
test.testO(n, 1)

# Create a node linked to a GridConnectivity
n = Internal.newZoneSubRegion(name='SubRegionGC', gcName='GC')
test.testO(n, 2)

# Create a node
n = Internal.newZoneSubRegion(name='SubRegion', pointList=[1, 2, 3, 4], gridLocation='FaceCenter')
test.testO(n, 3)

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10, 2, 2]], ztype='Structured')
Internal.newZoneSubRegion(name='SubRegion', pointRange=[1,10,2,2,1,1], gridLocation='FaceCenter', parent=z)
test.testO(n, 4)
