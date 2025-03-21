# - newZoneSubRegion (pyTree) -
import Converter.Internal as Internal

# Create a node linked to a BC
n = Internal.newZoneSubRegion(name='SubRegionBC', bcName='BC'); Internal.printTree(n)
#>> ['SubRegionBC',None,[2 sons],'ZoneSubRegion_t']
#>>    |_['BCRegionName',array('BC',dtype='|S1'),[0 son],'Descriptor_t']
#>>    |_['GridLocation',array('FaceCenter',dtype='|S1'),[0 son],'GridLocation_t']

# Create a node linked to a GridConnectivity
n = Internal.newZoneSubRegion(name='SubRegionGC', gcName='GC'); Internal.printTree(n)
#>> ['SubRegionGC',None,[2 sons],'ZoneSubRegion_t']
#>>    |_['GridConnectivityRegionName',array('GC',dtype='|S1'),[0 son],'Descriptor_t']
#>>    |_['GridLocation',array('FaceCenter',dtype='|S1'),[0 son],'GridLocation_t']

# Create a node
n = Internal.newZoneSubRegion(name='SubRegion', pointList=[1, 2, 3, 4], gridLocation='FaceCenter'); Internal.printTree(n)
#>> ['SubRegion',None,[2 sons],'ZoneSubRegion_t']
#>>    |_['PointList',array(shape=(4,),dtype='int32',order='F'),[0 son],'IndexArray_t']
#>>    |_['GridLocation',array('FaceCenter',dtype='|S1'),[0 son],'GridLocation_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10, 2, 2]], ztype='Structured')
Internal.newZoneSubRegion(name='SubRegion', pointRange=[1,10,2,2,1,1], gridLocation='FaceCenter', parent=z)
