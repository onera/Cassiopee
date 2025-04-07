# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivity1to1(name='Match', donorName='blk1', pointRange=[1,1,1,33,1,69], pointRangeDonor=[51,51,1,33,1,69], transform=[1,2,3]); Internal.printTree(n)
#>> ['Match',array('blk1',dtype='|S1'),[2 sons],'GridConnectivity1to1_t']
#>>    |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['PointRangeDonor',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['Transform',array(shape=(3,),dtype='int32',order='F'),[0 son],'"int[IndexDimension]"']

# Attach it to a parent node
d = Internal.newZoneGridConnectivity(name='ZoneGridConnectivity')
Internal.newGridConnectivity1to1(name='Match', donorName='blk1', pointRange=[1,1,1,33,1,69], pointRangeDonor=[51,51,1,33,1,69], transform=None, parent=d)
