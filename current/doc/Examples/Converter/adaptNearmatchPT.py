# - connectNearMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 3)); a1[0] = 'cart1'
a2 = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 3)); a2[0] = 'cart2'
a2 = T.oneovern(a2,(1,2,1))
t = C.newPyTree(['Base',a1,a2])
t = X.connectNearMatch(t)
Internal.printTree(Internal.getNodeFromName(t, 'nmatch1_0'))
#>> ['nmatch1_0',array('cart2',dtype='|S1'),[4 sons],'GridConnectivity_t']
#>>    |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['GridConnectivityType',array('Abutting',dtype='|S1'),[0 son],'GridConnectivityType_t']
#>>    |_['PointListDonor',array(shape=(3, 1),dtype='int32',order='C'),[0 son],'IndexArray_t']
#>>    |_['UserDefinedData',None,[3 sons],'UserDefinedData_t']
#>>        |_['PointRangeDonor',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'DataArray_t']
#>>        |_['Transform',array(shape=(3,),dtype='int32',order='F'),[0 son],'DataArray_t']
#>>        |_['NMRatio',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
tp = elsAProfile.adaptNearMatch(t)
Internal.printTree(Internal.getNodeFromName(tp, 'nmatch1_0'))
#>> ['nmatch1_0',array('cart2',dtype='|S1'),[4 sons],'GridConnectivity1to1_t']
#>>    |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['PointRangeDonor',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['Transform',array(shape=(3,),dtype='int32',order='F'),[0 son],'"int[IndexDimension]"']
#>>    |_['.Solver#Property',None,[6 sons],'UserDefinedData_t']
#>>        |_['jtype',array('nearmatch',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['type',array('join',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['matchside',array('fine',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['i_ratio',array([1],dtype='int32'),[0 son],'DataArray_t']
#>>        |_['j_ratio',array([2],dtype='int32'),[0 son],'DataArray_t']
#>>        |_['k_ratio',array([1],dtype='int32'),[0 son],'DataArray_t']
C.convertPyTree2File(tp, 'out.cgns')
