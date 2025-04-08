# - connectMatch (pyTree) -
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.elsAProfile as elsAProfile

a = G.cylinder((0.,0.,0.), 0.1, 1., 0., 90., 5., (11,11,11))
t = C.newPyTree(['Base',a])
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           translation=[0.,0.,5.])
t = X.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.],
                           rotationAngle=[0.,0.,90.])
Internal.printTree(Internal.getNodeFromName(t, 'match1_0'))
#>> ['match1_0',array('cylinder',dtype='|S1'),[4 sons],'GridConnectivity1to1_t']
#>>    |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['PointRangeDonor',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['Transform',array(shape=(3,),dtype='int32',order='F'),[0 son],'"int[IndexDimension]"']
#>>    |_['GridConnectivityProperty',None,[1 son],'GridConnectivityProperty_t']
#>>        |_['Periodic',None,[3 sons],'Periodic_t']
#>>            |_['RotationCenter',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>            |_['RotationAngle',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>            |_['Translation',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
tp = elsAProfile.adaptPeriodicMatch(t)
Internal.printTree(Internal.getNodeFromName(tp, 'match1_0'))
#>> ['match1_0',array('cylinder',dtype='|S1'),[5 sons],'GridConnectivity1to1_t']
#>>    |_['PointRange',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['PointRangeDonor',array(shape=(3, 2),dtype='int32',order='F'),[0 son],'IndexRange_t']
#>>    |_['Transform',array(shape=(3,),dtype='int32',order='F'),[0 son],'"int[IndexDimension]"']
#>>    |_['GridConnectivityProperty',None,[1 son],'GridConnectivityProperty_t']
#>>    |   |_['Periodic',None,[3 sons],'Periodic_t']
#>>    |       |_['RotationCenter',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |       |_['RotationAngle',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |       |_['Translation',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['.Solver#Property',None,[7 sons],'UserDefinedData_t']
#>>        |_['type',array('join',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['jtopo',array('periodic',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['jtype',array('match',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['ptype',array('tra',dtype='|S1'),[0 son],'DataArray_t']
#>>        |_['xtran',array([0.0],dtype='float64'),[0 son],'DataArray_t']
#>>        |_['ytran',array([0.0],dtype='float64'),[0 son],'DataArray_t']
#>>        |_['ztran',array([5.0],dtype='float64'),[0 son],'DataArray_t']
C.convertPyTree2File(tp, 'out.cgns')
