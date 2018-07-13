# - printTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a) # can print a zone (or any node or any list of nodes)
#>> ['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>   |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']
#>>   |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>>       |_['CoordinateX',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>       |_['CoordinateY',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>       |_['CoordinateZ',array(shape=(10, 10, 10),dtype='float64',order='F'),[0 son],'DataArray_t']

t = C.newPyTree(['Base',a])
Internal.printTree(t) # can print a tree
#>> ['CGNSTree',None,[2 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.1],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   ..

Internal.printTree(t, file='toto.txt') # in a file
f = open('toto.txt', 'a')
Internal.printTree(t, stdOut=f) # in a file object
Internal.printTree(t, editor='emacs') # with an editor
