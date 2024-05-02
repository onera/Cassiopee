# - convertFile2SkeletonTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Converter.Internal as Internal

if Cmpi.rank == 0:
    a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
    t = C.newPyTree(['Base', a])
    C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()

t1 = Cmpi.convertFile2SkeletonTree('in.cgns'); Internal.printTree(t1)
#>> ['CGNSTree',None,[2 sons],'CGNSTree_t']
#>>   |_['CGNSLibraryVersion',array([3.0999999046325684],dtype='float64'),[0 son],'CGNSLibraryVersion_t']
#>>   |_['Base',array(shape=(2,),dtype='int32',order='F'),[1 son],'CGNSBase_t']
#>>       |_['cart',array(shape=(3, 3),dtype='int32',order='F'),[2 sons],'Zone_t']
#>>           |_['ZoneType',array('Structured',dtype='|S1'),[0 son],'ZoneType_t']
#>>           |_['GridCoordinates',None,[3 sons],'GridCoordinates_t']
#>>               |_['CoordinateX',None,[0 son],'DataArray_t']
#>>               |_['CoordinateY',None,[0 son],'DataArray_t']
#>>               |_['CoordinateZ',None,[0 son],'DataArray_t']
