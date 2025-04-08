# - newGeometryReference (pyTree) -
import Converter.Internal as Internal

# Create a GeometryReference node
n = Internal.newGeometryReference(value='ICEM-CFD', file='MyCAD.tin'); Internal.printTree(n)
#>> ['GeometryReference',None,[2 sons],'GeometryReference_t']
#>>    |_['GeometryFormat',array('ICEM-CFD',dtype='|S1'),[0 son],'GeometryFormat_t']
#>>    |_['GeometryFile',array('MyCAD.tin',dtype='|S1'),[0 son],'GeometryFile_t']

# Create a GeometryReference node and attach it to a family
b = Internal.newFamily(name='FamWall')
n = Internal.newGeometryReference(value='NASA-IGES', file='MyCAD.iges', parent=b)
