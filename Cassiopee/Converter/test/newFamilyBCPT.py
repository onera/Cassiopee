# - newFamilyBC (pyTree) -
import Converter.Internal as Internal

# Create a FamilyBC node
n = Internal.newFamilyBC(value='BCWall'); Internal.printTree(n)
#>> ['FamilyBC',array('BCWall',dtype='|S1'),[0 son],'FamilyBC_t']

# Create a FamilyBC node and attach it to Family
b = Internal.newFamily(name='FamInjection')
n = Internal.newFamilyBC(value='UserDefined', parent=b)
