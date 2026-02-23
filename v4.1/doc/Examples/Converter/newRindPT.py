# - newRind (pyTree) -
import Converter.Internal as Internal

# Create a Rind node (Rind contains the number of ghost cells for a structured block)
n = Internal.newRind([0,0,0,0,1,1]); Internal.printTree(n)
#>> ['Rind',array(shape=(6,),dtype='int32',order='F'),[0 son],'Rind_t']

# Attach it to a parent node
d = Internal.newGridCoordinates()
Internal.newRind([0,0,0,0,1,1], parent=d)
