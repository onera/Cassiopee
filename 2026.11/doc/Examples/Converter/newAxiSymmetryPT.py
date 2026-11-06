# - newAxiSymmetry (pyTree) -
import Converter.Internal as Internal

# Create an axisymmetry node
n = Internal.newAxiSymmetry(referencePoint=[0.,0.,0.], axisVector=[0.,0.,0.]); Internal.printTree(n)
#>> ['AxiSymmetry',None,[2 sons],'AxiSymmetry_t']
#>>    |_['AxiSymmetryReferencePoint',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']
#>>    |_['AxiSymmetryAxisVector',array(shape=(3,),dtype='float64',order='F'),[0 son],'DataArray_t']

# Attach it to base
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newAxiSymmetry([0.,0.,0.], [0.,0.,0.], parent=b)
