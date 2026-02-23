# - newArbitraryGridMotion (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newArbitraryGridMotion(name='Motion', value='DeformingGrid'); Internal.printTree(n)
#>> ['Motion',None,[1 son],'ArbitraryGridMotion_t']
#>>    |_['ArbitraryGridMotion',array('DeformingGrid',dtype='|S1'),[0 son],'ArbitraryGridMotionType_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newArbitraryGridMotion(name='Motion', value='NonDeformingGrid', parent=z)
