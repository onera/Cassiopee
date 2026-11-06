# - newDataClass (pyTree) -
import Converter.Internal as Internal

# Create a DataClass node
n = Internal.newDataClass('Dimensional'); Internal.printTree(n)
#>> ['DataClass',array('Dimensional',dtype='|S1'),[0 son],'DataClass_t']

# Attach it to a parent node
d = Internal.newDiscreteData('DiscreteData')
Internal.newDataClass('Dimensional', parent=d)
