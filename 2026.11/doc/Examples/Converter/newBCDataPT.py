# - newBCData (pyTree) -
import Converter.Internal as Internal

# Create a BC data node
n = Internal.newBCData(name='BCData'); Internal.printTree(n)
#>> ['BCData',None,[0 son],'BCData_t']

# Attach it to a parent node
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined')
Internal.newBCData('BCNeumann', parent=d); Internal.printTree(d)
#>> ['BCDataSet',array('UserDefined',dtype='|S1'),[1 son],'BCDataSet_t']
#>>    |_['BCNeumann',None,[0 son],'BCData_t']
