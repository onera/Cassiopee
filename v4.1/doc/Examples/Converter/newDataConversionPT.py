# - newDataConversion (pyTree) -
import Converter.Internal as Internal

# Data(raw) = Data(nondimensional)*ConversionScale + ConversionOffset

# Create a DataConversion node
n = Internal.newDataConversion(conversionScale=1., conversionOffset=0.); Internal.printTree(n)
#>> ['DataConversion',array(shape=(2,),dtype='float64',order='F'),[0 son],'DataConversion_t']

# Attach it to a parent node
d = Internal.newDiscreteData('DiscreteData')
Internal.newDataClass('NondimensionalParameter', parent=d)
Internal.newDataConversion(conversionScale=1., conversionOffset=0., parent=d)
