# - array (array) -
import Converter as C

# Structured
b = C.array('x,y,z', 12, 9, 12); print(b)
#>> ['x,y,z', array(...), 12, 9, 12]

# Unstructured
a = C.array('x,y,z', 12, 9, 'QUAD'); print(a)
#>> ['x,y,z', array(...), array(..., dtype=int32), 'QUAD']
