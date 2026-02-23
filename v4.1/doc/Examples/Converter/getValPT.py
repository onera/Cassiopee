# - getVal (pyTree) -
import Converter.Internal as Internal

# getVal returns always numpys
z = Internal.newZone('Zone', zsize=[[10, 2, 0]], ztype='Structured')
print(Internal.getVal(z))
# >> [[10  2  0]]
n = Internal.getNodeFromName(z, 'ZoneType')
print(Internal.getVal(n))
#>> [b'S' b't' b'r' b'u' b'c' b't' b'u' b'r' b'e' b'd']
