# - isNamePresent (array) -
import KCore
import Converter as C

a = C.array('x,y,z , densitY',1,1,1)
pos = KCore.isNamePresent(a, 'Y')
print pos, 'must be -1'
pos = KCore.isNamePresent(a, 'densitY')
print pos, 'must be 3.'
