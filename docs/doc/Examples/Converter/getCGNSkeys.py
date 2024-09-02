# - getCGNSKey (pyTree) -
import Converter.elsAProfile as elsAProfile

for elsAkey in elsAProfile.keyselsA2CGNS:
    print('elsA %s corresponds to CGNS name %s'%(elsAkey,elsAProfile.getCGNSkeys(elsAkey)))
