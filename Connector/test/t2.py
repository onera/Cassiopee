# Cassiopee Connector installation test
import sys, os
from elsA import *

# Rm . from PYTHONPATH
try: del sys.path[sys.path.index('')]
except: pass
try: del sys.path[sys.path.index(os.getcwd())]
except: pass

try:
    import Connector.Cassiopee
    print("Connector.Cassiopee correctly installed.")
except Exception as inst:
    print("FAILED:",inst)
    print("FAILED: Connector.Cassiopee badly installed.")
