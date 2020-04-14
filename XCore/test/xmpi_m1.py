# - xmpi -
from XCore import xcore
import KCore.test as test
import Converter.Mpi as Cmpi

f = xcore.test_all()
test.testO(f, 1)
