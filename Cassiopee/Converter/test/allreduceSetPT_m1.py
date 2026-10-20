# Usage: kpython -n3 allreduceSetPT_m1.py
# - allreduceSet (pyTree) -
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

# Data set
if Cmpi.rank == 0: bcNames = {"inlet", "wall"}
elif Cmpi.rank == 1: bcNames = {"wall", "inlet", "outlet"}
else: bcNames = {"symmetry"}


# Number of unique boundary conditions
bcNames = Cmpi.allreduceSet(bcNames)
nbcs = len(bcNames)
if Cmpi.master: test.testO(nbcs, 1)
