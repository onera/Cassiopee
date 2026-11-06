# Usage: kpython -n3 allreduceDictPT.py
# - allreduceDict (pyTree) -
import Converter.Internal as Internal
import Converter.Mpi as Cmpi

# Data set 1
if Cmpi.rank == 0: bcNames = ["inlet", "wall"]
elif Cmpi.rank == 1: bcNames = ["wall", "inlet", "outlet", "symmetry"]
else: bcNames = []

# Data set 2
if Cmpi.rank == 0: bcNamesDict = {"inlet": [1., 3.], "wall": [5., 6., 10.]}
elif Cmpi.rank == 1: bcNamesDict = {"wall": [7.], "inlet": [2.], "outlet": [], "symmetry": [20., 21.]}
else: bcNamesDict = {"wall": [9., 8.], "outlet": [11.], "symmetry": [20., 21., 22.]}


# -- method 0: Return all unique values, then reduce with op
d = Cmpi.allreduceDict(bcNamesDict, op=Cmpi.SUM)
if Cmpi.master: print(d)
# >> 0: {'inlet': 6.0, 'wall': 45.0, 'outlet': 11.0, 'symmetry': 63.0}


# -- method 1: Return all values, then reduce with op
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d = Cmpi.allreduceDict(d, method=1, op=Cmpi.SUM)
if Cmpi.master: print(d)
# >> 0: {'inlet': 2, 'wall': 2, 'outlet': 1, 'symmetry': 1}


# -- method 2: Return all values, add fillValue if missing, then reduce with op
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d1 = Cmpi.allreduceDict(d, method=2, fillValue=False, op=Cmpi.LOR)
d2 = Cmpi.allreduceDict(d, method=2, fillValue=False, op=Cmpi.LAND)
if Cmpi.master: print(d1)  # all BC names are present on at least one rank
# >> 0: {'inlet': True, 'wall': True, 'outlet': True, 'symmetry': True}
if Cmpi.master: print(d2)  # no BC names are present on all ranks
# >> 0: {'inlet': False, 'wall': False, 'outlet': False, 'symmetry': False}
