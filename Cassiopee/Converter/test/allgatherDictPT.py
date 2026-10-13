# Usage: kpython -n3 allgatherDictPT.py
# - allgatherDict (pyTree) -
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


# -- method 0: Return all unique values
d = Cmpi.allgatherDict(bcNamesDict)
if Cmpi.master: print(d)
# >> 0: {'inlet': [1.0, 3.0, 2.0], 'wall': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'outlet': [11.0], 'symmetry': [20.0, 21.0, 22.0]}


# -- method 1: Return all values
d = Cmpi.allgatherDict(bcNamesDict, method=1)
if Cmpi.master: print(d)
# >> 0: {'inlet': [1.0, 3.0, 2.0], 'wall': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'outlet': [11.0], 'symmetry': [20.0, 21.0, 20.0, 21.0, 22.0]}


# -- method 2: Return all values, add fillValue if missing
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d = Cmpi.allgatherDict(d, method=2, fillValue=False)
if Cmpi.master: print(d)
# >> 0: {'inlet': [True, True, False], 'wall': [True, True, False], 'outlet': [False, True, False], 'symmetry': [False, True, False]}


# -- method 3: Return the single value, the first of the lowest rank
d = {name: Cmpi.rank for name in bcNames}
d = Cmpi.allgatherDict(d, method=3)
if Cmpi.master: print(d)
# >> 0: {'inlet': 0, 'wall': 0, 'outlet': 1, 'symmetry': 1}
