# Usage: kpython -n3 allgatherDictPT_m1.py
# - allgatherDict (pyTree) -
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

# Data set 1
if Cmpi.rank == 0: bcNames = ["inlet", "wall"]
elif Cmpi.rank == 1: bcNames = ["wall", "inlet", "outlet", "symmetry"]
else: bcNames = []

# Data set 2
if Cmpi.rank == 0: bcNamesDict = {"inlet": [1., 3.], "wall": [5., 6., 10.]}
elif Cmpi.rank == 1: bcNamesDict = {"wall": [7.], "inlet": [2.], "outlet": [], "symmetry": [20., 21.]}
else: bcNamesDict = {"wall": [9., 8.], "outlet": [11.], "symmetry": [20., 21., 22.]}


# -- method 0: Return all unique values
# Values are a non-iterable (int)
d = {name: Cmpi.rank for name in bcNames}
d = Cmpi.allgatherDict(d)
if Cmpi.master: test.testO(d, 1)
# {'inlet': [0, 1], 'wall': [0, 1], 'outlet': [1], 'symmetry': [1]}

# Values are an iterable (list)
d = {Cmpi.rank: bcNames}
d = Cmpi.allgatherDict(d)
if Cmpi.master: test.testO(d, 2)
# {0: ['inlet', 'wall'], 1: ['wall', 'inlet', 'outlet', 'symmetry'], 2: []}

d = Cmpi.allgatherDict(bcNamesDict)
if Cmpi.master: test.testO(d, 3)
# {'inlet': [1.0, 3.0, 2.0], 'wall': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'outlet': [11.0], 'symmetry': [20.0, 21.0, 22.0]}


# -- method 1: Return all values
# Values are a non-iterable (bool)
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d = Cmpi.allgatherDict(d, method=1)
if Cmpi.master: test.testO(d, 4)
# {'inlet': [True, True], 'wall': [True, True], 'outlet': [True], 'symmetry': [True]}

# Values are an iterable (list)
d = Cmpi.allgatherDict(bcNamesDict, method=1)
if Cmpi.master: test.testO(d, 5)
# {'inlet': [1.0, 3.0, 2.0], 'wall': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'outlet': [11.0], 'symmetry': [20.0, 21.0, 20.0, 21.0, 22.0]}


# -- method 2: Return all values, add fillValue if missing
# Values are a non-iterable (bool), default fillValue (None)
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d = Cmpi.allgatherDict(d, method=2)
if Cmpi.master: test.testO(d, 6)
# {'inlet': [True, True, None], 'wall': [True, True, None], 'outlet': [None, True, None], 'symmetry': [None, True, None]}

# Values are a non-iterable (bool), custom fillValue
isBCPresent = True
d = {name: isBCPresent for name in bcNames}
d = Cmpi.allgatherDict(d, method=2, fillValue=False)
if Cmpi.master: test.testO(d, 7)
# {'inlet': [True, True, False], 'wall': [True, True, False], 'outlet': [False, True, False], 'symmetry': [False, True, False]}

# Values are an iterable (tuple), custom fillValue (0.)
d = {name: (1.,) for name in bcNames}
d = Cmpi.allgatherDict(d, method=2, fillValue=0.)
if Cmpi.master: test.testO(d, 8)
# {'inlet': [1.0, 1.0, 0.0], 'wall': [1.0, 1.0, 0.0], 'outlet': [0.0, 1.0, 0.0], 'symmetry': [0.0, 1.0, 0.0]}


# -- method 3: Return the single value, the first of the lowest rank
# Values are a non-iterable (int)
d = {name: Cmpi.rank for name in bcNames}
d = Cmpi.allgatherDict(d, method=3)
if Cmpi.master: test.testO(d, 9)
# {'inlet': 0, 'wall': 0, 'outlet': 1, 'symmetry': 1}

# Values are an iterable (list)
d = Cmpi.allgatherDict(bcNamesDict, method=3)
if Cmpi.master: test.testO(d, 9)
# {'inlet': [1.0, 3.0], 'wall': [5.0, 6.0, 10.0], 'outlet': [], 'symmetry': [20.0, 21.0]}
