# Usage: kpython -n3 allgatherDictPT_m2.py
# - allgatherDict (pyTree) -
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import KCore.test as test

# Data set
if Cmpi.rank == 0:
    bcNamesDict = {
        "inlet": [1., 3.],
        "wall": {
            "flap1": {
                "upper": [5., 6., 10.],
                "lower": [5.2, 6.2, 10.2],
            },
            "flap2": [5.5, 6.5, 10.5]
        }
    }
elif Cmpi.rank == 1:
    bcNamesDict = {
        "wall": {
            "flap1": {
                "upper": [7.],
                "lower": [7.2],
            },
            "flap2": [7.5],
        },
        "inlet": [3.],
        "outlet": [],
        "symmetry": [20., 21.]
    }
else:
    bcNamesDict = {
        "wall": {
            "flap1": {
                "upper": [9., 8.],
                "lower": [9.2, 8.2],
            },
            "flap2": [9.5, 8.5],
        },
        "outlet": [11.],
        "symmetry": [20., 21., 22.]
    }


# -- method 0: Return all unique values
d = Cmpi.allgatherDict(bcNamesDict)
if Cmpi.master: print(d); #test.testO(d, 1)
# {'inlet': [1.0, 3.0], 'wall': {'flap1': {'upper': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'lower': [5.2, 6.2, 10.2, 7.2, 9.2, 8.2]}, 'flap2': [5.5, 6.5, 10.5, 7.5, 9.5, 8.5]}, 'outlet': [11.0], 'symmetry': [20.0, 21.0, 22.0]}

# -- method 1: Return all values
d = Cmpi.allgatherDict(bcNamesDict, method=1)
if Cmpi.master: print(d); #test.testO(d, 2)
# {'inlet': [1.0, 3.0, 3.0], 'wall': {'flap1': {'upper': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'lower': [5.2, 6.2, 10.2, 7.2, 9.2, 8.2]}, 'flap2': [5.5, 6.5, 10.5, 7.5, 9.5, 8.5]}, 'outlet': [11.0], 'symmetry': [20.0, 21.0, 20.0, 21.0, 22.0]}

# -- method 2: Return all values, add fillValue if missing
d = Cmpi.allgatherDict(bcNamesDict, method=2, fillValue=-1.)
if Cmpi.master: print(d); #test.testO(d, 3)
# {'inlet': [1.0, 3.0, 3.0, -1.0], 'wall': {'flap1': {'upper': [5.0, 6.0, 10.0, 7.0, 9.0, 8.0], 'lower': [5.2, 6.2, 10.2, 7.2, 9.2, 8.2]}, 'flap2': [5.5, 6.5, 10.5, 7.5, 9.5, 8.5]}, 'outlet': [-1.0, 11.0], 'symmetry': [-1.0, 20.0, 21.0, 20.0, 21.0, 22.0]}

# -- method 3: Return the single value, the first of the lowest rank
d = Cmpi.allgatherDict(bcNamesDict, method=3)
if Cmpi.master: print(d); #test.testO(d, 4)
# {'inlet': [1.0, 3.0], 'wall': {'flap1': {'upper': [5.0, 6.0, 10.0], 'lower': [5.2, 6.2, 10.2]}, 'flap2': [5.5, 6.5, 10.5]}, 'outlet': [], 'symmetry': [20.0, 21.0]}
