# Usage: kpython -n3 allreduceDictPT_m2.py
# - allreduceDict (pyTree) -
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


# -- method 0: Return all unique values, then reduce with op
d = Cmpi.allreduceDict(bcNamesDict, op=Cmpi.MAX)
if Cmpi.master: print(d); #test.testO(d, 1)
# {'inlet': 3.0, 'wall': {'flap1': {'upper': 10.0, 'lower': 10.2}, 'flap2': 10.5}, 'outlet': 11.0, 'symmetry': 22.0}

# -- method 1: Return all values, then reduce with op
d = Cmpi.allreduceDict(bcNamesDict, method=1, op=Cmpi.SUM)
if Cmpi.master: print(d); #test.testO(d, 2)
# {'inlet': 7.0, 'wall': {'flap1': {'upper': 45.0, 'lower': 46.2}, 'flap2': 48.0}, 'outlet': 11.0, 'symmetry': 104.0}

# -- method 2: Return all values, add fillValue if missing, then reduce with op
d = Cmpi.allreduceDict(bcNamesDict, method=2, fillValue=-1., op=Cmpi.SUM)
if Cmpi.master: print(d); #test.testO(d, 3)
# {'inlet': 6.0, 'wall': {'flap1': {'upper': 45.0, 'lower': 46.2}, 'flap2': 48.0}, 'outlet': 10.0, 'symmetry': 103.0}
