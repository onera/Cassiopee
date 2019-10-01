import Generator as G
import Converter as C
import Converter.expression as expr
import time

N = 300

print("------------------------  initVars -----------------------------")
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
tb1 = time.time()
C._initVars(crds, "{norm} = {x}*{x}+{y}*{y}+{z}*{z}-{x}*{y}-{x}*{z}-{y}*{z}")
tb2 = time.time()
print(crds[0])
print(crds[1][:, 0:10])
print("....")
print(crds[1][:, -10:])
print("Temps pris par initVars : {} secondes".format((tb2 - tb1)))

print("----------------------- expression  ----------------------------")
a = expr.ast("{norm} = {x}*{x}+{y}*{y}+{z}*{z}-{x}*{y}-{x}*{z}-{y}*{z}")
print(a)
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
ta1 = time.time()
C._addVars(crds, 'norm')
a.run(crds)
ta2 = time.time()
print(crds[0])
print(crds[1][:, 0:10])
print("....")
print(crds[1][:, -10:])
print("Temps pris par expression : {}".format((ta2 - ta1)))

print("----------------------- numpy  ----------------------------")
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
x = crds[1][0, :]
y = crds[1][1, :]
z = crds[1][2, :]
tc1 = time.time()
nrm = x * x + y * y + z * z - (x * y + x * z + y * z)
tc2 = time.time()
print(nrm[: 10])
print(nrm[-10:])
print("Temps pris par numpy : {}".format(tc2 - tc1))

del x
del y
del z
del nrm

print("------------------------  initVars -----------------------------")
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
tb1 = time.time()
C._initVars(crds, "{norm} = ({x}>0.7)*{x}+({y}<0.5)*{y}+({z}>0.1)*{z}-{x}*{y}-{x}*{z}-{y}*{z}")
tb2 = time.time()
print(crds[0])
print(crds[1][:, 0:10])
print("....")
print(crds[1][:, -10:])
print("Temps pris par initVars : {} secondes".format((tb2 - tb1)))

print("----------------------- expression  ----------------------------")
a = expr.ast("{norm} = ({x}>0.7)*{x}+({y}<0.5)*{y}+({z}>0.1)*{z}-{x}*{y}-{x}*{z}-{y}*{z}")
print(a)
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
ta1 = time.time()
C._addVars(crds, 'norm')
a.run(crds)
ta2 = time.time()
print(crds[0])
print(crds[1][:, 0:10])
print("....")
print(crds[1][:, -10:])
print("Temps pris par expression : {}".format((ta2 - ta1)))

print("----------------------- numpy  ----------------------------")
crds = G.cart((0, 0, 0), (1, 1, 1), (N, N, N), api=1)
x = crds[1][0, :]
y = crds[1][1, :]
z = crds[1][2, :]
tc1 = time.time()
nrm = (x > 0.7) * x + (y < 0.5) * y + (z > 0.1) * z - x * y - x * z - y * z  # x * x + y * y + z * z - (x * y + x * z + y * z)
tc2 = time.time()
print(nrm[: 10])
print(nrm[-10:])
print("Temps pris par numpy : {}".format(tc2 - tc1))
