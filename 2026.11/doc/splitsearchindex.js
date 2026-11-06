
file1 = open("search1.txt")
ret = file1.readline()
ret = ret.split(",")
out = []
for r in ret:
    r = r.split(":")
    out += r
file1.close()

for o in out:
    print(o)


file2 = open("split1.txt", "w")
for r in out:
    file2.write(r)
    file2.write('\n')
file2.close()
