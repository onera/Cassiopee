
# Driver utilities
import OCC

# identify points and return maps
def mapClosePoints(points, decimals=5):
    seen = {}
    old2New = {}
    new2Old = {}; j = 0
    for i, pt in enumerate(points):
        # Round coordinates to group close points
        key = tuple(round(c, decimals) for c in pt)
        if key in seen:
            old2New[i] = seen[key]
        else:
            seen[key] = j
            old2New[i] = j
            new2Old[j] = i; j += 1
    print('seen', seen)
    return old2New, new2Old

# convert a hook to a python CBL file
# IN: din: input dictionary
def convertCAD2CBL(hook, fileName, din):
    """Converter hook to python CBL file."""

    # open output file name
    __LOG__ = open(fileName, 'w')
    l = "import Roms.Driver as D\n"
    l += "T1  = D.Part('part1')\n"

    # create a dummy param
    l += "v1 = T1.Scalar('v1', 0.)\n"
    l += "v1.range = [0., 10.]\n"
    __LOG__.write(l)

    # query[no] -> query egde data
    nbEdges = OCC.getNbEdges(hook)
    print(nbEdges)
    query = {}
    for i in range(nbEdges):
        d = OCC.queryEdge(hook, i+1)
        query[i+1] = d

    # entity: edge query dictionary by din entries (flat list)
    # A SUPPRIMER
    entity = {}
    for k in din: # entry
        print("Info: convertCAD2CBL: reading %s."%k)
        if 'loft' in k or 'loftSet' in k or 'sweep' in k: # list of list of edges
            entity[k] = []
            for s in din[k]: # sketch
                for e in s: # entities
                    d = query[e]
                    entity[k].append(d)
        elif 'fill' in k or 'sketch' in k: # one list of edges
            entity[k] = []
            for e in din[k]: # edge/entity
                d = query[e]
                entity[k].append(d)

    # get all points in entities with a number
    points = []
    for i in query:
        d = query[i]
        dtype = d['type']
        if dtype == "Line":
            P0 = d['P0']
            points.append(P0)
            P1 = d['P1']
            points.append(P1)
        elif dtype == "Circle":
            center = d['center']
            points.append(center)
        elif dtype == "Arc":
            center = d['center']
            points.append(center)
        elif dtype == "Spline1":
            pts = d['Points']
            for p in range(pts.shape[1]):
                x = float(pts[0,p])
                y = float(pts[1,p])
                z = float(pts[2,p])
                points.append((x,y,z))
        elif dtype == "Ellipse":
            center = d['center']
            points.append(center)

    # Give a name to points, reidentify
    old2New, new2Old = mapClosePoints(points, 5)
    #print(points)
    #print("old2New: ", old2New)
    #print("new2Old: ", new2Old, flush=True)

    # write points to log
    l = ''
    for i in new2Old:
        p = points[new2Old[i]]
        l += f"P{i} = T1.Point('P{i}', ({p[0]}, {p[1]}, {p[2]}))\n"
    __LOG__.write(l)

    # write entities to log
    l = ''
    ip = 0 # old point index
    for i in query:
        d = query[i]
        dtype = d['type']
        if dtype == "Line":
            ip0 = old2New[ip]
            ip1 = old2New[ip+1]
            l += f"line{i} = T1.Line('line{i}', P{ip0}, P{ip1})\n"
            ip += 2
        elif dtype == "Circle":
            ip0 = old2New[ip]
            radius = d['radius']
            axis = d['axis']
            l += f"circle{i} = T1.Circle('circle{i}', P{ip0}, {radius})\n"
            ip += 1
        elif dtype == "Arc":
            ip0 = old2New[ip]
            radius = d['radius']
            axis = d['axis']
            l += f"arc{i} = T1.Circle('circle{i}', P{ip0}, {radius})\n"
            ip += 1
        elif dtype == "Spline1":
            pts = d['Points']
            npts = pts.shape[1]
            lp = ''
            for p in range(npts):
                ip0 = old2New[ip+p]
                lp += f"P{ip0}, "
            lp = lp[:-2]
            l += f"spline{i} = T1.Spline1('spline{i}', [{lp}])\n"
            ip += npts
        elif dtype == "Ellipse":
            ip0 = old2New[ip]
            axis = d['axis']
            R1 = d['R1']
            R2 = d['R2']
            l += f"ellipse{i} = T1.Ellipse('ellipse{i}', P{ip0}, {R1}, {R2})\n"
            ip += 1

    __LOG__.write(l)

    # write sketch to log
    isk = 0 # sketch index
    l = ''
    for k in din:
        if 'loft' in k or 'loftSet' in k or 'sweep' in k: # list of list of edges
            for edges in din[k]: # sketch
                lp = f"sketch{isk} = T1.Sketch('sketch{isk}', ["
                for i in edges: # entity
                    d = query[i]
                    dtype = d['type']
                    if dtype == "Line":
                        lp += f"line{i}, "
                    elif dtype == "Circle":
                        lp += f"circle{i}, "
                    elif dtype == "Arc":
                        lp += f"arc{i}, "
                    elif dtype == "Spline1":
                        lp += f"spline{i}, "
                    elif dtype == "Ellipse":
                        lp += f"ellipse{i}, "
                lp = lp[:-2]+'])\n'
                l += lp
                isk += 1
        else: # list of edges (one sketch)
            lp = f"sketch{isk} = T1.Sketch('sketch{isk}', ["
            for i in din[k]:
                d = query[i]
                dtype = d['type']
                if dtype == "Line":
                    lp += f"line{i}, "
                elif dtype == "Circle":
                    lp += f"circle{i}, "
                elif dtype == "Arc":
                    lp += f"arc{i}, "
                elif dtype == "Spline1":
                    lp += f"spline{i}, "
                elif dtype == "Ellipse":
                    lp += f"ellipse{i}, "
            lp = lp[:-2]+'])\n'
            l += lp
            isk += 1
    __LOG__.write(l)

    # write surface to log
    l = ''; isk = 0
    for k in din:
        if 'loft' in k:
            ns = len(din[k])
            l += f"{k} = T1.Loft('{k}', ["
            for i in range(ns):
                l += f"sketch{i+isk}, "
            isk += ns
            l = l[:-2]
            l += "])\n"
        elif 'sweep' in k:
            ns = len(din[k])
            l += f"{k} = T1.Sweep('{k}', ["
            for i in range(ns):
                l += f"sketch{i+isk}, "
            isk += ns
            l = l[:-2]
            l += "])\n"
        elif 'loftSet' in k:
            ns = len(din[k])
            l += f"{k} = T1.LoftSet('{k}', ["
            for i in range(ns):
                l += f"sketch{i+isk}, "
            isk += ns
            l = l[:-2]
            l += "])\n"

        elif 'fill' in k:
            l += f"{k} = T1.Fill('{k}', sketch{s})"
            isk += 1
    __LOG__.write(l)

    # write end of file
    l = "T1.solve()\n"
    l += "T1.instantiate({'v1':1.})\n"
    l += "T1.writeCAD('outCBL.step')\n"
    l += "D.DRIVER.plot()\n"
    __LOG__.write(l)

    __LOG__.close()

    # display output to screen (debug)
    f = open(fileName, "r")
    lines = f.readlines()
    for l in lines: print(l)
    f.close()
