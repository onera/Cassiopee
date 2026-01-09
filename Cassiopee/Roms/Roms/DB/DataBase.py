# Manage base of CFD data

# DataBase summary:
# name.db <dir>
#   name.sql (sql data base)
#   ref.cgns (common cgns data)
#   0001.cgns (key associated data)

# vocabulary:
# query(sqlstring) -> return a query (list of full sql data)
# query(point) -> return query for given parametric point
# fetch(query) -> return real data corresponding to query
# parametric point: a dict of parameter names/value (as in Driver walkDOE)

import sqlite3, shutil
import os, numpy, time, datetime
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter
import Compressor.PyTree as Compressor

class DataBase:
    def __init__(self, name, parameters=None, mode="a"):
        # mode
        if mode == 'a': # enable read/write
            self.mode = 'a'
        elif mode == 'r': # read only
            self.mode = 'r'
        elif mode == 'w': # start anew
            self.mode = 'w'
        else: raise ValueError('DataBase: invalid mode.')
        # database name
        if name[-3:] == ".db": name = name[:-3]
        if name[-4:] == ".db/": name = name[:-4]
        self.name = name
        # directory name
        self.dirName = name+'.db'
        if not os.path.exists(self.dirName):
            if mode == 'w':
                os.mkdir(self.dirName)
            elif mode == 'r':
                raise ValueError('DataBase: can not open file for reading.')
            else:
                os.mkdir(self.dirName)
        else:
            if mode == 'w':
                shutil.rmtree(self.dirName)
                os.mkdir(self.dirName)

        # pointer on sql db
        self.db = None
        self.cursor = None
        if os.path.exists("%s/%s.sql"%(self.dirName, self.name)):
            self.open()
            p = self.queryColumnNames()[5:]
            self.parameters = p
        else:
            if parameters is None:
                raise ValueError("DataBase: can not create data base with no parameters.")
            self.parameters = parameters
            self.open()
            self.createTable()
            print("DataBase: creating %s"%self.dirName)
        # column name list
        self.columns = ['id','descp','date','reference','variables']+self.parameters

    # open sql db
    def open(self):
        """Open sql db."""
        self.db = sqlite3.connect('%s/%s.sql'%(self.dirName,self.name))
        self.db.execute("PRAGMA journal_mode=WAL;")
        self.cursor = self.db.cursor()
        return self.cursor

    # close sql db
    def close(self):
        """Close the db."""
        self.db.close()

    # create sql table
    def createTable(self):
        """Create the table."""
        columns = [("id", "INTEGER PRIMARY KEY"),
                   ("descp", "TEXT"), # description
                   ("date", "TEXT"), # date of insertion
                   ("reference", "TEXT"), # reference mesh name
                   ("variables", "TEXT")] # variables list (comma separated)
        for p in self.parameters:
            columns += [("%s"%p, "REAL")]
        columnsSql = ", ".join([f"{col} {dtype}" for col, dtype in columns])
        self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {self.name} ({columnsSql})")

    # create reference (all without fields)
    def registerReference(self, data, name):
        """Write reference file used for coordinates and connectivity."""
        if self.mode == 'r': raise ValueError('register: can not write in read only mode.')
        cgnsName = self.dirName+'/%s'%name+'.cgns'
        tp, _ = Internal.node2PyTree(data)
        tp = Internal.copyRef(tp)
        Internal._rmNodesFromType(tp, 'FlowSolution_t')
        Compressor._compressAll(tp) # lossless
        C.convertPyTree2File(tp, cgnsName)
        return None

    # enable concurrent write to database
    def concurrentExecute(self, com1, com2):
        """Enables multiple tries to commit to db."""
        for i in range(5):  # retry up to 5 times
            try:
                self.cursor.execute("BEGIN;")
                self.cursor.execute(com1, com2)
                self.db.commit()
            except sqlite3.OperationalError as e:
                if "locked" in str(e).lower():
                    time.sleep(0.5)  # wait before retry
                else:
                    print("DataBase: commit failed.")

    # insert a sample in data base
    # IN: descp: description string
    # IN: point: parametric point (dict)
    # IN: ref: reference mesh if your sample has the same topology
    # IN: data: the tree, zone, data
    # IN: compressionTol: if None: lossless else compression accuracy
    def register(self, descp, point, ref=None, data=None, compressionTol=None):
        """Register data in db."""
        if self.mode == 'r': raise ValueError('register: can not write in read only mode.')

        if len(point) != len(self.parameters):
            raise ValueError("register: must have all parameters: "+str(self.parameters))
        if ref is None: ref = "None"

        varString = ''
        if data is not None:
            if isinstance(data, numpy.ndarray):
                data = ['data', data, [], 'DataArray_t']
                vars = ['data']
            elif isinstance(data, list):
                if Internal.typeOfNode(data) == -1:
                    raise ValueError("register: data is invalid.")
                vars = C.getVarNames(data, excludeXYZ=True)
                variables = vars[0]
                if ref != "None" and Internal.getNodeFromName(data, 'CoordinateX') is not None:
                    variables += ['dx', 'dy', 'dz']
            else:
                raise ValueError("register: data is invalid.")
            varString = ''
            for v in variables: varString += v+','
            if len(varString) > 2: varString = varString[:-1]

        # check if parameters already exists in db
        # and return id
        com = ''
        for p in point:
            com += "%s <= %g+1.e-10 AND %s >= %g-1.e-10 AND "%(p, point[p], p, point[p])
        if len(com) >= 4: com = com[:-4]
        com1 = 'SELECT * FROM %s WHERE '%self.name
        com1 = com1+com
        self.cursor.execute(com1)
        q = self.cursor.fetchall()
        if q != []: id = q[0][0]
        else:
            com = "SELECT id FROM %s ORDER BY id DESC LIMIT 1"%self.name
            self.cursor.execute(com)
            row = self.cursor.fetchone()
            if row is not None: id = row[0]+1
            else: id = 1 # this is an error

        # get current date
        now = datetime.datetime.now()
        dateString = now.strftime("%Y-%m-%dT%H:%M:%S")

        # register in sql
        com = f'REPLACE INTO {self.name}'
        com += '(id, descp, date, reference, variables'
        com2 = ' VALUES (?, ?, ?, ?, ?'
        com3 = [id, descp, dateString, ref, varString]
        for p in self.parameters:
            com += ', '+p
            com2 += ', ?'
            if p in point:
                com3.append(point[p])
            else: raise ValueError("register: parameter %s not found in input."%p)
        com += ')'
        com2 += ')'
        self.concurrentExecute(com+com2, com3)
        #self.cursor.execute(com+com2, com3)
        #self.db.commit()

        # register dcoords/fields cgns
        if data is not None:
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            tp, ntype = Internal.node2PyTree(data)
            tp = Internal.copyRef(tp)
            zones = Internal.getZones(tp)
            for z in zones:
                dcoords = None
                FC = Internal.getNodeFromType1(z, 'GridCoordinates_t')
                if FC is not None: px = Internal.getNodeFromName1(FC, 'CoordinateX')
                else: px = None
                if px is not None: # if coordinates in zone
                    # check reference if possible
                    refCgnsName = self.dirName+'/%s'%ref+'.cgns'
                    refFC = None
                    if os.path.exists(refCgnsName): # ok if ref=None
                        refFC = Filter.readNodesFromPaths(refCgnsName, ['%s/%s'%(Internal.getPath(tp,z),FC[0])])[0]
                        for i in Internal.getNodesFromType1(refFC, 'DataArray_t'):
                            Compressor._unpackNode(i)
                    if refFC is not None:
                        px = Internal.getNodeFromName1(FC, 'CoordinateX')
                        refx = Internal.getNodeFromName1(refFC, 'CoordinateX')
                        py = Internal.getNodeFromName1(FC, 'CoordinateY')
                        refy = Internal.getNodeFromName1(refFC, 'CoordinateY')
                        pz = Internal.getNodeFromName1(FC, 'CoordinateZ')
                        refz = Internal.getNodeFromName1(refFC, 'CoordinateZ')
                        dx = px[1] - refx[1]
                        dy = py[1] - refy[1]
                        dz = pz[1] - refz[1]
                        #retx = numpy.allclose(dx, 0., atol=1.e-10)
                        #rety = numpy.allclose(dy, 0., atol=1.e-10)
                        #retz = numpy.allclose(dz, 0., atol=1.e-10)
                        #if retx == False or rety == False or retz == False:
                        dxn = Internal.newDataArray('dx', dx)
                        dyn = Internal.newDataArray('dy', dy)
                        dzn = Internal.newDataArray('dz', dz)
                        dcoords = [dxn,dyn,dzn]

                # Fields
                ZT = Internal.getNodeFromType1(z, 'ZoneType_t')
                FS1 = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                FS2 = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)

                # if ref, clean zone and attach only fields and dcoords
                if ref != "None":
                    z[2] = []
                    if ZT is not None: z[2] += [ZT]
                    if FS1 is not None: z[2] += [FS1]
                    if FS2 is not None: z[2] += [FS2]

                    # attach dx in any case
                    if dcoords is not None:
                        if FS1 is None:
                            Internal.newFlowSolution(Internal.__FlowSolutionNodes__, 'Vertex', parent=z)
                            FS1 = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
                        # add dx,dy,dz to flow solution
                        FS1[2] += dcoords
            if compressionTol is None: Compressor._compressAll(tp) # lossless
            else:
                Compressor._compressFields(tp, tol=compressionTol, ctype=0) # approx
                Compressor._compressAll(tp) # lossless
            C.convertPyTree2File(tp, cgnsName)

    # get catalog
    def queryCatalog(self):
        """Return all rows of db."""
        self.cursor.execute('SELECT * FROM %s'%self.name)
        q = self.cursor.fetchall()
        return q

    # only for check, use self.columns instead
    def queryColumnNames(self):
        """Return column names of db."""
        self.cursor.execute("PRAGMA table_info(%s);"%self.name)
        q = self.cursor.fetchall()
        columnNames = [info[1] for info in q]
        return columnNames

    # return a query
    def query(self, com=None):
        """Return a query."""
        if com is None:
            # query catalog
            self.cursor.execute('SELECT * FROM %s'%self.name)
            rows = self.cursor.fetchall()
            return rows
        elif isinstance(com, str): # sql string command
            # query com
            com1 = 'SELECT * FROM %s WHERE '%self.name
            com1 = com1+com
            self.cursor.execute(com1)
            q = self.cursor.fetchall()
            return q
        elif isinstance(com, dict): # parametric point query
            com1 = ''
            for p in com:
                com1 += "%s <= %g+1.e-10 AND %s >= %g-1.e-10 AND "%(p, com[p], p, com[p])
            if len(com1) >= 4: com1 = com1[:-4]
            com2 = 'SELECT * FROM %s WHERE '%self.name
            com1 = com2+com1
            self.cursor.execute(com1)
            q = self.cursor.fetchall()
            return q

    # return True if parametric point exists in db
    def exist(self, point=None):
        if point is not None:
            q = self.query(point)
            if q == []: return False
            else: return True
        else: return False

    # load data from query
    # return list of trees
    def fetchTree(self, q):
        """Fetch a query and return tree."""
        ts = []; tref = None; refn = None
        for r in q:
            id = r[0]
            ref = r[3]
            if refn is None: refn = ref
            if ref != refn:
                ref = refn; tref = None
            # load ref mesh
            if ref != "None" and tref is None:
                cgnsName = self.dirName+'/%s'%ref+'.cgns'
                tref = C.convertFile2PyTree(cgnsName)
            # load data
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            t = C.convertFile2PyTree(cgnsName)
            if tref is not None:
                t = Internal.merge([tref, t])
                for z in Internal.getZones(t):
                    dx = Internal.getNodeFromName2(z, 'dx')
                    dy = Internal.getNodeFromName2(z, 'dy')
                    dz = Internal.getNodeFromName2(z, 'dz')
                    if dx is not None and dy is not None and dz is not None:
                        px = Internal.getNodeFromName2(z, 'CoordinateX')
                        px[1] += dx[1]
                        py = Internal.getNodeFromName2(z, 'CoordinateY')
                        py[1] += dy[1]
                        pz = Internal.getNodeFromName2(z, 'CoordinateZ')
                        pz[1] += dz[1]
            ts.append(t)
        return ts

    # fetch all parameters of q as a vector
    def fetchPointVector(self, q, exportJax=False):
        """Fetch a query and return param vector."""
        # build param points vector
        if exportJax: import jax.numpy as jnp
        nparam = len(q[0])-5
        nrows = len(q)
        if exportJax:
            param = jnp.zeros((nrows,nparam), dtype=numpy.float64)
            for c, r in enumerate(q): param.at[c,:] = r[5:]
        else:
            param = numpy.zeros((nrows,nparam), dtype=numpy.float64)
            for c, r in enumerate(q): param[c,:] = r[5:]
        return param

    # fetch all parameters of q as a point dict of DOE
    def fetchPoints(self, q):
        """Fetch a query and return points."""
        out = []
        for r in q:
            rp = r[5:]
            d = {}
            for c in range(len(rp)):
                d[self.parameters[c]] = rp[c]
            out.append(d)
        return out

    # fetch all parameters of q as a matrix
    def fetchMatrix(self, q, variables, exportJax=False):
        """Fetch a query and return a matrix."""
        if exportJax: import jax.numpy as jnp
        # sizes
        nq = len(q) # number of parametric points
        nv = len(variables) # number of variables
        pt0 = None
        pt1 = None
        matrix = None

        for c, r in enumerate(q): # columns: parametric points
            id = r[0]
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            h = Filter.Handle(cgnsName)
            a = h.loadSkeleton()
            h._loadVariables(a, var=variables)
            #vars = h.getVariables()
            zones = Internal.getZones(a)
            nz = len(zones)
            # Check that variables have been loaded
            for v in variables:
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                for z in zones:
                    p = Internal.getNodeFromName2(z, v)
                    if p is None:
                        raise ValueError('fetchMatrix: variable not found in data set.')

            if pt0 is None:
                sizetot = 0
                pt0 = numpy.zeros((nv,nz), dtype=numpy.int32)
                pt1 = numpy.zeros((nv,nz), dtype=numpy.int32)

                for x, v in enumerate(variables):
                    v = v.split(':')
                    if len(v) == 2: v = v[1]
                    else: v = v[0]
                    for n, z in enumerate(zones):
                        p = Internal.getNodeFromName(z, v)
                        nf = p[1].size
                        pt0[x,n] = sizetot
                        sizetot += nf
                        pt1[x,n] = sizetot

            for x, v in enumerate(variables): # rows : field variables
                v = v.split(':')
                if len(v) == 2: v = v[1]
                else: v = v[0]
                for n, z in enumerate(zones): # rows : field variables per zone
                    p = Internal.getNodeFromName(z, v)
                    nf = p[1].size
                    if exportJax:
                        if matrix is None:
                            matrix = jnp.zeros((sizetot, nq), dtype=numpy.float64)
                        jp = jnp.array(p.ravel('k'))
                        matrix.at[pt0[x,n]:pt1[x,n], c].set(jp)
                    else:
                        if matrix is None:
                            matrix = numpy.zeros((sizetot, nq), dtype=numpy.float64)
                        matrix[pt0[x,n]:pt1[x,n], c] = p[1].ravel('k')

        return matrix

    # fetch mass matrix of reference
    def fetchW(self, ref, exportJax=False):
        """Return the mass vector for ref."""
        import Generator.PyTree as G
        if exportJax: import jax.numpy as jnp
        cgnsName = self.dirName+'/%s'%ref+'.cgns'
        tref = C.convertFile2PyTree(cgnsName)
        G._getVolumeMap(tref)
        zones = Internal.getZones(tref)
        sizetot = 0
        for z in zones:
            p = Internal.getNodeFromName2(z, "vol")
            sizetot += p[1].size
        if exportJax:
            v = jnp.zeros((sizetot, 1), dtype=numpy.float64)
        else:
            v = numpy.zeros((sizetot, 1), dtype=numpy.float64)
        pt0 = 0
        for n, z in enumerate(zones):
            p = Internal.getNodeFromName2(z, "vol")
            size = p[1].size
            if exportJax: v.at[pt0:pt0+size, 0] = p[1].ravel('k')
            else: v[pt0:pt0+size, 0] = p[1].ravel('k')
            pt0 += size
        return v

    # delete rows corresponding to q
    def delete(self, q):
        """Delete queries from data base."""
        if self.mode == 'r': raise ValueError('delete: can not delete in read only mode.')

        for c, r in enumerate(q):
            # remove row in sql
            com1 = 'DELETE FROM %s WHERE '%self.name
            #com = 'id = %d AND descp = "%s" AND date = "%s" AND reference = "%s" AND variables = "%s"'%(r[0],r[1],r[2],r[3],r[4])
            #for c in range(5, len(r)):
            #    com += ' AND %s >= %g-1.e-10 AND %s <= %g+1.e-10'%(self.columns[c],r[c], self.columns[c],r[c])
            com = 'id = %d'%r[0]
            com1 = com1+com+';'
            self.cursor.execute(com1)
            self.db.commit()
            # remove file
            id = r[0]
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            try: os.remove(cgnsName)
            except: pass

    # join this data base with another one
    def join(self, db2):
        if self.mode == 'r': raise ValueError('join: can not join in read only mode.')
        # check that parameters are identicals
        for p in self.parameters:
            if p not in db2.parameters:
                raise ValueError('join: parameter %s is not in both dbs.'%p)

        # get the last id of self
        com = "SELECT id FROM %s ORDER BY id DESC LIMIT 1"%self.name
        self.cursor.execute(com)
        q = self.cursor.fetchone()
        if q is not None: id = q[0]+1
        else: id = 1

        # get all reference names in self and db2
        refNames1 = {}
        self.cursor.execute('SELECT * FROM %s'%self.name)
        rows = self.cursor.fetchall()
        for q in rows:
            ref = q[3]
            refNames1[ref] = 1
        refNames2 = {}
        db2.cursor.execute('SELECT * FROM %s'%db2.name)
        rows = db2.cursor.fetchall()
        for q in rows:
            ref = q[3]
            refNames2[ref] = 1

        # get all recordings in db2, add them and copy .cgns
        db2.cursor.execute('SELECT * FROM %s'%db2.name)
        rows = db2.cursor.fetchall()
        pts = db2.fetchPoints(rows)
        for c, pt in enumerate(pts):
            q = rows[c]
            id = q[0]
            print("join: adding", pt, "to "+self.dirName)
            self.register(q[1], pt, ref=q[3], data=None)
            id2 = self.cursor.lastrowid
            shutil.copy(db2.dirName+"/%05d"%id+".cgns", self.dirName+"/%05d"%id2+".cgns")

        # copy references if different names
        for r in refNames2:
            if r not in refNames1:
                shutil.copy(db2.dirName+"/"+r+".cgns", self.dirName)

    # monitor: write string to a log file monitor.txt
    def monitor(self, text):
        """Write text to log file in db."""
        logName = self.dirName+'/monitor.txt'
        now = datetime.datetime.now()
        dateString = now.strftime("%Y-%m-%dT%H:%M:%S")
        fp = open(logName, "a")
        fp.write(dateString+':'+text)
        fp.close()

    # print query to screen
    # mode=0: horizontal
    # mode=1: vertical
    def print(self, q, mode=0, fullListVar=False):
        """Print query to screen."""
        if mode == 0: # horizontal
            print(len(q),'entries.')
            size  = 20
            size2 = 15
            print('='*size*len(self.columns))

            txt = ''; txt2 = ''
            for c, i in enumerate(self.columns):
                if c != 3 and c != 4:
                    txt += str(i)[0:size2-1].ljust(size)
                else: txt2 += str(i)[0:size2-1].ljust(size)
            print(txt+txt2)

            print('='*size*len(self.columns))
            for r in q:
                txt = ''; txt2 = ''
                for c, i in enumerate(r):
                    txtAdd = str(i)[0:size2-1]
                    if len(str(i))>size2: txtAdd += '... '
                    if c != 3 and c != 4: txt += txtAdd.ljust(size)
                    else: txt2 += txtAdd.ljust(size)
                print(txt+txt2)
        else: # vertical
            print(len(q),'entries.')
            maxLength = 165
            size = max(maxLength//len(q),20)
            size2 = size-5
            print('='*((size+1)*len(q)+13))
            for c, i in enumerate(self.columns):
                txt = str(i)[0:12-1].ljust(12)+ '|'
                for r in q:
                    txtAdd = str(r[c])[0:size2-1]
                    if len(str(r[c]))>size2: txtAdd += '... '
                    txt += txtAdd.ljust(size) + '|'
                print(txt)
            print('='*((size+1)*len(q)+13))
            if fullListVar:
                print('\n \n')
                for i, r in enumerate(q):
                    print('=============================')
                    print('Full list of variables: id %d'%r[0])
                    print(r[4])
                    print('============================= \n')
