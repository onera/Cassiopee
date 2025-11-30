# Manage base of CFD data

# DataBase sumary:
# name <dir>
#   name.db (sql data base)
#   ref.cgns (common cgns data)
#   0001.cgns (key associated data)

# vocabulary:
# query(sqlstring) -> return a query (list of full sql data)
# fetch(query) -> return real data corresponding to query
# point: a dict of parameter names/value (as in DOE)

import sqlite3
import os, numpy, time, datetime
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter
import Compressor.PyTree as Compressor

class DataBase:
    def __init__(self, name, parameters=None):
        if name[-3:] == ".db": name = name[:-3]
        # database name
        self.name = name
        # directory name
        self.dirName = name+'.db'
        if not os.path.exists(self.dirName): os.mkdir(self.dirName)
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

    # open sql data base
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
    def registerReference(self, t, name):
        """Write reference file used for coordinates and connectivity."""
        cgnsName = self.dirName+'/%s'%name+'.cgns'
        tp = Internal.copyRef(t)
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
    def register(self, descp, point, ref=None, variables=[], data=None, tol=-1.):
        """Register data in db."""
        if len(point) != len(self.parameters):
            raise ValueError("register: must have all parameters: "+str(self.parameters))
        if ref is None: ref = "ref1"

        if isinstance(variables, list):
            varString = ''
            for v in variables: varString += ','+v
        elif isinstance(variables, str):
            varString = variables
        else: raise ValueError("register: variables must be a string or a list of strings.")

        if data is not None:
            if isinstance(data, numpy.ndarray):
                data = ['data', data, [], 'DataArray_t']
                vars = ['data']
            elif isinstance(data, list):
                if Internal.typeOfNode(data) == -1:
                    raise ValueError("register: data is invalid.")
                vars = C.getVarNames(data, excludeXYZ=True)
                if variables == []: variables = vars[0]
            else:
                raise ValueError("register: data is invalid.")
            varString = ''
            for v in variables: varString += v+','
            if len(varString) > 2: varString = varString[:-1]

        # check if parameters already exists in db
        # and return id
        com = ''
        for p in point:
            com += "%s = %g AND "%(p, point[p])
        if len(com) >= 4: com = com[:-4]
        com1 = 'SELECT * FROM %s WHERE '%self.name
        com1 = com1+com
        self.cursor.execute(com1)
        q = self.cursor.fetchall()
        if q != []: id = q[0][0]
        else: id = self.cursor.lastrowid+1

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

        # register fields cgns
        if data is not None:
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            tp = Internal.copyRef(data)
            zones = Internal.getZones(tp)
            dcoords = None
            for z in zones:
                # Coordinates
                FC = Internal.getNodeFromType1(z, 'GridCoordinates_t')
                if FC is not None:
                    # check reference if possible
                    refCgnsName = self.dirName+'/%s'%ref+'.cgns'
                    refFC = None
                    if os.path.exists(refCgnsName): # of if ref=None
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
                        retx = numpy.allclose(dx, 0., atol=1.e-10)
                        rety = numpy.allclose(dy, 0., atol=1.e-10)
                        retz = numpy.allclose(dz, 0., atol=1.e-10)
                        if retx == False or rety == False or retz == False:
                            dxn = Internal.newDataArray('dx', dx)
                            dyn = Internal.newDataArray('dy', dy)
                            dzn = Internal.newDataArray('dz', dz)
                            dcoords = [dxn,dyn,dzn]
                # Fields
                ZT = Internal.getNodeFromType1(z, 'ZoneType_t')
                FS = Internal.getNodesFromType1(z, 'FlowSolution_t')
                z[2] = []
                if ZT is not None: z[2] += [ZT]
                if FS is not None: z[2] += FS
                if dcoords is not None:
                    # add dx,dy,dz to flow solution
                    b = Internal.getNodeFromType1(z, 'FlowSolution_t')
                    if b is not None: b[2] += dcoords
            if tol <= 0: Compressor._compressAll(tp) # lossless
            else:
                Compressor._compressFields(tp, tol=tol, ctype=0) # approx
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
        # query com
        com1 = 'SELECT * FROM %s WHERE '%self.name
        com1 = com1+com
        self.cursor.execute(com1)
        q = self.cursor.fetchall()
        return q

    # return True if parameters exist in db
    def exist(self, point=None):
        com = ''
        if point is not None:
            for p in point:
                com += "%s = %g AND "%(p, point[p])
        if len(com) >= 4: com = com[:-4]
        com1 = 'SELECT * FROM %s WHERE '%self.name
        com1 = com1+com
        self.cursor.execute(com1)
        q = self.cursor.fetchall()
        if q == []: return False
        else: return True

    # load sample from query
    # return list of trees
    # mode=0: liste of trees
    # mode=1: list of FS
    # mode=2: one tree + multiple FS (not implemented)
    def fetchTree(self, q, mode=0):
        """Fetch a query and return tree."""
        ts = []; tref = None; refn = None
        for r in q:
            id = r[0]
            ref = r[3]
            if refn is None: refn = ref
            if ref != refn:
                ref = refn; tref = None
            # load ref mesh
            if tref is None and mode == 0:
                cgnsName = self.dirName+'/%s'%ref+'.cgns'
                tref = C.convertFile2PyTree(cgnsName)
            # load data
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            t = C.convertFile2PyTree(cgnsName)
            if mode == 0:
                t = Internal.merge([tref, t])
                for z in Internal.getZones(t):
                    dx = Internal.getNodeFromName2(z, 'dx')
                    dy = Internal.getNodeFromName2(z, 'dy')
                    dz = Internal.getNodeFromName2(z, 'dz')
                    if dx is not None and dy is not None and dz is not None:
                        px = Internal.getNodeFromName2(z, 'CoordinateX')
                        px = px + dx
                        py = Internal.getNodeFromName2(z, 'CoordinateX')
                        py = py + dy
                        pz = Internal.getNodeFromName2(z, 'CoordinateX')
                        pz = pz + dz
            ts.append(t)
        return ts

    # fetch all parameters as a vector
    def fetchParam(self, q):
        """Fetch a query and return param vector."""
        # build param vector
        nparam = len(q[0])-5
        nrows = len(q)
        param = numpy.zeros((nrows,nparam), dtype=numpy.float64)
        for c, r in enumerate(q):
            param[c,:] = r[5:]
        return param

    # fetch all parameters of q as a point dict of DOE
    def fetchPoint(self, q):
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
    def fetchMatrix(self, q, variable):
        """Fetch a query and return matrix."""
        # build param vector
        nrows = len(q)
        matrix = None

        for c, r in enumerate(q):
            id = r[0]
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            h = Filter.Handle(cgnsName)
            a = h.loadSkeleton()
            #h._loadZonesWoVars(a)
            h._loadVariables(a, var=[variable])
            zones = Internal.getZones(a)
            z = zones[0]
            p = Internal.getNodeFromName(z, variable)
            nf = p[1].size
            if matrix is None:
                matrix = numpy.zeros((nrows,nf), dtype=numpy.float64)
            matrix[c,:] = p[1].ravel('k')
        return matrix

    # delete rows corresponding to q
    def delete(self, q):
        """Delete queries from data base."""
        for c, r in enumerate(q):
            # remove row in sql
            com1 = 'DELETE FROM %s WHERE '%self.name
            com = 'id = %d AND descp = "%s" AND date = "%s" AND reference = "%s" AND variables = "%s"'%(r[0],r[1],r[2],r[3],r[4])
            for c in range(5, len(r)):
                com += ' AND %s = %g'%(self.columns[c],r[c])
            com1 = com1+com+';'
            self.cursor.execute(com1)
            # remove file
            id = r[0]
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            os.remove(cgnsName)

    # monitor: write string to a log file monitor.txt
    def monitor(self, text):
        """Write text to log file in db."""
        logName = self.dirName+'/monitor.txt'
        now = datetime.datetime.now()
        dateString = now.strftime("%Y-%m-%dT%H:%M:%S")
        fp = open(logName, "a")
        fp.write(dateString+':'+text)
        fp.close()