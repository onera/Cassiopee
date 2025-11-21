# Manage base of data

# DataBase sumary:
# name <dir>
#   name.db (sql data base)
#   ref.cgns (common cgns data)
#   0001.cgns (key data)

import sqlite3
import os, numpy
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter
import Compressor.PyTree as Compressor

class DataBase:
    def __init__(self, name, parameters):
        if name[-3:] == ".db": name = name[-3:]
        # database name
        self.name = name
        # dir name
        self.dirName = name+'.db'
        # parameters name list
        self.parameters = parameters
        if not os.path.exists(self.dirName): os.mkdir(self.dirName)
        # pointer on sql db
        self.db = None
        self.cursor = None
        if os.path.exists("%s/%s.sql"%(self.dirName, self.name)):
            self.open()
        else:
            self.open()
            self.createTable()
        print("DataBase: creating %s"%self.dirName)

    # open sql data base
    def open(self):
        self.db = sqlite3.connect('%s/%s.sql'%(self.dirName,self.name))
        self.cursor = self.db.cursor()
        return self.cursor

    # close sql db
    def close(self):
        self.db.close()

    # create sql table
    def createTable(self):
        columns = [("id", "INTEGER PRIMARY KEY"),
                   ("descp", "TEXT"),
                   ("reference", "TEXT"),
                   ("variables", "TEXT")]

        for p in self.parameters:
            columns += [("%s"%p, "REAL")]

        columnsSql = ", ".join([f"{col} {dtype}" for col, dtype in columns])
        self.cursor.execute(f"CREATE TABLE IF NOT EXISTS {self.name} ({columnsSql})")

    # create reference (all without fields)
    def registerReference(self, t, name):
        cgnsName = self.dirName+'/%s'%name+'.cgns'
        tp = Internal.copyRef(t)
        Internal._rmNodesFromType(tp, 'FlowSolution_t')
        Compressor._compressAll(tp) # lossless
        C.convertPyTree2File(tp, cgnsName)
        return None

    # insert a sample in data base
    def register(self, id, descp, parameters, ref=None, variables=[], data=None):
        if len(parameters) != len(self.parameters):
            raise ValueError("register: must have all parameters: "+str(self.paramters))
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
            elif isinstance(data, list):
                if Internal.typeOfNode(data) == -1:
                    raise ValueError("register: data is invalid.")
                vars = C.getVarNames(data)
                if variables == []: variables = vars[0]
            else:
                raise ValueError("register: data is invalid.")

        # register in sql
        com = f'REPLACE INTO {self.name}'
        com += '(id, descp, reference, variables'
        com2 = ' VALUES (?, ?, ?, ?'
        com3 = [id, descp, ref, varString]
        for p in self.parameters:
            com += ', '+p
            com2 += ', ?'
            if p in parameters:
                com3.append(parameters[p])
            else: raise ValueError("register: parameter %s not found in input."%p)
        com += ')'
        com2 += ')'
        self.cursor.execute(com+com2, com3)
        self.db.commit()

        # register fields cgns
        if data is not None:
            cgnsName = self.dirName+'/%05d'%id+'.cgns'
            tp = Internal.copyRef(data)
            zones = Internal.getZones(tp)
            for z in zones:
                ZT = Internal.getNodeFromType1(z, 'ZoneType_t')
                FS = Internal.getNodesFromType1(z, 'FlowSolution_t')
                z[2] = [ZT]
                z[2] += FS
            Compressor._compressAll(tp) # lossless
            C.convertPyTree2File(tp, cgnsName)

    # get catalog
    def queryCatalog(self):
        self.cursor.execute('SELECT * FROM %s'%self.name)
        q = self.cursor.fetchall()
        return q

    # return a query
    def query(self, com=None):
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

    # load sample from query
    # return list of trees
    # mode=0: liste of trees
    # mode=1: list of FS
    # mode=2: one tree + multiple FS (not implemented)
    def fetchTree(self, q, mode=0):
        ts = []; tref = None; refn = None
        for r in q:
            id = r[0]
            ref = r[2]
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
            ts.append(t)
        return ts

    def fetchParam(self, q):
        # build param vector
        nparam = len(q[0])-4
        nrows = len(q)
        param = numpy.zeros((nrows,nparam), dtype=numpy.float64)
        for c, r in enumerate(q):
            param[c,:] = r[4:]
        return param

    def fetchMatrix(self, q, variable):
        # build param vector
        nparam = len(q[0])-4
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