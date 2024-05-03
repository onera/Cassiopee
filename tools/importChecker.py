#!/usr/bin/env python
##############################################################################
#
# Copyright (c) 2003 Zope Corporation and Contributors.
# All Rights Reserved.
#
# This software is subject to the provisions of the Zope Public License,
# Version 2.1 (ZPL).  A copy of the ZPL should accompany this distribution.
# THIS SOFTWARE IS PROVIDED "AS IS" AND ANY AND ALL EXPRESS OR IMPLIED
# WARRANTIES ARE DISCLAIMED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF TITLE, MERCHANTABILITY, AGAINST INFRINGEMENT, AND FITNESS
# FOR A PARTICULAR PURPOSE.
#
##############################################################################
"""Import checker

This utility finds unused imports in Python modules.  Its output is
grep-like and thus emacs-friendly.

$Id: importchecker.py 123992 2012-01-09 11:13:03Z janwijbrand $
"""

try: import compiler
except: import ast as compiler
import os, os.path
import sys

def _findDottedNamesHelper(node, result):
    more_node = node
    name = node.__class__.__name__
    if name == 'Getattr':
        dotted = []
        while name == 'Getattr':
            dotted.append(node.attrname)
            node = node.expr
            name = node.__class__.__name__
        if name == 'Name':
            dotted.append(node.name)
            dotted.reverse()
            for i in range(1, len(dotted)):
                result.append('.'.join(dotted[:i]))
            result.append('.'.join(dotted))
            return
    elif name == 'Name':
        result.append(node.name)
        return
    elif name == 'AssAttr':
        # Can be on an import as well.
        # for instance
        # from x import y
        # y.k = v
        expr = node.expr
        result.append(getattr(expr, 'name', ''))
        return
    for child in more_node.getChildNodes():
        _findDottedNamesHelper(child, result)


def findDottedNames(node):
    """Find dotted names in an AST tree node
    """
    result = []
    _findDottedNamesHelper(node, result)
    return result


class ImportFinder:
    """An instance of this class will be used to walk over a compiler AST
    tree (a module). During that operation, the appropriate methods of
    this visitor will be called
    """

    def __init__(self):
        self._map = {}

    def visitFrom(self, stmt):
        """Will be called for 'from foo import bar' statements
        """
        # XXX take first two items of statement list as in Python 2.5 this
        # list contains more information items.
        module_name, names = stmt.asList()[:2]
        if module_name == '__future__':
            # we don't care what's imported from the future
            return
        names_dict = {}
        for orig_name, as_name in names:
            # we don't care about from import *
            if orig_name == '*':
                continue
            if as_name is None:
                name = orig_name
            else:
                name = as_name
            names_dict[name] = orig_name
        self._map.setdefault(module_name, {'names': names_dict,
                                           'lineno': stmt.lineno})

    def visitImport(self, stmt):
        """Will be called for 'import foo.bar' statements
        """
        for orig_name, as_name in stmt.names:
            if as_name is None:
                name = orig_name
            else:
                name = as_name
            self._map.setdefault(orig_name, {'names': {name: orig_name},
                                             'lineno': stmt.lineno})

    def getMap(self):
        return self._map


def findImports(mod):
    """Find import statements in module and put the result in a mapping.
    """
    visitor = ImportFinder()
    compiler.walk(mod, visitor)
    return visitor.getMap()


class Module:
    """This represents a python module.
    """

    def __init__(self, path):
        mod = compiler.parseFile(path)
        self._path = path
        self._map = findImports(mod)
        dottednames = {}
        self._dottednames = findDottedNames(mod)

    def getPath(self):
        """Return the path to this module's file.
        """
        return self._path

    def getImportedModuleNames(self):
        """Return the names of imported modules.
        """
        return self._map.keys()

    def getImportNames(self):
        """Return the names of imports; add dottednames as well.
        """
        result = []
        map = self._map
        for module_name in map.keys():
            for usedname, originalname in map[module_name]['names'].items():
                result.append((originalname, module_name))
                # add any other name that we could be using
                for dottedname in self._dottednames:
                    usednamedot = usedname + '.'
                    if dottedname.startswith(usednamedot):
                        attrname = dottedname[len(usednamedot):].split('.')[0]
                        result.append((attrname, module_name))
        return result

    def getUnusedImports(self):
        """Get unused imports of this module (the whole import info).
        """
        result = []
        for value in self._map.values():
            for usedname, originalname in value['names'].items():
                if usedname not in self._dottednames:
                    result.append((originalname, value['lineno']))
        return result


class ModuleFinder:

    def __init__(self):
        self._files = []

    def visit(self, arg, dirname, names):
        """This method will be called when we walk the filesystem
        tree. It looks for python modules and stored their filenames.
        """
        for name in names:
            # get all .py files that aren't weirdo emacs droppings
            if name.endswith('.py') and not name.startswith('.#'):
                self._files.append(os.path.join(dirname, name))

    def getModuleFilenames(self):
        return self._files


def findModules(path):
    """Find python modules in the given path and return their absolute
    filenames in a sequence.
    """
    finder = ModuleFinder()
    os.walk(path, finder.visit, ())
    return finder.getModuleFilenames()


class ImportDatabase:
    """This database keeps tracks of imports.

    It allows to NOT report cases where a module imports something
    just so that another module can import it (import dependencies).
    """

    def __init__(self, root_path):
        self._root_path = root_path
        self._modules = {}
        self._names = {}

    def resolveDottedModuleName(self, dotted_name, module):
        """Return path to file representing module, or None if no such
        thing. Can do this relative from module.
        """
        dotted_path = dotted_name.replace('.', '/')
        # try relative import first
        path = os.path.join(os.path.dirname(module.getPath()), dotted_path)
        path = self._resolveHelper(path)
        if path is not None:
            return path
        # absolute import (assumed to be from this tree)
        if os.path.isfile(os.path.join(self._root_path, '__init__.py')):
            startpath, dummy = os.path.split(self._root_path)
        else:
            startpath = self._root_path
        return self._resolveHelper(os.path.join(startpath, dotted_path))

    def _resolveHelper(self, path):
        if os.path.isfile(path + '.py'):
            return path + '.py'
        if os.path.isdir(path):
            path = os.path.join(path, '__init__.py')
            if os.path.isfile(path):
                return path
        return None

    def findModules(self):
        """Find modules in the given path.
        """
        for modulepath in findModules(self._root_path):
            module = Module(modulepath)
            self.addModule(module)

    def addModule(self, module):
        """Add information about a module to the database. A module in
        this case is not a python module object, but an instance of
        the above defined Module class.w
        """
        self_path = module.getPath()
        # do nothing if we already know about it
        if self._modules.has_key(self_path):
            return

        self._modules[self_path] = module

        # add imported names to internal names mapping; this will
        # allow us identify dependent imports later
        names = self._names
        for name, from_module_name in module.getImportNames():
            path = self.resolveDottedModuleName(from_module_name, module)
            t = (path, name)
            modulepaths = names.get(t, {})
            if not modulepaths.has_key(self_path):
                modulepaths[self_path] = 1
            names[t] = modulepaths

    def getUnusedImports(self):
        """Get unused imports of all known modules.
        """
        result = {}
        for path, module in self._modules.items():
            result[path] = self.getUnusedImportsInModule(module)
        return result

    def getUnusedImportsInModule(self, module):
        """Get all unused imports in a module.
        """
        result = []
        for name, lineno in module.getUnusedImports():
            if not self.isNameImportedFrom(name, module):
                result.append((name, lineno))
        return result

    def isNameImportedFrom(self, name, module):
        """Return true if name is imported from module by another module.
        """
        return self._names.has_key((module.getPath(), name))

    def getModulesImportingNameFrom(self, name, module):
        """Return list of known modules that import name from module.
        """
        result = []
        for path in self._names.get((module.getPath(), name), {}).keys():
            result.append(self._modules[path])
        return result


def main(path=None):
    cwd = os.getcwd()
    lencwd = len(cwd)+1

    try:
        path = path or sys.argv[1]
    except IndexError:
        print("No path supplied")
        sys.exit(1)

    fullpath = os.path.abspath(path)
    path = fullpath
    if not os.path.isdir(fullpath):
        path = os.path.dirname(fullpath)

    db = ImportDatabase(path)
    if os.path.isdir(fullpath):
        db.findModules()
    else:
        db.addModule(Module(fullpath))
    unused_imports = db.getUnusedImports()
    module_paths = unused_imports.keys()
    module_paths.sort()
    for path in module_paths:
        info = unused_imports[path]
        if path.startswith(cwd):
            path = path[lencwd:]
        if not info:
            continue
        line2names = {}
        for name, line in info:
            names = line2names.get(line, [])
            names.append(name)
            line2names[line] = names
        lines = line2names.keys()
        lines.sort()
        for line in lines:
            names = ', '.join(line2names[line])
            print("%s :%s: %s" % (path, line, names))

if __name__ == "__main__":
    main()
