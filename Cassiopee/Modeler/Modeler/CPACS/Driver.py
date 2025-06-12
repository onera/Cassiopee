# Driver for cpacs
# IN: param. cpacs
import sympy
import Modeler.modeler

# parameter: attached to a cpacs element or global
class DesignParameter:
    """Design parameter."""
    def __init__(self, name, expr="0.", range=["0.", "1."], path=None):
        # parameter name
        self.name = name
        # sympy symbol
        self.symbol = sympy.symbols(name)
        # parameter path in cpacs (if any)
        self.path = path
        # range
        self.rangeMin = range[0]
        self.rangeMax = range[1]
        # expression (string)
        self.expr = expr
        # sympy equation
        self.eq = sympy.Eq(self.symbol, eval(expr))
        # instantiated value
        self.value = 0.

    def set(self, value):

        return

class Driver:
    """Driver"""
    def __init__(self, cpacsFile):
        # associated cpacs file
        self.cpacsFile = cpacsFile
        xmlFile = self.cpacsFile.replace(".cpacs", ".xml")
        self.xmlFile = xmlFile
        # all design parameters
        self.designParameters = {}
        # control parameters
        self.controlParameters = {}
        # geometric parameters
        self.geomParameters = {}

    def register(self, param):
        """Register a param."""
        self.designParameters[param.name] = param

    def findParametersInCPACS(self):
        return

    def instantiateParameters(self):
        """Eval parameters."""
        eqs = []; vars = []
        for k in self.designParameters:
            eqs.append(self.designParameters[k].eq)
        for k in self.designParameters:
            vars.append(self.designParameters[k].symbol)
        solution = sympy.solve(eqs, vars)
        print(solution)
        for s in solution: print(s.name)

        # check if all parameters are instantiated
        inst = 0
        for s in solution:
            if solution[s].is_Float: # instantiated
                self.designParameters[s.name].value = float(solution[s])
                print("Parameter %s set to %f."%(s.name, self.designParameters[s.name].value))
                inst += 1
        if inst != len(self.designParameters.keys()):
            print("Warning: some parameters were not instantiated.")
            print("Warning: please instantiate %d more.", inst - len(self.designParameters.keys()))
        return

    def checkValuesInRange(self):
        for k in self.designParameters:
            p = self.designParameters[k]
            value = p.value
            rangeMin = p.rangeMin
            rangeMax = p.rangeMax
            if value < rangeMin: print("Warning: value out of range.")
            if value < rangeMax: print("Warning: value out of range.")

    def instantiateCPACS(self):
        f = open(self.cpacsFile, "r")
        g = open(self.xmlFile, "w")
        for line in f:
            line = self.replace(line)
            print(line)
            g.write(line)
        f.close()
        g.close()
        return

    def replace(self, line):
        for k in self.designParameters:
            p = self.designParameters[k]
            value = str(p.value)
            line = line.replace("#"+p.name, value)
        return line

    def instantiateCAD(self, fileName):
        Modeler.modeler.exportCAD(self.xmlFile, fileName, "fmt_step")
        return

    def instantiateMesh(self):
        return
