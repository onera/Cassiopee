"""Restricted pickle loader that safely deserializes objects by allowing only
a whitelist of approved modules and classes."""
import pickle
import importlib
import io

class RestrictedUnpickler(pickle.Unpickler):
    SAFE_GLOBALS = {
        "builtins": {"dict", "list", "set", "tuple", "str", "int", "float", "bool", "bytes"},
        "collections": {"OrderedDict"},
        "numpy": {"dtype"},
        "numpy.core.multiarray": {"_reconstruct", "scalar"},
        "numpy.core.numeric": {"_frombuffer"}
    }

    def __init__(self, file, *, encoding="ASCII", errors="strict"):
        super().__init__(file, encoding=encoding, errors=errors)

    def find_class(self, module, name):
        if module in self.SAFE_GLOBALS and name in self.SAFE_GLOBALS[module]:
            mod = importlib.import_module(module)
            return getattr(mod, name)
        raise pickle.UnpicklingError(f"Forbidden class: {module}.{name}")

def restrictedPickleLoad(file, encoding="ASCII", errors="strict"):
    return RestrictedUnpickler(file, encoding=encoding, errors=errors).load()

def restrictedPickleLoads(data):
    return RestrictedUnpickler(io.BytesIO(data)).load()
