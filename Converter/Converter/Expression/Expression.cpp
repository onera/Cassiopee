#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

#define K_ARRAY_UNIQUE_SYMBOL
#include "Expression/ast.hpp"
#include "Expression/symbol_table.hpp"
#include "Expression/math_function.hpp"
#include "Memory/vector_view.hpp"
#include "converter.h"
using namespace K_FLD;

static struct python_parameters_dictionnary {
    std::size_t         nb_vars;
    std::vector<char *> m_param_names;
} py_params;

struct py_ast_handler {
    PyObject_HEAD Expression::ast *pt_ast;
};

PyAPI_DATA(PyTypeObject) py_ast_handler_type;

// Namespace anonyme, remplace avantageusement le static
namespace {
    std::vector<const char *> list_of_symbols() {
        auto &                    st = Expression::symbol_table::get();
        std::vector<const char *> symbols;
        symbols.reserve(st.size());
        for (auto &vk : st)
            symbols.push_back(vk.first.data());
        return symbols;
    }
    // .........................................................................................
    /*
    PyObject *py_list_of_symbols(PyObject *self) {
        auto      lst_symbols = list_of_symbols();
        PyObject *py_lst      = PyList_New(Py_ssize_t(lst_symbols.size()));
        if (py_lst == nullptr) {
            PyErr_SetString(PyExc_RuntimeError, "Failed to build list of symbols !");
            return NULL;
        }
        for (std::size_t index = 0; index < lst_symbols.size(); ++index) {
            PyList_SetItem(py_lst, Py_ssize_t(index), PyString_FromString(lst_symbols[index]));
        }
        return py_lst;
    }
    */
    // -----------------------------------------------------------------------------------------
    //py_ast_handler *_py_ast_handler_new() { return PyObject_NewVar(py_ast_handler, &py_ast_handler_type, 0); }
    // =========================================================================================
    void py_ast_handler_dealloc(py_ast_handler *self) {
        if (self->pt_ast != nullptr) delete self->pt_ast;
        self->ob_type->tp_free((PyObject *)self);
    }
    // -----------------------------------------------------------------------------------------
    PyObject *py_ast_handler_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
        py_ast_handler *self;
        self         = (py_ast_handler *)type->tp_alloc(type, 0);
        self->pt_ast = nullptr;
        return (PyObject *)self;
    }
    // -----------------------------------------------------------------------------------------
    PyObject *py_ast_handler_str(py_ast_handler *self) {
        return PyString_FromString(std::string(*self->pt_ast).data());
    }
    // -----------------------------------------------------------------------------------------
    struct data_array {
        PyObject * array;
        FldArrayF *f;
        FldArrayI *cn;
        E_Int      res;
    };
    PyObject *py_ast_handler_call(py_ast_handler *self, PyObject *args, PyObject *kwds) {
        std::unordered_map<std::string, vector_view<double>> cpp_dico;
        std::vector<double>                                  scal_variables;
        // Preparation de la table des symboles :
        // .......................................
        auto &st          = Expression::symbol_table::get();
        auto  symbol_kwds = list_of_symbols();

        // Recherche du nombre d'arguments passes pour l'expression avec verification
        // ..........................................................................
        if (not PyTuple_Check(args)) std::cerr << "args would be a python tuple. Strange..." << std::endl;
        assert((PyTuple_Check(args) != 0) && "args would be a python tuple. Strange... Call serial killer !");
        int nb_args_vars = PyTuple_Size(args);
        int nb_dict_vars = (kwds == nullptr ? 0 : PyDict_Size(kwds));
        int nb_vars      = std::max(nb_args_vars, 0) + std::max(nb_dict_vars, 0);
        if (nb_vars == 0) // Pas d'arguments de passe !
        {
            // Cela ne peut etre qu'une expression constante qui renvoie un double dans l'etat de l'arbre ast
            std::vector<double> result = (*self->pt_ast)(cpp_dico); // Le dictionnaire est vide...
            return PyFloat_FromDouble(result[0]);                   // Et on retourne un simple double
        }
        scal_variables.reserve(nb_vars);
        // =====================================================================
        // ==                   PARCOURTS DES ARGUMENTS                       ==
        // =====================================================================
        E_Int      ni, nj, nk;
        FldArrayF *f;
        FldArrayI *cn;
        char *     varString;
        char *     eltType;
        // On va prendre le premier argument de type array ( dans args ou dans kwds si args est vide ) comme modele
        // d'array Et on verifie sa validite en tant que tableau cassiopee Si la fonction ne contient que des scalaires
        // en parametre, nfld = 1 et npts = 1. res = 0
        E_Int     res, npts, nfld;
        PyObject *parray;
        if (nb_args_vars > 0) {
            parray = PyTuple_GetItem(args, 0);
            res    = K_ARRAY::getFromArray2(parray, varString, f, ni, nj, nk, cn, eltType);
            npts   = f->getSize();
            nfld   = f->getNfld();
        } else {
            assert((nb_dict_vars > 0) && "Logical error : if args is void, dicts would be not void !");
            PyObject *values = PyDict_Values(kwds);
            parray           = PyList_GetItem(values, 0);
            int ind          = 1;
            while (PyFloat_Check(parray) && (ind < nb_dict_vars)) {
                parray = PyList_GetItem(values, ind);
                ind += 1;
            }
            if (PyFloat_Check(parray)) {
                res  = 0;
                npts = 1;
                nfld = 1;
            } else {
                res                 = 1;
                PyObject * mem_view = PyMemoryView_FromObject(parray);
                Py_buffer *py_buf   = PyMemoryView_GET_BUFFER(mem_view);
                npts                = py_buf->len / sizeof(double);
                nfld                = 1;
                ni                  = npts;
                nj                  = 1;
                nk                  = 1;
            }
        }
        if ((res < 0) || (res > 2)) {
            std::string s_error = std::string(varString) +
                                  " is not valid. Wrong kind of array : " + std::to_string(res) + " doesn't exist !";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        if (nfld <= 0 || npts <= 0) {
            std::string s_error = std::string(varString) +
                                  " : No field defined in array => nfld = " + std::to_string(nfld) +
                                  " and npts = " + std::to_string(npts);
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            RELEASESHAREDB(res, parray, f, cn);
            return NULL;
        }
        if (nb_args_vars > 0) RELEASESHAREDB(res, parray, f, cn);
        // Parcourt des arguments
        std::vector<data_array> arrays(nb_args_vars);
        for (int iargs = 0; iargs < nb_args_vars; ++iargs) {
            // Tous les objets de la liste doivent etre des arrays de style
            // cassiopee compatible avec le tableau de reference :
            PyObject *array2 = PyTuple_GetItem(args, iargs);
            assert(array2 != nullptr);
            arrays[iargs].array = array2;
            E_Int ni2, nj2, nk2;
            char *varString2;
            char *eltType2;
            arrays[iargs].res =
                K_ARRAY::getFromArray2(array2, varString2, arrays[iargs].f, ni2, nj2, nk2, arrays[iargs].cn, eltType2);
            E_Int npts2 = arrays[iargs].f->getSize();
            // Verification de la coherence des tableaux avec le tableau resultat :
            if (arrays[iargs].res != res) {
                std::string s_error = std::string("Incompatible kind of array with result array : result array is a") +
                                      (res == 1 ? " structured array" : "n unstructured array") + " and " +
                                      std::string(varString2) + " is a" +
                                      (arrays[iargs].res == 1 ? " structured array" : "n unstructured array");
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 0; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                return NULL;
            }
            // Le nombre de champs peut etre different mais pas la taille de chaque champs !
            if (npts2 != npts) {
                std::string s_error = std::string("Incompatible dimension with result array :  result array has ") +
                                      std::to_string(npts) + " elements and " + std::string(varString2) + " has" +
                                      std::to_string(npts2) + " elements.";
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 0; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                return NULL;
            }
            // On extrait le nom de toutes les variables definies dans le
            // tableau
            std::vector<char *> vars;
            K_ARRAY::extractVars(varString2, vars);
            for (size_t ivar = 0; ivar < vars.size(); ++ivar) {
                std::string key_var = vars[ivar];
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_var) != symbol_kwds.end()) {
                    // La variable est bien utilisee dans l'expression, on
                    // rajoute le "tableau" au dictionnaire C++
                    cpp_dico[key_var] =
                        vector_view<double>((double *)arrays[iargs].f->begin(ivar + 1), arrays[iargs].f->getSize());
                    // Sinon on ignore simplement cette variable
                }
            }
            for (auto &v : vars)
                delete[] v;
        }
        if (nb_dict_vars > 0) {
            // Parcourt du dictionnaire,
            // Le dictionnaire ne doit contenir que des valeurs scalaires ou des objets acceptant un buffer memoire
            PyObject *py_keys = PyDict_Keys(kwds);
            for (int idict = 0; idict < nb_dict_vars; ++idict) {
                PyObject *  py_str = PyList_GetItem(py_keys, idict);
                std::string key_str(PyString_AsString(py_str));
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_str) == symbol_kwds.end()) {
                    std::string s_error = key_str + " is not a variable of the expression";
                    PyErr_SetString(PyExc_NameError, s_error.c_str());
                    for (int jargs = 0; jargs < nb_args_vars; jargs++)
                        RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                    return NULL;
                }
                PyObject *py_val = PyDict_GetItem(kwds, py_str);
                if (PyFloat_Check(py_val)) // C'est un double de passe
                {
                    double x = PyFloat_AsDouble(py_val);
                    scal_variables.push_back(x);
                    cpp_dico[key_str] = vector_view<double>(&scal_variables.back(), 1);
                } else {
                    // On vérifie si on a bien un objet de type tableau ( compatible avec un memoryview )
                    PyObject *mem_view = PyMemoryView_FromObject(py_val);
                    if (mem_view == nullptr) {
                        std::string s_error = key_str + " has wrong kind of value ( not a float or an array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        for (int jargs = 0; jargs < nb_args_vars; jargs++)
                            RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                        return NULL;
                    }
                    Py_buffer *py_buf = PyMemoryView_GET_BUFFER(mem_view);
                    // On vérifie qu'on a bien des "doubles"
                    if (py_buf->format[0] != 'd') {
                        std::string s_error =
                            key_str + " has wrong kind of value ( not double precision floats in array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        for (int jargs = 0; jargs < nb_args_vars; jargs++)
                            RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                        return NULL;
                    }
                    cpp_dico[key_str] = vector_view<double>((double *)py_buf->buf, py_buf->len / sizeof(double));
                }
            }
        }
        // std::cerr << "evaluating" << std::endl;
        PyObject *result;
        if (npts == 1) {
            std::vector<double> cpp_res = (*self->pt_ast)(cpp_dico);
            return PyFloat_FromDouble(cpp_res[0]);
        } else {
            const char *resultStr = "Result";
            if (res == 1) // Si structure :
            {
                //        std::cerr << "Je construis un tableau de dimension " << ni << "x" << nj << "x" << nk <<
                //        std::endl;
                result = K_ARRAY::buildArray(1, resultStr, ni, nj, nk);
            } else {
                // A VERIFIER AVEC CHRISTOPHE
                E_Int csize = arrays[0].cn->getSize();
                std::cerr << "Je construis un tableau de dimension " << csize << std::endl;
                result = K_ARRAY::buildArray(1, resultStr, npts, csize, -1, eltType, false, csize);
            }
            E_Float *fnp;
            if (res == 1) {
                fnp = K_ARRAY::getFieldPtr(result);
            } else {
                fnp = K_ARRAY::getFieldPtr(result);
                E_Int *cnp = K_ARRAY::getConnectPtr(result);
                E_Int* cnpp = K_ARRAY::getConnectPtr(arrays[0].array);
                E_Int  size = arrays[0].cn->getSize();
                E_Int  i;
#pragma omp parallel for shared(size, cnpp, cnp) private(i)
                for (i = 0; i < size; i++)
                    cnp[i] = cnpp[i];
            }
            auto vout = K_MEMORY::vector_view<E_Float>(fnp, npts);
            Py_BEGIN_ALLOW_THREADS;
            self->pt_ast->eval(cpp_dico, vout);
            Py_END_ALLOW_THREADS;
            for (int jargs = 0; jargs < nb_args_vars; jargs++) {
                RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
            }
            return result;
        }
        Py_RETURN_NONE;
    }
    // -------------------------------------------------------------------------
    const char *eval_doc =
        R"DOC(
    Evaluate the expression, storing the result in the variable outName stored in the output array.
    The first argument of the method is the array where stores the result and the second argument the name of the variable where the result will be stored.
    Example :
    --------
      >>> expr = Expression.ast("{x}**2+{y}**2+{z}**2")
      >>> param = ["x,y,z", [ numpy.array([1,2,3,4]), numpy.array([4,3,2,1]), numpy.array([5,6,7,8]) ] ]  # Verifier avec christophe sur la syntaxe de l'array....
      >>> result = ["norm1,norm2", [ numpy.array(4), numpy.array(4) ] ]
      >>> expr.eval(result, "norm1", param)
      >>> expr.eval(result, "norm2", x = [-1,-2,-3,-4], y = [ 2, 4, 6, 8], z = [1,3,5,7] )
    )DOC";
    PyObject *py_ast_handler_eval(py_ast_handler *self, PyObject *args, PyObject *kwds) {
        std::unordered_map<std::string, vector_view<double>> cpp_dico;
        std::vector<double>                                  scal_variables;
        // Preparation de la table des symboles :
        // .......................................
        auto &st          = Expression::symbol_table::get();
        auto  symbol_kwds = list_of_symbols();

        // Recherche du nombre d'arguments passes pour l'expression avec verification
        // ..........................................................................
        assert((PyTuple_Check(args) != 0) && "args would be a python list. Strange... Call serial killer !");
        assert(((kwds == nullptr) || (PyDict_Check(kwds) != 0)) &&
               "kwds would be a python dictionnary. Strange... Call serial killed !");
        int nb_args_vars = PyTuple_GET_SIZE(args);
        if (nb_args_vars == 0) {
            std::string s_error = "An array must be given as first argument of the eval method";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        int nb_dict_vars = (kwds == nullptr ? 0 : PyDict_Size(kwds));
        int nb_vars      = nb_args_vars + nb_dict_vars;
        scal_variables.reserve(nb_vars);
        // =====================================================================
        // ==                   PARCOURTS DES ARGUMENTS                       ==
        // =====================================================================
        E_Int      ni, nj, nk;
        FldArrayF *f;
        FldArrayI *cn;
        char *     varString;
        char *     eltType;
        // On va prendre le premier argument de args  comme modele
        // d'array Et on verifie sa validite en tant que tableau cassiopee.
        E_Int     res, npts, nfld;
        PyObject *parray = PyTuple_GetItem(args, 0);
        assert(parray != nullptr);
        res  = K_ARRAY::getFromArray2(parray, varString, f, ni, nj, nk, cn, eltType);
        npts = f->getSize();
        nfld = f->getNfld();
        if ((res < 1) || (res > 2)) {
            std::string s_error = std::string(varString) +
                                  " is not valid. Wrong kind of array : " + std::to_string(res) + " doesn't exist !";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        if (nfld <= 0 || npts <= 0) {
            std::string s_error = std::string(varString) +
                                  " : No field defined in array => nfld = " + std::to_string(nfld) +
                                  " and npts = " + std::to_string(npts);
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            RELEASESHAREDB(res, parray, f, cn);
            return NULL;
        }
        PyObject *py_name = PyTuple_GetItem(args, 1);
        char* outName;
        if (PyString_Check(py_name))
        {
         outName = PyString_AsString(py_name);   
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(py_name))
        {
            outName = PyBytes_AsString(PyUnicode_AsUTF8String(py_name));
        }
#endif
        else
        {
            std::string s_error = std::string("Second argument of the method must be the name of a variable");
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            RELEASESHAREDB(res, parray, f, cn);
            return NULL;
        }
        // Parcourt des arguments ( commence a deux car le zero est reserve pour le tableau de sortie et le un pour le
        // nom de la variable de sortie )
        bool must_release = true;
        std::vector<data_array> arrays(nb_args_vars);
        for (int iargs = 2; iargs < nb_args_vars; ++iargs) {
            // Tous les objets de la liste doivent etre des arrays de style
            // cassiopee compatible avec le tableau de reference :
            PyObject *array2    = PyTuple_GetItem(args, iargs);
            arrays[iargs].array = array2;
            E_Int ni2, nj2, nk2;
            char *varString2;
            char *eltType2;
            arrays[iargs].res =
                K_ARRAY::getFromArray2(array2, varString2, arrays[iargs].f, ni2, nj2, nk2, arrays[iargs].cn, eltType2);
            must_release = must_release and (arrays[iargs].f->begin() != f->begin());
            E_Int npts2 = arrays[iargs].f->getSize();
            // Verification de la coherence des tableaux avec le tableau resultat :
            if (arrays[iargs].res != res) {
                std::string s_error = std::string("Incompatible kind of array with result array : result array is a") +
                                      (res == 1 ? " structured array" : "n unstructured array") + " and " +
                                      std::string(varString2) + " is a" +
                                      (arrays[iargs].res == 1 ? " structured array" : "n unstructured array");
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 2; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                RELEASESHAREDB(res, parray, f, cn);
                return NULL;
            }
            // Le nombre de champs peut etre different mais pas la taille de chaque champs !
            if (npts2 != npts) {
                std::string s_error = std::string("Incompatible dimension with result array : result array has ") +
                                      std::to_string(npts) + " elements and " + std::string(varString2) + " has" +
                                      std::to_string(npts2) + " elements.";
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 2; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                RELEASESHAREDB(res, parray, f, cn);
                return NULL;
            }
            // On extrait le nom de toutes les variables definies dans le
            // tableau
            std::vector<char *> vars;
            K_ARRAY::extractVars(varString2, vars);
            for (size_t ivar = 0; ivar < vars.size(); ++ivar) {
                std::string key_var = vars[ivar];
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_var) != symbol_kwds.end()) {
                    // La variable est bien utilisee dans l'expression, on
                    // rajoute le "tableau" au dictionnaire C++
                    cpp_dico[key_var] =
                        vector_view<double>((double *)arrays[iargs].f->begin(ivar + 1), arrays[iargs].f->getSize());
                    // Sinon on ignore simplement cette variable
                }
            }
            for (auto &v : vars)
                delete[] v;
        }
        // Parcourt du dictionnaire,
        // Le dictionnaire ne doit contenir que des valeurs scalaires ou des objets acceptant un buffer memoire
        if (nb_dict_vars > 0) {
            PyObject *py_keys = PyDict_Keys(kwds);
            for (int idict = 0; idict < nb_dict_vars; ++idict) {
                PyObject *  py_str = PyList_GetItem(py_keys, idict);
                std::string key_str(PyString_AsString(py_str));
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_str) == symbol_kwds.end()) {
                    std::string s_error = key_str + " is not a variable of the expression";
                    PyErr_SetString(PyExc_NameError, s_error.c_str());
                    RELEASESHAREDB(res, parray, f, cn);
                    return NULL;
                }
                PyObject *py_val = PyDict_GetItem(kwds, py_str);
                if (PyFloat_Check(py_val)) // C'est un double de passe
                {
                    double x = PyFloat_AsDouble(py_val);
                    scal_variables.push_back(x);
                    cpp_dico[key_str] = vector_view<double>(&scal_variables.back(), 1);
                } else {
                    // On vérifie si on a bien un objet de type tableau ( compatible avec un memoryview )
                    PyObject *mem_view = PyMemoryView_FromObject(py_val);
                    if (mem_view == nullptr) {
                        std::string s_error = key_str + " has wrong kind of value ( not a float or array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        RELEASESHAREDB(res, parray, f, cn);
                        return NULL;
                    }
                    Py_buffer *py_buf = PyMemoryView_GET_BUFFER(mem_view);
                    // On vérifie qu'on a bien des "doubles"
                    if (py_buf->format[0] != 'd') {
                        std::string s_error =
                            key_str + " has wrong kind of value ( not double precision floats in array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        RELEASESHAREDB(res, parray, f, cn);
                        return NULL;
                    }
                    cpp_dico[key_str] = vector_view<double>((double *)py_buf->buf, py_buf->len / sizeof(double));
                }
            }
        }
        // On cherche tout d'abord si la variable demandee existe dans le array :
        E_Int    posvar = K_ARRAY::isNamePresent(outName, varString);
        E_Int res2;
        if (posvar == -1) // Ce nom n'existe pas, on rajoute le champs a l'array
        {
            // On va prendre le premier argument de args  comme modele
            // d'array Et on verifie sa validite en tant que tableau cassiopee.
            K_ARRAY::addFieldInArray(parray,outName);
            res2  = K_ARRAY::getFromArray2(parray, varString, f, ni, nj, nk, cn, eltType);
            posvar = K_ARRAY::isNamePresent(outName, varString);
        }
        assert(posvar >=0);
        E_Float* fnp = f->begin(posvar+1);
        assert(fnp != nullptr);
        auto vout = K_MEMORY::vector_view<E_Float>(fnp, npts);
        Py_BEGIN_ALLOW_THREADS;
        self->pt_ast->eval(cpp_dico, vout);
        Py_END_ALLOW_THREADS;

        for ( int jargs = 2; jargs < nb_args_vars; jargs++ ) {
            RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
        }

        if ( must_release ) {
            RELEASESHAREDB(res2, parray, f, cn);
        }

        Py_RETURN_NONE;
    }
    // -------------------------------------------------------------------------
    const char *run_doc =
        R"DOC(
    Evaluate the expression, without returning result. The expression must be have an assignment operator if one want get the result.

    Example :
    ---------
      >>> expr = Expression.ast("{norm} = {x}**2+{y}**2+{z}**2")
      >>> param = C.array("x,y,z", 4, 1, 1)
      >>> nrm   = C.array("norm", 4, 1, 1)
      >>> param2= C.array("x,y,z,norm", 4, 1, 1)
      >>> expr.run(param, nrm)
      >>> expr.run(param2)
    )DOC";
    PyObject *py_ast_handler_run(py_ast_handler *self, PyObject *args, PyObject *kwds) {
        std::unordered_map<std::string, vector_view<double>> cpp_dico;
        std::vector<double>                                  scal_variables;
        // Preparation de la table des symboles :
        // .......................................
        auto &st          = Expression::symbol_table::get();
        auto  symbol_kwds = list_of_symbols();

        // Recherche du nombre d'arguments passes pour l'expression avec verification
        // ..........................................................................
        if (args == nullptr)
        {
            std::string s_error = "None args ( NULL pointer ). Is it possible ?";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        assert((PyTuple_Check(args) != 0) && "args would be a python list. Strange... Call serial killer !");
        assert(((kwds == nullptr) || (PyDict_Check(kwds) != 0)) &&
               "kwds would be a python dictionnary. Strange... Call serial killed !");
        int nb_args_vars = PyTuple_GET_SIZE(args);
        if (nb_args_vars == 0) {
            std::string s_error = "An array must be given as first argument of the eval method";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        int nb_dict_vars = (kwds == nullptr ? 0 : PyDict_Size(kwds));
        int nb_vars      = nb_args_vars + nb_dict_vars;
        scal_variables.reserve(nb_vars);
        // =====================================================================
        // ==                   PARCOURTS DES ARGUMENTS                       ==
        // =====================================================================
        E_Int      ni, nj, nk;
        FldArrayF *f;
        FldArrayI *cn;
        char *     varString;
        char *     eltType;
        // On va prendre le premier argument de args  comme modele
        // d'array Et on verifie sa validite en tant que tableau cassiopee.
        E_Int     res, npts, nfld;

        PyObject *parray = PyTuple_GetItem(args, 0);
        assert(parray != nullptr);
        res  = K_ARRAY::getFromArray2(parray, varString, f, ni, nj, nk, cn, eltType);
        npts = f->getSize();
        nfld = f->getNfld();
        if ((res < 1) || (res > 2)) {
            std::string s_error = std::string(varString) +
                                  " is not valid. Wrong kind of array : " + std::to_string(res) + " doesn't exist !";
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            return NULL;
        }
        if (nfld <= 0 || npts <= 0) {
            std::string s_error = std::string(varString) +
                                  " : No field defined in array => nfld = " + std::to_string(nfld) +
                                  " and npts = " + std::to_string(npts);
            PyErr_SetString(PyExc_ValueError, s_error.c_str());
            RELEASESHAREDB(res, parray, f, cn);
            return NULL;
        }
        RELEASESHAREDB(res, parray, f, cn);
        std::vector<data_array> arrays(nb_args_vars);
        for (int iargs = 0; iargs < nb_args_vars; ++iargs) {
            // Tous les objets de la liste doivent etre des arrays de style
            // cassiopee compatible avec le tableau de reference :
            PyObject *array2    = PyTuple_GetItem(args, iargs);
            arrays[iargs].array = array2;
            E_Int ni2, nj2, nk2;
            char *varString2;
            char *eltType2;
            arrays[iargs].res =
                K_ARRAY::getFromArray2(array2, varString2, arrays[iargs].f, ni2, nj2, nk2, arrays[iargs].cn, eltType2);
            E_Int npts2 = arrays[iargs].f->getSize();
            // Verification de la coherence des tableaux avec le tableau resultat :
            if (arrays[iargs].res != res) {
                std::string s_error = std::string("Incompatible kind of array with result array : result array is a") +
                                      (res == 1 ? " structured array" : "n unstructured array") + " and " +
                                      std::string(varString2) + " is a" +
                                      (arrays[iargs].res == 1 ? " structured array" : "n unstructured array");
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 2; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                RELEASESHAREDB(res, parray, f, cn);
                return NULL;
            }
            // Le nombre de champs peut etre different mais pas la taille de chaque champs !
            if (npts2 != npts) {
                std::string s_error = std::string("Incompatible dimension with result array : result array has ") +
                                      std::to_string(npts) + " elements and " + std::string(varString2) + " has" +
                                      std::to_string(npts2) + " elements.";
                PyErr_SetString(PyExc_ValueError, s_error.c_str());
                for (int jargs = 2; jargs <= iargs; jargs++)
                    RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
                RELEASESHAREDB(res, parray, f, cn);
                return NULL;
            }
            // On extrait le nom de toutes les variables definies dans le
            // tableau
            std::vector<char *> vars;
            K_ARRAY::extractVars(varString2, vars);
            for (size_t ivar = 0; ivar < vars.size(); ++ivar) {
                std::string key_var = vars[ivar];
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_var) != symbol_kwds.end()) {
                    // La variable est bien utilisee dans l'expression, on
                    // rajoute le "tableau" au dictionnaire C++
                    cpp_dico[key_var] =
                        vector_view<double>((double *)arrays[iargs].f->begin(ivar + 1), arrays[iargs].f->getSize());
                    // Sinon on ignore simplement cette variable
                }
            }
            for (auto &v : vars)
                delete[] v;
        }
        // Parcourt du dictionnaire,
        // Le dictionnaire ne doit contenir que des valeurs scalaires ou des objets acceptant un buffer memoire
        if (nb_dict_vars > 0) {
            PyObject *py_keys = PyDict_Keys(kwds);
            for (int idict = 0; idict < nb_dict_vars; ++idict) {
                PyObject *  py_str = PyList_GetItem(py_keys, idict);
                std::string key_str(PyString_AsString(py_str));
                if (std::find(symbol_kwds.begin(), symbol_kwds.end(), key_str) == symbol_kwds.end()) {
                    std::string s_error = key_str + " is not a variable of the expression";
                    PyErr_SetString(PyExc_NameError, s_error.c_str());
                    RELEASESHAREDB(res, parray, f, cn);
                    return NULL;
                }
                PyObject *py_val = PyDict_GetItem(kwds, py_str);
                if (PyFloat_Check(py_val)) // C'est un double de passe
                {
                    double x = PyFloat_AsDouble(py_val);
                    scal_variables.push_back(x);
                    cpp_dico[key_str] = vector_view<double>(&scal_variables.back(), 1);
                } else {
                    // On vérifie si on a bien un objet de type tableau ( compatible avec un memoryview )
                    PyObject *mem_view = PyMemoryView_FromObject(py_val);
                    if (mem_view == nullptr) {
                        std::string s_error = key_str + " has wrong kind of value ( not a float or array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        RELEASESHAREDB(res, parray, f, cn);
                        return NULL;
                    }
                    Py_buffer *py_buf = PyMemoryView_GET_BUFFER(mem_view);
                    // On vérifie qu'on a bien des "doubles"
                    if (py_buf->format[0] != 'd') {
                        std::string s_error =
                            key_str + " has wrong kind of value ( not double precision floats in array )";
                        PyErr_SetString(PyExc_NameError, s_error.c_str());
                        RELEASESHAREDB(res, parray, f, cn);
                        return NULL;
                    }
                    cpp_dico[key_str] = vector_view<double>((double *)py_buf->buf, py_buf->len / sizeof(double));
                }
            }
        }
        // On cherche tout d'abord si la variable demandee existe dans le array :
        Py_BEGIN_ALLOW_THREADS;
        self->pt_ast->eval(cpp_dico);
        Py_END_ALLOW_THREADS;

        for ( int jargs = 1; jargs < nb_args_vars; jargs++ ) {
            RELEASESHAREDB(arrays[jargs].res, arrays[jargs].array, arrays[jargs].f, arrays[jargs].cn);
        }
        
        Py_RETURN_NONE;
    }
    // ====================================================================================================>
    static PyMethodDef ast_methods[] = {
        {"eval", (PyCFunction)py_ast_handler_eval, METH_VARARGS | METH_KEYWORDS, eval_doc},  // Sentinelle
        {"run", (PyCFunction)py_ast_handler_run, METH_VARARGS | METH_KEYWORDS, run_doc}, {NULL} // Sentinelle
    };
    //
    static int ast_init(py_ast_handler *self, PyObject *args) {
        const char *expression;
        if (!PyArg_ParseTuple(args, "s", &expression)) return -1;
        try {
            self->pt_ast = new Expression::ast(expression);
        } catch(std::exception& e)
        {
            std::string s_error = std::string(expression) + " : " + e.what();
            PyErr_SetString(PyExc_SyntaxError, s_error.c_str());
            return -1;
        }
        catch(...)
        {
            std::string s_error = std::string(expression) + " : unexcepted exception !";
            PyErr_SetString(PyExc_RuntimeError, s_error.c_str());
            return -1;
        }
        return 0;
    }

} // namespace
// ========================================================================
// Definition of new type ast for python language
static const char ast_doc[]           = "Abstract syntax tree documentation to do !";
#if defined (_WIN32) || defined (_WIN64)
__declspec(dllexport)
#endif
PyTypeObject      py_ast_handler_type = {
    PyObject_HEAD_INIT(NULL) 0,               // Object size (minus type size)
    "Expression.ast",                         // Name of the python type
    sizeof(py_ast_handler),                   // Basic object size
    0,                                        // Item size
    (destructor)py_ast_handler_dealloc,       // Destructor
    (printfunc)0,                             // Print function
    0,                                        // get attr.
    0,                                        // set attr.
    0,                                        // Compare function
    (reprfunc)0,                              // Representation function
    0,                                        // Object as number
    0,                                        // Object as a sequence
    0,                                        // Object as a dictionnary
    0,                                        // Object as a hash table
    (ternaryfunc)py_ast_handler_call,         // Object as a function
    (reprfunc)py_ast_handler_str,             // Object as a string
    0,                                        // Get attr. Obj.
    0,                                        // Set attr. Obj.
    0,                                        // Object as a buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // Object caracteristics
    ast_doc,                                  // Documentation of the type
    0,                                        // How traverse this object for garbage collextor
    0,                                        // How clear the object (for garbage collector if needed)
    0,                                        // Rich compare function
    0,                                        // Weak list of set
    0,                                        // Iterator on the object
    0,                                        // Iterator next function
    ast_methods,                              // Methods for this object
    0,                                        // Members of this object
    0,                                        // getters and setters
    0,                                        // Base object
    0,                                        // Object dictionnary
    0,                                        // Descr. get (?)
    0,                                        // Descr. set (?)
    0,                                        // Dictionnary offset
    (initproc)ast_init,                       // Initialization method for instance of this type
    0,                                        // Object static allocation
    py_ast_handler_new,                       // Instance creation function
};
// ===========================================================================================================
namespace {
    const char *derivate_doc =
        R"DOC(
        Derivate the expression inside an expression tree.

        This function derivate an expression with all variables defined in this expression.
        New variables are created for derivate of the initial variables, pre-appending a "d_"
        at the left of each initial variables.

        Example :
        ========
            >>> import Generator as G
            >>> import Converter as C
            >>> import Converter.expression as expr
            >>> import numpy as np
            >>> crds = G.cart((0, 0, 0), (1, 1, 1), (3, 3, 1), api=1)
            >>> C._addVars(crds, 'norm')
            >>> C._addVars(crds, 'd_norm')
            >>> a = expr.ast("{norm} = {x}**2+{y}**2+{z}**2")
            >>> da = expr.derivate(a)
            >>> shp = crds[1][0].shape
            >>> da.run(crds, d_x=np.ones(shp), d_y=np.zeros(shp), d_z=np.zeros(shp))
            >>> print da
            (d_norm=(((2.000000*d_x)*x)+(((2.000000*d_y)*y)+((2.000000*d_z)*z))))
            >>> print crds
            ['x,y,z,norm,d_norm', array([[ 0.,  1.,  2.,  0.,  1.,  2.,  0.,  1.,  2.],
                                         [ 0.,  0.,  0.,  1.,  1.,  1.,  2.,  2.,  2.],
                                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                                         [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                                         [ 0.,  2.,  4.,  0.,  2.,  4.,  0.,  2.,  4.]]), 3, 3, 1]
    )DOC";
    PyObject *py_derivate(PyObject *self, PyObject *args) {
        py_ast_handler* py_ast;
        if ( not PyArg_ParseTuple(args, "O!", &py_ast_handler_type, &py_ast) )
        {
            return NULL;
        }
        auto da = py_ast->pt_ast->derivate();
        py_ast_handler* py_dast = (py_ast_handler*)py_ast_handler_new(&py_ast_handler_type, nullptr, nullptr);
        py_dast->pt_ast = new Expression::ast(da);
        return (PyObject*)py_dast;
    }
}
// ===========================================================================================================
static PyMethodDef expression_methods[] = {
    {"derivate", (PyCFunction)py_derivate, METH_VARARGS, derivate_doc},
    {NULL, NULL}
};
static char        expression_mod_doc[] = R"DOC(
This module is intented to evaluate and derivate some expressions given by the user as string litteral. The expression can be evaluated on vectors.
In this case, the expression is evaluated coefficient per coefficient to avoid some intermediate vector to copy which slow the global evaluation.

The choosen syntax is close to Tecplot and Cassiopee syntax ( look at _initVars method in Converter module ) :

- Each variable is enclosed with curly braces. By example : {x}, {hitcher guide to the galaxy}, 
)DOC";

PyMODINIT_FUNC initexpression() {
    PyObject *m;

    if (PyType_Ready(&py_ast_handler_type) < 0) return;
    m = Py_InitModule4("expression", expression_methods, expression_mod_doc, (PyObject *)NULL, PYTHON_API_VERSION);
    /* Très important : initialise numpy afin de pouvoir l'utiliser ici !!!!
     */
    import_array();

    Expression::init_math_functions();

    if (m == NULL) return;

    /** XXXX Add constants here */

    /** */
    Py_INCREF(&py_ast_handler_type);
    PyModule_AddObject(m, "ast", (PyObject *)&py_ast_handler_type);
}
