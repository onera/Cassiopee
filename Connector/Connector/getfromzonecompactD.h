        E_Int loc = 1; //on recupere le noeud Flowsolution
# include "getfromzoneDcompact.h"
       int nvariables = PyList_Size(pyVariables);
       if (nvariables > 0)
       {
         PyObject* tpl0 = PyList_GetItem(pyVariables, 0);
         char* varname  = NULL;
         if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
         else if (PyUnicode_Check(tpl0)) varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0)); 
#endif
         strcpy(varStringOut, varname);
         for (int i = 1; i < nvariables; i++)
          {
            PyObject* tpl0 = PyList_GetItem(pyVariables, i);
            char* varname = NULL;
            if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
            else if (PyUnicode_Check(tpl0)) varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0)); 
#endif
            strcat(varStringOut,","); strcat(varStringOut,varname);
          }
       }
       if( vartype <= 3 &&  vartype >= 1) nvars =5;
       else                               nvars =6;
