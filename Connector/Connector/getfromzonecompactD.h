        E_Int loc = 1; //on recupere le noeud Flowsolution
# include "getfromzoneDcompact.h"
       int nvariables = PyList_Size(pyVariables);
       if (nvariables > 0)
       {
         PyObject* tpl0 = PyList_GetItem(pyVariables, 0);
         char* varname  = PyString_AsString(tpl0);
         strcpy(varStringOut,varname);
         for (int i = 1; i < nvariables; i++)
          {
            PyObject* tpl0 = PyList_GetItem(pyVariables, i);
            char* varname = PyString_AsString(tpl0);        
            strcat(varStringOut,","); strcat(varStringOut,varname);
          }
       }
       if( vartype <= 3 &&  vartype >= 1) nvars =5;
       else                               nvars =6;
