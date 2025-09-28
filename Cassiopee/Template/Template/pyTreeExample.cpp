/* Exemple de fonction utilisant l'interface pyTree */

#include "template.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Recupere une zone (les pointeurs sont partages avec python)
//=============================================================================
PyObject* K_TEMPLATE::pyTreeExample(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone;
  if (!PYPARSETUPLE_(args, O_ SSS_, &zone, &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  E_Int ni, nj, nk, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;
  E_Int res = K_PYTREE::getFromZone(zone, 1, 2, varString,
                                    fields, locs, ni, nj, nk,
                                    cn, cnSize, cnNfld,
                                    eltType, hook,
                                    GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  /* Plus d'info dans KCore/PyTree/PyTree.h */
  if (res == 1)
  {
    printf("Zone structuree (%d %d %d).\n", ni, nj, nk);
  }
  else if (res == 2)
  {
    printf("Zone non structuree de type %s.\n", eltType);
  }
  else
  {
    printf("Zone invalide.\n");
    Py_INCREF(Py_None);
    return Py_None;
  }

  /* Exemple d'acces aux champs */
  if (res == 1) /* structure */
  {
    E_Int nfld = fields.size();
    printf("La zone contient %d champs.\n", nfld);
    printf("qui sont %s.\n", varString);
    // Pointeur sur les donnees ordonnees i+j*ni+k*ni*nj pour chaque champ
    // Ce pointeur est partage avec python.
    // Si on modifie les valeurs du tableau ici, elles seront
    // modifiees dans python
    for (E_Int n = 0; n < nfld; n++)
    {
      E_Float* p = fields[n];
      if (locs[n] == 0) // champs en noeuds
      {
        printf("Champs no %d en noeuds.\n", n);
        for (E_Int k = 0; k < nk; k++)
          for (E_Int j = 0; j < nj; j++)
            for (E_Int i = 0; i < ni; i++)
            {
              E_Int ind = i+j*ni+k*ni*nj;
              printf("%d %d %d -> %f\n", i,j,k, p[ind]);
            }
      }
      else if (locs[n] == 1) // champ en centres
      {
        printf("Champs no %d en centres.\n", n);
        for (E_Int k = 0; k < nk-1; k++)
          for (E_Int j = 0; j < nj-1; j++)
            for (E_Int i = 0; i < ni-1; i++)
            {
              E_Int ind = i+j*(ni-1)+k*(ni-1)*(nj-1);
              printf("%d %d %d -> %f\n", i,j,k, p[ind]);
            }
      }
    }
  }
  else if (res == 2) /* no structure */
  {

  }

  RELEASESHAREDZ(hook, varString, eltType);

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Recupere un arbre (les pointeurs sont partages avec python)
//=============================================================================
PyObject* K_TEMPLATE::pyTreeExample1(PyObject* self, PyObject* args)
{
  PyObject* t;
  if (!PYPARSETUPLE_(args, O_, &t)) return NULL;

  /* Exemple de parcours d'arbre en C */
  vector<PyArrayObject*> hook;

  /* Recupere les noeuds base */
  vector<PyObject*> bases;
  K_PYTREE::getNodesFromType1(t, "CGNSBase_t", bases);

  for (size_t i = 0; i < bases.size(); i++)
  {
    char* baseName = K_PYTREE::getNodeName(bases[i]);
    printf("Detected bases: %s\n", baseName);
  }

  /* Recupere les noeuds zones pour chaque base */
  for (size_t i = 0; i < bases.size(); i++)
  {
    vector<PyObject*> zones;
    PyObject* base = bases[i];
    K_PYTREE::getNodesFromType1(base, "Zone_t", zones);
    for (size_t j = 0; i < zones.size(); j++)
    {
      char* zoneName = K_PYTREE::getNodeName(zones[j];
      printf("Detected zones: %s\n", zoneName);

      // Recuperation du noeud GridCoordinates
      PyObject* node = K_PYTREE::getNodeFromName1(zones[j], "GridCoordinates");
      PyObject* x = K_PYTREE::getNodeFromName1(node, "CoordinateX");
      E_Float* ptx = K_PYTREE::getValueAF(x, hook);
    }
  }

  RELEASEHOOK(hook);
  Py_INCREF(Py_None);
  return Py_None;
}
