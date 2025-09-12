/* Exemple de passage d'un Cassiopee array de python au C */

#include "template.h"

using namespace K_FLD;

//=============================================================================
// Recupere un Cassioppee array et le modifie
//=============================================================================
PyObject* K_TEMPLATE::arrayExample(PyObject* self, PyObject* args)
{
  // Recupere le pointeur sur l'objet python
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Recupere les pointeurs sur les tableaux contenus dans l'array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                               cn, eltType);

  /* More info on arrays in KCore/Array/Array.h */
  
  if (res == 1) // C'est un array structure
  {
    E_Int nfld = f->getNfld();
    printf("C'est un array structure (%d %d %d).\n", ni, nj, nk);
    printf("Il contient %d champs.\n", nfld);
    printf("qui sont %s.\n", varString);
    // Pointeur sur les donnees ordonnees i+j*ni+k*ni*nj pour chaque champ
    // Ce pointeur est partage avec python.
    // Si on modifie les valeurs du tableau ici, elles seront
    // modifiees dans python
    E_Float* pt = f->begin();

    // Exemple pour ecrire le contenu d'un champ
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* p = f->begin(n);
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            E_Int ind = i+j*ni+k*ni*nj;
            printf("%d %d %d -> %f\n", i,j,k, p[ind]);
          }
    }
  }
  else if (res == 2) // C'est un array non structure
  {
    E_Int nfld = f->getNfld();
    printf("C'est un array non structure de type %s.\n", eltType);
    printf("Il contient %d champs.\n", nfld);
    printf("qui sont %s.\n", varString);

    if (strcmp(eltType, "NGON") == 0 || strcmp(eltType, "NGON*") == 0)
    {
      // NGON
    }
    else
    {
      // Elements basiques (TRI, TETRA, ...)
      E_Int npts = f->getSize();
      E_Int nelts = cn->getSize();
      printf("Il contient %d noeuds et %d elements.\n", npts, nelts);
      // Pointeur sur la connectivite elts->noeuds
      // Recuperer le nbre de ptrs?
      E_Int* cn1 = cn->begin(1);
      
      
    }
  }
  else // Ce n'est pas un array
  {
    printf("Ce n'est pas un array.\n");
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  // Release memory
  RELEASESHAREDB(res, array, f, cn);

  Py_INCREF(Py_None);
  return Py_None;
}
