/*
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

# include "transform.h"
using namespace K_FLD;
using namespace K_SEARCH;

//=============================================================================
/* Déclarations et Fonctions */
//=============================================================================

inline bool isSamePoint(E_Float p1x, E_Float p1y, E_Float p1z, E_Float p2x, E_Float p2y, E_Float p2z) 
{
  return (std::abs(p1x - p2x) < 1e-10 && std::abs(p1y - p2y) < 1e-10 && std::abs(p1z - p2z) < 1e-10);
}

// ============================================================================
/* Conservative smoothing (Kuprat)*/
// IN: array  : maillage (structuré i-array ou non-structuré 1D dans le plan (x,y))
// IN: sweeps : nombre de passes de lissage
// IN: retour : 0 = lissage aller seulement, 1 = lissage aller + retour (default à 0, donne un lissage assymétrique) 
// IN: pas    : incrément de parcours (1, 2 ou 3) mis à 1 par default
// OUT: tableau lissé (même connectivité, coordonnées modifiées, type conservé)
// ============================================================================
PyObject* K_TRANSFORM::consSmooth(PyObject* self, PyObject* args)
{
  PyObject* array;

  E_Int sweeps = 1;
  E_Int pas = 2;
  E_Int retour = 1;
  
  if (!PYPARSETUPLE_(args, O_ III_, &array, &sweeps,  &retour, &pas))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =  K_ARRAY::getFromArray3(array, varString,
                                      f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: array is invalid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (retour != 0 && retour != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: int retour is invalid.");
    return NULL;
  }

  if (pas < 1 && pas > 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "consSmooth: int step is invalid.");
    return NULL;
  }

  // Build output
  PyObject* tpl;
  FldArrayF* fo;
  FldArrayI* cno = NULL;

  if (res == 1) // structure
  {
    tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, f->getApi());
    K_ARRAY::getFromArray3(tpl, fo);

    E_Float* x = fo->begin(posx);
    E_Float* y = fo->begin(posy);
    E_Float* z = fo->begin(posz);

    E_Int npts = im;
    
    E_Int idx0;
    E_Int idx1;
    E_Int idx2;
    E_Int idx3;

    // Ouvert ou fermé ?
    E_Float dx = x[1]-x[0], dy = y[1]-y[0], dz = z[1]-z[0];
    E_Float refL = std::sqrt(dx*dx + dy*dy + dz*dz); // longueur du 1er segment

    E_Float distBouts = std::sqrt((x[0]-x[npts-1])*(x[0]-x[npts-1]) + (y[0]-y[npts-1])*(y[0]-y[npts-1]) + (z[0]-z[npts-1])*(z[0]-z[npts-1]));

    E_Int ouvert = (distBouts < 1e-6 * refL) ? 0 : 1;
    E_Int nUnique = (ouvert == 0) ? npts - 1 : npts; 
    E_Int start = 0;
    E_Int end   = (ouvert == 0) ? nUnique : nUnique - 3; 
    
    for (E_Int k = 0; k < sweeps; k++)
    {

      for (E_Int j = 0; j < retour + 1; j++)
      {
        
        for (E_Int i = start; i < end; i = i + pas)
        {
          
          if (ouvert == 0)

          {
            // Cyclique, toujours dans [0, nUnique-1]
            idx0 = (j == 0) ? i % nUnique : (nUnique - i % nUnique) % nUnique;
            idx1 = (j == 0) ? (i + 1) % nUnique : (nUnique - ((i + 1)   % nUnique)) % nUnique;
            idx2 = (j == 0) ? (i + 2) % nUnique : (nUnique - ((i + 2)   % nUnique)) % nUnique;
            idx3 = (j == 0) ? (i + 3) % nUnique : (nUnique - ((i + 3)   % nUnique)) % nUnique;
          }
          else 
          {
            idx0 = (j == 0) ? i  : (nUnique-1) - i;
            idx1 = (j == 0) ? (i + 1) : (nUnique-1) - (i + 1);
            idx2 = (j == 0) ? (i + 2)  : (nUnique-1) - (i + 2);
            idx3 = (j == 0) ? (i + 3)  : (nUnique-1) - (i + 3) ;
          }

          /* On récupère les coordonnées des points i, i+1, i+2, i+3 */
          E_Float xi = x[idx0]; E_Float yi = y[idx0]; E_Float zi = z[idx0]; 
          E_Float xip1 = x[idx1]; E_Float yip1 = y[idx1]; E_Float zip1 = z[idx1]; 
          E_Float xip2 = x[idx2]; E_Float yip2 = y[idx2]; E_Float zip2 = z[idx2]; 
          E_Float xip3 = x[idx3]; E_Float yip3 = y[idx3]; E_Float zip3 = z[idx3];

          /* On calcule les différences de points qui nous interessent (i+3 - i), (i+2 - i) et (i+1 - i) */
          E_Float dv1x = xip1 - xi; E_Float dv1y = yip1 - yi; E_Float dv1z = zip1 - zi; /* xi+1 - xi */
          E_Float dv2x = xip2 - xi; E_Float dv2y = yip2 - yi; E_Float dv2z = zip2 - zi; /* xi+2 - xi */
          E_Float dv3x = xip3 - xi; E_Float dv3y = yip3 - yi; E_Float dv3z = zip3 - zi; /* xi+3 - xi */

          /* On définit uNormal = unit normal to baseline (i+3;i) */
          E_Float normeDv3 = (dv3x * dv3x) + (dv3y * dv3y) + (dv3z * dv3z); /*  ||{xi+3-xi}**||² = ||xi+3-xi||² */
          
          if (normeDv3 < 1e-12) continue;
          E_Float divNorme = 1.0 / normeDv3 ;/*  1 / ||xi+3-xi||² */

          E_Float uNormalx = divNorme * (-dv3y); E_Float uNormaly = divNorme * (dv3x); E_Float uNormalz = 0;  /* {xi+3-xi}** / ||xi+3-xi||² */

          /* Calcul de l'aire signée */
          E_Float aire = 0.5 * (dv3x * dv2y - dv3y * dv2x) + \
                        0.5 * (dv2x * dv1y - dv2y * dv1x) ;

          E_Float h = 1.5 * aire;

          x[idx1] = 2.0/3.0 * xi + 1.0 /3.0 * xip3 + h * uNormalx; y[idx1] = 2.0/3.0 * yi + 1.0 /3.0 * yip3 + h * uNormaly; z[idx1] = 2.0/3.0 * zi + 1.0 /3.0 * zip3 + h * uNormalz;
          x[idx2] = 1.0/3.0 * xi + 2.0 /3.0 * xip3 + h * uNormalx; y[idx2] = 1.0/3.0 * yi + 2.0 /3.0 * yip3 + h * uNormaly; z[idx2] = 1.0/3.0 * zi + 2.0 /3.0 * zip3 + h * uNormalz;
          
        }
      }
    
    }

    if (ouvert == 0)
    {
      x[nUnique] = x[0]; 
      y[nUnique] = y[0]; 
      z[nUnique] = z[0];
    }

    RELEASESHAREDS(tpl, fo); // Libération Structurée
  
  }
  else if (res == 2) // non structuré
  {
     
    tpl = K_ARRAY::buildArray3(f->getNfld(), varString,
                               f->getSize(), cn->getSize(),
                               eltType, false, f->getApi());
    K_ARRAY::getFromArray3(tpl, fo, cno);

    E_Float* x = fo->begin(posx);
    E_Float* y = fo->begin(posy);
    E_Float* z = fo->begin(posz);

    for (E_Int v = 1; v <= f->getNfld(); v++) 
    {
      E_Float* fp1 = f->begin(v);
      E_Float* fp2 = fo->begin(v);
      for (E_Int i = 0; i < f->getSize(); i++) fp2[i] = fp1[i];
    }
    
    for (E_Int n = 1; n <= cno->getNfld(); n++) 
    {
      for (E_Int i = 0; i < cno->getSize(); i++) 
      {
        (*cno)(i,n) = (*cn)(i,n);
      }
    }

    E_Int npts = f->getSize(); 

    std::vector< std::vector<E_Int> > nodeAdj(npts);// vertex/elt adjacents
    K_CONNECT::connectEV2VE(*cn, nodeAdj);
  
    // Ouvert ou fermé ?? determine les noeuds frontières
    E_Int n1 = -1; E_Int n2 = -1;
    E_Int ouvert = 0;

    for(E_Int p = 0; p < npts; p++) 
    {
      if((E_Int)nodeAdj[p].size() == 1) 
      {
        if (n1 == -1) n1 = p;
        else n2 = p;
      } 
    }

    if (n1 != -1 && n2 != -1)
    {
      if (isSamePoint(x[n1], y[n1], z[n1], x[n2], y[n2], z[n2])) 
      {
        //CAS 2 : Maillage topologiquement OUVERT, mais geometriquement FERME. On soude alors virtuellement 
        ouvert = 0;
        printf("Géométrie fermé avec point doublons : %d et %d\n", n1, n2);
      }
      else  // CAS 3 : Vraie courbe OUVERTE (Les extremites vont rester fixes)
      {
        ouvert = 1;
        printf("Géométrie ouverte : frontières %d et %d\n", n1, n2);
      }
    }

    // On crée une chaine ordonnée de noeuds
    std::vector<E_Int> bar;
    
    E_Int firstNode = (ouvert == 1) ? n1
                    : (n1 != -1) ? n1 // CAS 2 : fermé avec doublon, on part du premier doublon
                    : ((*cn)(0, 1) - 1); //CAS 1 : fermé pur, on part de n'importe quel noeud 
    E_Int cur = firstNode, prev = -1;
    bar.push_back(cur);

    while (true)
    {
      E_Int next = -1;
      
      for (E_Int elt : nodeAdj[cur])  // elt = indice d'élément
      {
        E_Int nA = (*cn)(elt, 1) - 1;
        E_Int nB = (*cn)(elt, 2) - 1;
        E_Int neighbor = (nA == cur) ? nB : nA; // le voisin = l'autre nœud du segment
        if (neighbor != prev) { next = neighbor; break; }
      }

      if (next == -1) break;                              // CAS 3 : nœud avec 1 seul voisin atteint = fin d'une chaîne ouverte
      if (ouvert == 0 && next == firstNode) break;        // CAS 1 : retour au départ
      if (ouvert == 0 && n1 != -1 && next == n2) break;   // CAS 2 : on ne garde pas le doublon
      
      bar.push_back(next);
      prev = cur;
      cur  = next;
      if (ouvert == 1 && cur == n2) break; // on a atteint l'autre extrémité ouverte
    }

    E_Int nUnique = (E_Int)bar.size();

    E_Int start = 0;
    E_Int end   = (ouvert == 0) ? nUnique : nUnique - 3;
  
    for (E_Int k = 0; k < sweeps; k++)
    {
      for (E_Int j = 0; j < retour + 1; j++)
      {
        for (E_Int i = start; i < end; i += pas)
        {
          E_Int c0, c1, c2, c3; // indices dans le bar
          if (ouvert == 0)
          {
            // CAS fermé
            c0 = (j == 0) ? i % nUnique : (nUnique - i % nUnique) % nUnique;
            c1 = (j == 0) ? (i + 1) % nUnique : (nUnique - ((i + 1)   % nUnique)) % nUnique;
            c2 = (j == 0) ? (i + 2) % nUnique : (nUnique - ((i + 2)   % nUnique)) % nUnique;
            c3 = (j == 0) ? (i + 3) % nUnique : (nUnique - ((i + 3)   % nUnique)) % nUnique;
          }
          else
          {
            // CAS ouvert : 
            c0 = (j == 0) ? i : (nUnique - 1) - i;
            c1 = (j == 0) ? (i + 1) : (nUnique - 1) - (i + 1);
            c2 = (j == 0) ? (i + 2) : (nUnique - 1) - (i + 2);
            c3 = (j == 0) ? (i + 3) : (nUnique - 1) - (i + 3);
          }

          // indices réels des nœuds
          E_Int idx0 = bar[c0];
          E_Int idx1 = bar[c1];
          E_Int idx2 = bar[c2];
          E_Int idx3 = bar[c3];

          E_Float xi = x[idx0]; E_Float yi = y[idx0]; E_Float zi = z[idx0]; 
          E_Float xip1 = x[idx1]; E_Float yip1 = y[idx1]; E_Float zip1 = z[idx1]; 
          E_Float xip2 = x[idx2]; E_Float yip2 = y[idx2]; E_Float zip2 = z[idx2]; 
          E_Float xip3 = x[idx3]; E_Float yip3 = y[idx3]; E_Float zip3 = z[idx3];

          /* On calcule les différences de points qui nous interessent (i+3 - i), (i+2 - i) et (i+1 - i) */
          E_Float dv1x = xip1 - xi; E_Float dv1y = yip1 - yi; E_Float dv1z = zip1 - zi; /* xi+1 - xi */
          E_Float dv2x = xip2 - xi; E_Float dv2y = yip2 - yi; E_Float dv2z = zip2 - zi; /* xi+2 - xi */
          E_Float dv3x = xip3 - xi; E_Float dv3y = yip3 - yi; E_Float dv3z = zip3 - zi; /* xi+3 - xi */

          /* On définit uNormal = unit normal to baseline (i+3;i) */
          E_Float normeDv3 = (dv3x * dv3x) + (dv3y * dv3y) + (dv3z * dv3z); /*  ||{xi+3-xi}**||² = ||xi+3-xi||² */
          
          if (normeDv3 < 1e-12) continue;
          E_Float divNorme = 1.0 / normeDv3 ;/*  1 / ||xi+3-xi||² */

          E_Float uNormalx = divNorme * (-dv3y); E_Float uNormaly = divNorme * (dv3x); E_Float uNormalz = 0;  /* {xi+3-xi}** / ||xi+3-xi||² */

          /* Calcul de l'aire signée */
          E_Float aire = 0.5 * (dv3x * dv2y - dv3y * dv2x) + \
                        0.5 * (dv2x * dv1y - dv2y * dv1x) ;

          E_Float h = 1.5 * aire;

          x[idx1] = 2.0/3.0 * xi + 1.0 /3.0 * xip3 + h * uNormalx; y[idx1] = 2.0/3.0 * yi + 1.0 /3.0 * yip3 + h * uNormaly; z[idx1] = 2.0/3.0 * zi + 1.0 /3.0 * zip3 + h * uNormalz;
          x[idx2] = 1.0/3.0 * xi + 2.0 /3.0 * xip3 + h * uNormalx; y[idx2] = 1.0/3.0 * yi + 2.0 /3.0 * yip3 + h * uNormaly; z[idx2] = 1.0/3.0 * zi + 2.0 /3.0 * zip3 + h * uNormalz;
          
          

        }
      }
    }

    if (ouvert == 0 && n1 != -1 && n2 != -1)
    {
      x[n2] = x[n1];
      y[n2] = y[n1];
      z[n2] = z[n1];
    }

    RELEASESHAREDU(tpl, fo, cno); // Libération Non-Structurée
  }

  RELEASESHAREDB(res, array, f, cn); 
  return tpl;

}



