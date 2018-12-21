/*    
    Copyright 2013-2019 Onera.

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
# include "converter.h"
# include "Search/BbTree.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Conformisation topologique d'un NGON
   Le jeu de faces correspondant a une plus grande face la remplace 
   IN: f: champ NGON
   IN: posx, posy, posy: position x,y,z dans f
   IN: cn: connectivite NGON
   IN: tol: tolerance
   OUT: cno: nouvelle connectivite NGON */
//=============================================================================
void K_CONVERTER::conformizeNGon(
  FldArrayF& f, 
  E_Int posx, E_Int posy, E_Int posz,
  FldArrayI& cn, E_Float tol, FldArrayI*& cno)
{
  // On enregistre les BB des faces dans un BbTree
  E_Float* px = f.begin(posx);
  E_Float* py = f.begin(posy);
  E_Float* pz = f.begin(posz);

  // Calcul de posFaces
  E_Int* cnp = cn.begin();
  E_Int nfaces = cnp[0];
  FldArrayI posFaces;
  K_CONNECT::getPosFaces(cn, posFaces);
 
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  K_SEARCH::BbTree3D* bbtree;
  vector<BBox3DType*> boxes(nfaces);

#pragma omp parallel default(shared)
  {
    E_Int* ptr; E_Int np; E_Int ind;
    E_Float minB[3]; E_Float maxB[3];
    E_Float x, y, z;
    for (E_Int i = 0; i < nfaces; i++)
    {
      ptr = cnp+posFaces[i];
      np = ptr[0];
      minB[0] = K_CONST::E_MAX_FLOAT; minB[1] = K_CONST::E_MAX_FLOAT; minB[2] = K_CONST::E_MAX_FLOAT;
      maxB[0] = -K_CONST::E_MAX_FLOAT; maxB[1] = -K_CONST::E_MAX_FLOAT; maxB[2] = -K_CONST::E_MAX_FLOAT;
      for (E_Int j = 0; j < np; j++)
      {
        ind = ptr[j+1]-1;
        x = px[ind]; y = py[ind]; z = pz[ind];
        minB[0] = std::min(minB[0], x);
        minB[1] = std::min(minB[1], y);
        minB[2] = std::min(minB[2], z);
        maxB[0] = std::max(maxB[0], x);
        maxB[1] = std::max(maxB[1], y);
        maxB[2] = std::max(maxB[2], z);
      }
      boxes[i] = new BBox3DType(minB, maxB);
    }
  }
  
  bbtree = new K_SEARCH::BbTree3D(boxes, 1.e-12);
  //printf("building the BB Tree\n");

  // Recherche des matching faces (topogeometrique)
  
  E_Int** indir = new E_Int* [nfaces];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < nfaces; i++) indir[i] = NULL;

  //printf("input mesh has %d faces.\n", nfaces);

#pragma omp parallel default(shared)
  {
  E_Int size, initFrontSize;
  E_Float* edge;
  E_Int nof1, nof2, match, pos;
  list<E_Float*>::iterator it; 
  list<E_Float*>::iterator it1; list<E_Float*>::iterator it2;
  list<E_Int>::iterator iti;
  E_Int frontSize;
  E_Float minB[3]; E_Float maxB[3];
  E_Int* ptr; E_Int* ptr2;
  E_Int np, np2, ind, ret3, ret4;
  E_Float x, y, z, xp, yp, zp;
  E_Float pt1[3]; E_Float pt2[3]; E_Float pt3[3]; E_Float pt4[3];

#pragma omp for schedule(dynamic)
  for (nof1 = 0; nof1 < nfaces; nof1++)
  {
    //printf("face %d over %d faces.\n", nof1, nfaces);
    // BB de la face courante
    ptr = cnp+posFaces[nof1];
    np = ptr[0];
    /*
    minB[0] = K_CONST::E_MAX_FLOAT; minB[1] = K_CONST::E_MAX_FLOAT; minB[2] = K_CONST::E_MAX_FLOAT;
    maxB[0] = -K_CONST::E_MAX_FLOAT; maxB[1] = -K_CONST::E_MAX_FLOAT; maxB[2] = -K_CONST::E_MAX_FLOAT;
    for (E_Int j = 0; j < np; j++)
    {
      ind = ptr[j+1]-1;
      x = px[ind]; y = py[ind]; z = pz[ind];
      minB[0] = std::min(minB[0], x);
      minB[1] = std::min(minB[1], y);
      minB[2] = std::min(minB[2], z);
      maxB[0] = std::max(maxB[0], x);
      maxB[1] = std::max(maxB[1], y);
      maxB[2] = std::max(maxB[2], z);
    }
    */
    minB[0] = boxes[nof1]->minB[0];
    minB[1] = boxes[nof1]->minB[1];
    minB[2] = boxes[nof1]->minB[2];
    maxB[0] = boxes[nof1]->maxB[0];
    maxB[1] = boxes[nof1]->maxB[1];
    maxB[2] = boxes[nof1]->maxB[2];

    /*
    if (np == 4 && 
        minB[0] > 2.214 - 1.e-1 &&
        minB[1] > -2.76 - 1.e-1 &&
        minB[2] > -3.88 - 1.e-1 &&
        maxB[0] < 2.28 + 1.e-1 &&
        maxB[1] < -2.251 + 1.e-1 &&
        maxB[2] < -3.374 + 1.e-1) printf("ind=%d\n", nof1); // 16649
    */

    // Formation du front de F1 (front)
    list<E_Float*> front; list<E_Float*> newFront;
    list<E_Int> replacement; // indices des faces remplacant F1
    ind = ptr[np]-1;
    xp = px[ind]; yp = py[ind]; zp = pz[ind];
    for (E_Int j = 0; j < np; j++)
    {
      ind = ptr[j+1]-1;
      x = px[ind]; y = py[ind]; z = pz[ind];
      edge = new E_Float[7];
      edge[0] = x; edge[1] = y; edge[2] = z;
      edge[3] = xp; edge[4] = yp; edge[5] = zp;
      edge[6] = 0; // tag
      front.push_back(edge);
      xp = x; yp = y; zp = z;
    }
    initFrontSize = front.size();

    if (nof1 == 16649)
    {
      for (it1 = front.begin(); it1 != front.end(); it1++)
      {
        edge = *it1;
        /* printf("QUAD edge %f %f %f -> %f %f %f\n", edge[0],edge[1],edge[2],edge[3],edge[4],edge[5]); */
      }
    }

    // Recherche des faces intersectant la face (au sens BB Tree)
    // indicesBB est la liste des facettes candidates (contient elle-meme)
    vector<E_Int> indicesBB;
    bbtree->getOverlappingBoxes(minB, maxB, indicesBB);

    for (size_t k = 0; k < indicesBB.size(); k++) // faces candidates
    {
      nof2 = indicesBB[k];
      if (nof1 != nof2)
      {
        // liste des edges de F2
        list<E_Float*> F2edges;
        pos = posFaces[nof2];
        ptr2 = cnp+pos;
        np2 = ptr2[0];
        /*
        if (nof1 == 16649 && np2 == 3) printf("Testing face nof1=%d with nof2=%d\n", nof1, nof2);
        */
        ind = ptr2[np2]-1;
        xp = px[ind]; yp = py[ind]; zp = pz[ind];
        for (E_Int j = 0; j < np2; j++)
        {
          ind = ptr2[j+1]-1;
          x = px[ind]; y = py[ind]; z = pz[ind];
          edge = new E_Float[7];
          edge[0] = x; edge[1] = y; edge[2] = z;
          edge[3] = xp; edge[4] = yp; edge[5] = zp;
          edge[6] = 0; // tag
          F2edges.push_back(edge);
          xp = x; yp = y; zp = z;
        }
        
        /*
        if (nof1 == 16649 && np2 == 3)
        {
          for (it1 = F2edges.begin(); it1 != F2edges.end(); it1++)
          {
            edge = *it1;
            printf("TRI edge %f %f %f -> %f %f %f\n", edge[0],edge[1],edge[2],edge[3],edge[4],edge[5]);
          }
        }
        */
        //printf("done list F2: %d %d\n", (int)F2edges.size(), (int)front.size());

        // recherche les edges de F2 qui matchent le front
        match = 0;
        for (it1 = front.begin(); it1 != front.end(); it1++)
        {
          edge = (*it1);
          pt1[0] = edge[0]; pt1[1] = edge[1]; pt1[2] = edge[2];
          pt2[0] = edge[3]; pt2[1] = edge[4]; pt2[2] = edge[5];
          for (it2 = F2edges.begin(); it2 != F2edges.end(); it2++)
          {
            // is edge2 in edge1?
            edge = (*it2);
            pt3[0] = edge[0]; pt3[1] = edge[1]; pt3[2] = edge[2];
            pt4[0] = edge[3]; pt4[1] = edge[4]; pt4[2] = edge[5];

            ret3 = K_COMPGEOM::pointInSegment(pt1, pt2, pt3, tol);
            ret4 = K_COMPGEOM::pointInSegment(pt1, pt2, pt4, tol);
            //printf("1: %f %f %f ; %f %f %f\n", pt1[0], pt1[1],pt1[2], pt2[0], pt2[1],pt2[2]);
            //printf("2: %f %f %f ; %f %f %f\n", pt3[0], pt3[1],pt3[2], pt4[0], pt4[1],pt4[2]);
            //printf("ret %d %d\n", ret3, ret4);

            if (ret3 > 0 && ret4 > 0) match++;
          }
        }
        /*
        if (nof1 == 16649 && np2 == 3) printf("Edges match = %d\n", match);
        */
        list<E_Float*> added;
        if (match < 2) goto skip; // 2 edges sont requis au moins pour matcher
        
        for (it1 = front.begin(); it1 != front.end(); it1++)
        {
          edge = (*it1);
          pt1[0] = edge[0]; pt1[1] = edge[1]; pt1[2] = edge[2];
          pt2[0] = edge[3]; pt2[1] = edge[4]; pt2[2] = edge[5];
          for (it2 = F2edges.begin(); it2 != F2edges.end(); it2++)
          {
            // is edge2 in edge1?
            edge = (*it2);
            pt3[0] = edge[0]; pt3[1] = edge[1]; pt3[2] = edge[2];
            pt4[0] = edge[3]; pt4[1] = edge[4]; pt4[2] = edge[5];

            ret3 = K_COMPGEOM::pointInSegment(pt1, pt2, pt3, tol);
            ret4 = K_COMPGEOM::pointInSegment(pt1, pt2, pt4, tol);
            if (ret3 == 1 && ret4 == 2) { // full edge match
              /*
              if (nof1 == 16649 && np2 == 3) printf("Full edge\n");
              */
              (*it1)[6] = 1; (*it2)[6] = 1; //printf("full edge\n");
            }
            else if (ret3 == 2 && ret4 == 1) { // full edge match
              /*
              if (nof1 == 16649 && np2 == 3) printf("Full edge\n");
              */
              (*it1)[6] = 1; (*it2)[6] = 1; //printf("full edge\n");
            }
            else if (ret3 == 1 && ret4 == 3) { // partial edge match
              (*it1)[6] = 1; (*it2)[6] = 1;
              x = pt4[0]; y = pt4[1]; z = pt4[2];
              xp = pt2[0]; yp = pt2[1]; zp = pt2[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
            }
            else if (ret3 == 3 && ret4 == 2) { // partial match
              (*it1)[6] = 1; (*it2)[6] = 1;
              x = pt1[0]; y = pt1[1]; z = pt1[2];
              xp = pt3[0]; yp = pt3[1]; zp = pt3[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
            }
            else if (ret3 == 3 && ret4 == 3) { // partial match
              (*it1)[6] = 1; (*it2)[6] = 1;
              x = pt1[0]; y = pt1[1]; z = pt1[2];
              xp = pt3[0]; yp = pt3[1]; zp = pt3[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
              x = pt4[0]; y = pt4[1]; z = pt4[2];
              xp = pt2[0]; yp = pt2[1]; zp = pt2[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
            }
            else if (ret4 == 1 && ret3 == 3) { // partial match
              (*it1)[6] = 1; (*it2)[6] = 1;
              x = pt3[0]; y = pt3[1]; z = pt3[2];
              xp = pt2[0]; yp = pt2[1]; zp = pt2[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
            }
            else if (ret4 == 3 && ret3 == 2) { // partial match
              (*it1)[6] = 1; (*it2)[6] = 1;
              x = pt1[0]; y = pt1[1]; z = pt1[2];
              xp = pt4[0]; yp = pt4[1]; zp = pt4[2];
              edge = new E_Float[7];
              edge[0] = x; edge[1] = y; edge[2] = z;
              edge[3] = xp; edge[4] = yp; edge[5] = zp;
              edge[6] = 0; // tag
              added.push_back(edge);
            }
          }
        } // identification
        front.insert(front.end(), added.begin(), added.end());

        // nettoyage du front
        newFront.clear();
        for (it = front.begin(); it != front.end(); it++) 
        {
          if ((*it)[6] == 1) { delete [] *it; } 
          else newFront.push_back(*it);
        }
        front = newFront;
        
        replacement.push_back(nof2);
        //if (nof1 == 16649 && np2 == 3) printf("replaced with %d. Front restant=%d\n", nof2, front.size());
        skip:;
        // Clear les edges de F2
        for (it = F2edges.begin(); it != F2edges.end(); it++) delete [] *it;
        
      } // nof1 != nof2
    } // faces candidates

    frontSize = front.size();
    if (frontSize != initFrontSize && frontSize != 0) // partial match
    {
      ;
      //printf("face %d: partial match (front=%d) - no replace\n", nof1, front.size());
    }
    else if (frontSize != initFrontSize && frontSize == 0) // full match
    {
      /*
      if (nof1 == 16649 && np2 == 3)
      {
        printf("face %d: full match - replace\n", nof1);
        printf("face %d is replaced by faces : ", nof1);
        for (iti = replacement.begin(); iti != replacement.end(); iti++)
          printf("%d ", *iti);
        printf("\n");
      }
      */
      // replace
      size = replacement.size();
      E_Int* a = new E_Int [size+1];
      E_Int c = 1;
      a[0] = size;
      for (iti = replacement.begin(); iti != replacement.end(); iti++)
      {
        a[c] = (*iti); c++;
      }
      indir[nof1] = a;
    }
    else  ; //printf("face %d: no match\n", nof1);

    // Clear le front (si il en reste)
    for (it = front.begin(); it != front.end(); it++) delete [] *it;
    front.clear();
  }
  }

  // free boxes
  for (E_Int i = 0; i < nfaces; i++) delete boxes[i];
  delete bbtree;

  // reconstruction de la connectivite elt->faces
  E_Int* ptr = cn.begin();
  E_Int sizeEF = 0;
  E_Int sizeFN = ptr[1];
  ptr += 2+sizeFN;
  E_Int nelts = ptr[0];
  E_Int nf, nf2; 
  E_Int* addr;
  ptr += 2;
  for (E_Int i = 0; i < nelts; i++)
  {
    nf = ptr[0]; sizeEF++;
    for (E_Int j = 0; j < nf; j++) 
    { 
      addr = indir[ptr[j+1]-1];
      if (addr == NULL) sizeEF++;
      else sizeEF += addr[0];
    }
    ptr += nf+1;
  }
  //printf("sizeEF=%d, nelts=%d\n", sizeEF, nelts);
  
  cno = new FldArrayI(4 + sizeFN + sizeEF);
  E_Int* cnop = cno->begin();
  
  ptr = cn.begin();
  cnp = cn.begin();
  
  // Recopie la connectivite faces/noeuds
  for (E_Int i = 0; i < sizeFN+2; i++) cnop[i] = ptr[i];
  
  ptr = cnp+4+sizeFN; // pointe sur EF
  E_Int* ptr2 = cnop+sizeFN+2;
  ptr2[0] = nelts;
  ptr2[1] = sizeEF;
  ptr2 += 2;
 
  for (E_Int i = 0; i < nelts; i++)
  {
    nf = ptr[0]; nf2 = 0; 
    for (E_Int j = 0; j < nf; j++) 
    {
      addr = indir[ptr[j+1]-1];
      if (addr == NULL) { ptr2[nf2+1] = ptr[j+1]; nf2++; }
      else 
      {
        for (E_Int k = 0; k < addr[0]; k++)
        { ptr2[nf2+1] = addr[k+1]+1; nf2++; }
      }
    }
    ptr2[0] = nf2;
    ptr2 += nf2+1; ptr += nf+1;
  }
  
  // free indir
  for (E_Int i = 0; i < nfaces; i++)
  {
    if (indir[i] != NULL) delete [] indir[i];
  }
  delete [] indir;

  // DEBUG
  /*
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cno, cFE);
  printf("final %d %d\n", cFE(16649,1), cFE(16649,2));
  printf("final2 %d %d\n", cFE(13971,1), cFE(13971,2));
  printf("final3 %d %d\n", cFE(13969,1), cFE(13969,2));
  */
  
}
