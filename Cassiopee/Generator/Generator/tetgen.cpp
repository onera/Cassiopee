/*    
    Copyright 2013-2024 Onera.

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
# include "generator.h"
using namespace std;
using namespace K_FLD;
# include "Tetgen/tetgen.h"
int __Error__ = 0;
jmp_buf __env__;

//=========================================================================
/* Generation de maillage tetra a partir d'un maillage surfacique (tetgen) 
   Le maillage surfacique n'est pas modifie */
//=========================================================================
PyObject* K_GENERATOR::tetgen(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float grading, maxh;
  E_Int remeshBoundaries;
  char* optionString;
  PyObject* holes; // coord. d'un point a l'interieur
  if (!PYPARSETUPLE_(args, O_ RR_ I_ O_ S_, 
                    &array, &maxh, &grading, &remeshBoundaries, &holes, &optionString))
    return NULL;

  // grading is not used
  // ambiguite sur maxh (max volume ou volume impose)

  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res <= 0) return NULL;
  if (res == 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: input must be TRI.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: input must be TRI.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "tetgen: coordinates not found in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  // Getting hole list
  E_Int nbholes = 0;
  E_Float* pholes = NULL;
  if (PyList_Check(holes) == true)
  {
    E_Int n = PyList_Size(holes);
    for (E_Int i = 0; i < n; i++)
    {
      PyObject* o = PyList_GetItem(holes, i);
      if (PyList_Check(o) == true && PyList_Size(o) == 3)
        nbholes++;
    }
    pholes = new E_Float [nbholes*3];
    for (E_Int i = 0; i < n; i++)
    {
      PyObject* o = PyList_GetItem(holes, i);
      for (E_Int j = 0; j < 3; j++)
        pholes[3*i+j] = PyFloat_AsDouble(PyList_GetItem(o, j));
    }
  }

  PyObject* pout = NULL;

  tetgenio in, out;
  tetgenbehavior b;

  // si les options sont passees par la ligne d'options
  // elles priment sur celles passees par arguments
  E_Int l = strlen(optionString); 
  if (l > 0) // options par optionString
  {
    b.parse_commandline(optionString);
  }
  else // options par arguments
  { 
    // Maille l'interieur d'un PLC
    b.plc = 1;
    // modification du maillage surfacique?
    if (remeshBoundaries == 0) b.nobisect = 1; // -Y (conserve la frontiere ext)
    else b.nobisect = 0;
    // b.nobisect_param = 1; // 0, 1 ou 2
    // b.addsteiner_algo = 1; // 1 ou 2
    // Ajoute des pts pour la qualite
    b.quality = 1; // -q (ajoute de points)

    //b.minratio = 1.2; // radius-edge ratio >= 0.612 - better control on aspect ratio
    
    //b.diagnose = 1;

    //b.regionattrib = 1; // met des attributs de region

    // options en ligne
    // option -q : ratio edge/radius (anisotropie?) default:2.0 (input)
    // option -a : volmax
    // option -r : remesh a previous mesh
    // option -Y : empeche le split de facettes
    // option -M : empeche la fusion des facettes coplanaires
    // option -T : tolerance pour les tests coplanaires 1.e-8
  
    if (maxh > 0)
    {
      b.fixedvolume = 1;
      b.maxvolume = (maxh*maxh*maxh)*0.25; // -a (max volume constraint)
    }
    else b.varvolume = 1;

    //b.varvolume = 0; // il faut garder 0 - CB

    // coplanar test
    b.epsilon = 1.e-10;

    // forced simple
    //b.nobisect = 0;
    //b.quality = 0;
    //b.fixedvolume = 0;
    //b.epsilon = 1.e-8;

    // pour du remaillage coarsening
    // b.refine = 1;
    // b.coarsen = 1;
    // b.coarsen_param = 1; // ??
    // b.coarsen_percent = 0.5; // x % less points
  }

  // print b
  /*
  printf("plc %d\n", b.plc);                                      // '-p', 0.
  printf("psc %d\n", b.psc);                                      // '-s', 0.
  printf("refine %d\n", b.refine);                                // '-r', 0.
  printf("quality %d\n", b.quality);                              // '-q', 0.
  printf("nobisect %d\n", b.nobisect);                            // '-Y', 0.
  printf("coarsen %d\n", b.coarsen);                              // '-R', 0.
  printf("weighted %d\n", b.weighted);                            // '-w', 0.
  printf("brio_hilbert %d\n", b.brio_hilbert);                    // '-b', 1.
  printf("incrflip %d\n", b.incrflip);                            // '-l', 0.
  printf("flipinsert %d\n", b.flipinsert);                        // '-L', 0.
  printf("metric %d\n", b.metric);                                // '-m', 0.
  printf("varvolume %d\n", b.varvolume);                          // '-a', 0.
  printf("fixedvolume %d\n", b.fixedvolume);                      // '-a', 0.
  printf("regionattrib %d\n", b.regionattrib);                    // '-A', 0.

  printf("conforming %d\n", b.conforming);                        // '-D', 0.
  printf("insertaddpoints %d\n", b.insertaddpoints);              // '-i', 0.
  printf("diagnose %d\n", b.diagnose);                            // '-d', 0.
  printf("convex %d\n", b.convex);                                // '-c', 0.
  printf("nomergefacet %d\n", b.nomergefacet);                    // '-M', 0.
  printf("nomergevertex %d\n", b.nomergevertex);                  // '-M', 0.
  printf("noexact %d\n", b.noexact);                              // '-X', 0.
  printf("nostaticfilter %d\n", b.nostaticfilter);                // '-X', 0.
  printf("zeroindex %d\n", b.zeroindex);                          // '-z', 0.
  printf("facesout %d\n", b.facesout);                            // '-f', 0.
  printf("edgesout %d\n", b.edgesout);                            // '-e', 0.
  printf("neighout %d\n", b.neighout);                            // '-n', 0.
  printf("voroout %d\n", b.voroout);                              // '-v', 0.
  printf("meditview %d\n", b.meditview);                          // '-g', 0.
  printf("vtkview %d\n", b.vtkview);                              // '-k', 0.
  printf("nobound %d\n", b.nobound);                              // '-B', 0.
  printf("nonodewritten %d\n", b.nonodewritten);                  // '-N', 0.
  printf("noelewritten %d\n", b.noelewritten);                    // '-E', 0.
  printf("nofacewritten %d\n", b.nofacewritten);                  // '-F', 0.
  printf("noiterationnum %d\n", b.noiterationnum);                // '-I', 0.
  printf("nojettison %d\n", b.nojettison);                        // '-J', 0.
  printf("reversetetori %d\n", b.reversetetori);                  // '-R', 0.
  printf("docheck %d\n", b.docheck);                              // '-C', 0.
  printf("quiet %d\n", b.quiet);                                  // '-Q', 0.
  printf("verbose %d\n", b.verbose);                              // '-V', 0.

  // Parameters of TetGen. 
  printf("vertexperblock %d\n", b.vertexperblock);                 // '-x', 4092.
  printf("tetrahedraperblock %d\n", b.tetrahedraperblock);         // '-x', 8188.
  printf("shellfaceperblock %d\n", b.shellfaceperblock);           // '-x', 2044.
  printf("nobisect_param %d\n", b.nobisect_param);                 // '-Y', 2.
  printf("addsteiner_algo %d\n", b.addsteiner_algo);               // '-Y/', 1.
  printf("coarsen_param %d\n", b.coarsen_param);                   // '-R', 0.
  printf("weighted_param %d\n", b.weighted_param);                 // '-w', 0.
  printf("fliplinklevel %d\n", b.fliplinklevel);                   // -1.
  printf("flipstarsize %d\n", b.flipstarsize);                     // -1.
  printf("fliplinklevelinc %d\n", b.fliplinklevelinc);             //  1.
  printf("reflevel %d\n", b.reflevel);                             // '-D', 3.
  printf("optlevel %d\n", b.optlevel);                             // '-O', 2.
  printf("optscheme %d\n", b.optscheme);                           // '-O', 7.
  printf("delmaxfliplevel %d\n", b.delmaxfliplevel);               // 1.
  printf("order %d\n", b.order);                                   // '-o', 1.
  printf("steinerleft %d\n", b.steinerleft);                       // '-S', 0.
  printf("no_sort %d\n", b.no_sort);                               // 0.
  printf("hilbert_order %d\n", b.hilbert_order);                   // '-b///', 52.
  printf("hilbert_limit %d\n", b.hilbert_limit);                   // '-b//'  8.
  printf("brio_threshold %d\n", b.brio_threshold);                 // '-b' 64.
  printf("brio_ratio %g\n", b.brio_ratio);                         // '-b/' 0.125.
  printf("facet_ang_tol %g\n", b.facet_ang_tol);                   // '-p', 179.9.
  printf("maxvolume %g\n", b.maxvolume);                           // '-a', -1.0.
  printf("minratio %g\n", b.minratio);                             // '-q', 0.0.
  printf("mindihedral %g\n", b.mindihedral);                       // '-q', 5.0.
  printf("optmaxdihedral %g\n", b.optmaxdihedral);                 // 165.0.
  printf("optminsmtdihed %g\n", b.optminsmtdihed);                 // 179.0.
  printf("optminslidihed %g\n", b.optminslidihed);                 // 179.0.  
  printf("epsilon %g\n", b.epsilon);                               // '-T', 1.0e-8.
  printf("minedgelength %g\n", b.minedgelength);                   // 0.0.
  printf("coarsen_percent %g\n", b.coarsen_percent);               // -R1/#, 1.0.
  */

  // Remplissage de in a partir de array
  in.mesh_dim = 3;
  in.numberofpointattributes = 0;

  E_Int n = f->getSize();
  E_Float* ptx = f->begin(posx);
  E_Float* pty = f->begin(posy);
  E_Float* ptz = f->begin(posz);
  in.numberofpoints = n;
  in.pointlist = new E_Float [n*3];
  E_Float* ptin = in.pointlist;
  for (E_Int i = 0; i < n; i++)
  {
    ptin[0] = ptx[i]; ptin[1] = pty[i]; ptin[2] = ptz[i]; ptin += 3;
  }
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  E_Int ne = cn->getSize();
  in.numberoffacets = ne;
  in.facetlist = new tetgenio::facet[ne];
  for (E_Int i = 0; i < ne; i++)
  {
    tetgenio::facet* f = &(in.facetlist[i]);
    tetgenio::init(f);
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[1];
    tetgenio::polygon* p = &f->polygonlist[0];
    tetgenio::init(p);
    p->numberofvertices = 3;
    p->vertexlist = new int[p->numberofvertices];
    
    p->vertexlist[0] = cn1[i]-1;
    p->vertexlist[1] = cn2[i]-1;
    p->vertexlist[2] = cn3[i]-1;
  }

  // Ajoute des holes
  // in.numberofholes = 1;
  // in.holelist = new double [3* in.numberofholes];
  // in.holelist[0] = 0.5; // x
  // in.holelist[1] = 0.5; // y
  // in.holelist[2] = 0.5; // z

  in.numberofholes = nbholes;
  in.holelist = new double [3* in.numberofholes];
  for (E_Int i = 0; i < nbholes; i++)
  {
    in.holelist[3*i+0] = pholes[3*i+0];
    in.holelist[3*i+1] = pholes[3*i+1];
    in.holelist[3*i+2] = pholes[3*i+2];
  }
  delete [] pholes;
  

  // Ajoute les regions lists
  //in.numberofregions = 1;
  //in.regionlist = new double [5* in.numberofregions];
  //in.regionlist[0] = 0.5; // x
  //in.regionlist[1] = 0.5; // y 
  //in.regionlist[2] = 0.5; // z
  //in.regionlist[3] = 1; // entier reperant la region (region attrib)
  //in.regionlist[4] = 0.01; // ?

  RELEASESHAREDU(array, f, cn);
  
  // in: surface mesh or previous tetra mesh
  // out: output mesh
  // addin: if != NULL, list of constraint points
  // bgmin: if != NULL, background mesh with size function
  int val;
  val = setjmp(__env__); // enregistre l'environnement pour retour d'erreur
  if (val == 0) // premier passage
    tetrahedralize(&b, &in, &out, NULL, NULL);

  if (val == 1) // erreur
  { 
    PyErr_SetString(PyExc_TypeError, "tetgen: failed.");
    return NULL;
  }

  /* Build output */
  E_Int np = out.numberofpoints;
  ne = out.numberoftetrahedra;
  printf("INFO: tetgen: generate %d points and %d elements.\n",np,ne);
  pout = K_ARRAY::buildArray(3, "x,y,z", np, ne, 4, NULL, false, 0);
  E_Float* fp = K_ARRAY::getFieldPtr(pout);
  E_Float* fx = fp; E_Float* fy = fp+np; E_Float* fz = fp+2*np;
  E_Float* pt = out.pointlist;
  for (E_Int i = 0; i < np; i++)
  {
    fx[i] = pt[0]; fy[i] = pt[1]; fz[i] = pt[2]; pt += 3;
  }
  E_Int* cp = K_ARRAY::getConnectPtr(pout);
  E_Int* cp1 = cp; E_Int* cp2 = cp+ne; 
  E_Int* cp3 = cp+2*ne; E_Int* cp4 = cp+3*ne;
  E_Int ind = 0;
  int* tl = out.tetrahedronlist;
  for (E_Int i = 0; i < ne; i++)
  { 
    cp1[i] = tl[ind]+1;
    cp2[i] = tl[ind+1]+1;
    cp3[i] = tl[ind+2]+1;
    cp4[i] = tl[ind+3]+1;
    ind += 4;
  }

  return pout;
}
