# How to write universal array3 code

## Type of arrays

Arrays can be:
- structured
- unstructured:
  - BE/ME
  - NGON

## In place functions

Use getFromArray3 to pass from python object to FldArrays:
  
    PyObject* o;
    if (!PYPARSETUPLE_(args, O_, &o)) return NULL;
    E_Int ni, nj, nk;
    K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
    char* varString; char* eltType;
    E_Int ret = K_ARRAY::getFromArray3(o, varString, f, ni, nj, nk, c, eltType);

if ret = 1, it is a structured array
if ret = 2, it is an unstructured array.

Release memory. Dont forget at the end:

    RELEASESHAREDB(ret, o, f, c);
  
Universal access on field size:

    E_Int nfld = f->getNfld(); // nbre de champs
    E_Int npts = f->getSize(); // nbre de pts
  
Universal access on field pointers:

    E_Float* x = f->begin(1); 
    for (E_Int i = 0; i < 5; i++) printf(" " SF_F_ " ", x[i]);
    // modification
    x[0] = -0.05;

Universal access on field with operator()

    FldArrayF& fr = (*f);
    fr(1, 1) = +0.05;

Universal access on NGON connectivity (ret=2, eltType="NGON" or "NGON*"):

    E_Int nfaces = c->getNFaces();
    E_Int nelts = c->getNElts();

    E_Int* ngon = c->getNGon();
    E_Int* nface = c->getNFace();
    E_Int* indPG = c->getIndPG(); // always exists, must be called
    E_Int* indPH = c->getIndPH();

    // Universal acces on face 0
    E_Int size;
    E_Int* face = c->getFace(0, size, ngon, indPG);
    printf("face " SF_D_ ":", E_Int(0));
    for (E_Int i = 0; i < size; i++) printf(" " SF_D_ " ", face[i]);
    printf("\n");
    
    // Universal access on element 0
    E_Int* elt = c->getElt(0, size, nface, indPH);
    printf("elt " SF_D_ ":", E_Int(0));
    for (E_Int i = 0; i < size; i++) printf(" " SF_D_ " ", elt[i]);
    printf("\n");

Universal access on NGON sizes (may be usefull for allocation):

    E_Int sizeNGon = c->getSizeNGon(); // for NGONv3, contains face number and so is greater than for CGNSv4
    E_Int sizeNFace = c->getSizeNFace();
    

Universal access on BE/ME connectivity (ret=2, eltType="TRI","QUAD", "TRI,QUAD", "TRI*")

    E_Int nc = c->getNConnect();
    // in case of BE, nc is 1
    
    // universal eltType split
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    // access all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(c->getConnect(ic));
      E_Int nelts = cm.getSize(); // number of elements
      E_Int nvpe = cm.getNfld(); // number of vertices by element
      char* eltTypel = eltTypes[i]; // type of elements

      for (E_Int i = 0; i < nelts; i++)
      {
        cm(i,1) = 1; cm(i,2) = 2; cm(i,3) = 3; // dont use begin on connects
      }
    }

    // delete eltType split
    for (size_t i = 0; i < eltTypes.size(); i++) delete [] eltTypes[i];


Getting array api (internal storage information, rare usage):

    E_Int apif = f->getApi();

if apif=1, it is an array1
if apif=2, it is an array2 or an array3
For structured array, there is no difference between array2 and array3.

For NGONs:

    E_Int isNGon = c->getNGonType();
    // isNGon=1: NGON, NFACE CGNSv3 array1 compact
    // isNGON=2: NGON, NFACE, [indPG], [indPF] rake CGNSv3
    // isNGON=3: NGON, NFACE, indPG, indPF rake CGNSv4
  

## Create new array

Build empty arrays. Pass an api arg in your function.

Structured:

    PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, api);

Then get pointers to fill values:

    K_ARRAY::getFromArray3(tpl, f);

BE:

    PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", npts, ncells, 
                                         eltType, false, api);

Then get pointers to fill values:

    K_ARRAY::getFromArray3(tpl, f, cn);

ME:

    PyObject* K_ARRAY::buildArray3(3, "x,y,z", npts,
                                   std::vector<E_Int>& neltsPerConnect,
                                   "TRI,QUAD", false, api)

NGON:

    PyObject* K_ARRAY::buildArray3(3, "x,y,z", npts, ncells, nfaces,
                                  "NGON", sizeNGon, sizeNFace, ngonType,
                                  false, api);


## Copy functions

If input is api=1, build api=1 else build api=3.

Structured: build identical array to f. Varstring can change. If api=-1, use input api.

    PyObject* K_ARRAY::buildArray3(f, "x,y,z",
                                   ni, nj, nk, api=-1)

BE/ME/NGON: build identical array to f and cn. Varstring can change. If api=-1, use input api.

    PyObject* K_ARRAY::buildArray3(f, "x,y,z", varString, cn,
                                   "TRI", api=-1)