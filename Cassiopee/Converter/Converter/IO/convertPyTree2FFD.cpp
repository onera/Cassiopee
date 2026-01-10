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
#include "PyTree/PyTree.h"
#include "converter.h"
#include <iostream>
extern "C"
{
  void writeffdfilef_(E_Int& n, E_Int* int_kode, E_Int& nnode, E_Int& nelmt,
                      E_Float* real_state1, E_Float* real_state2,
            		      E_Float* Coord_x, E_Float* Coord_y, E_Float* Coord_z,
		                  E_Float* Density, E_Float* MomentumX, E_Float* MomentumY, E_Float* MomentumZ, E_Float* Energy, E_Float* MutsMu, E_Float* TurbulentDistance,
  		                E_Int& nlimt, E_Float* Var_l,
		                  E_Int& nmtch, E_Int* ielmtmtch1, E_Int* ielmtmtch2, E_Int* izoneznzn, E_Int* ielmtznzn,
                      E_Int& npoint, E_Int* nodmtch, E_Int* kpoinmtch);
}
# include <string.h>
using namespace std;
using namespace K_FLD;

E_Int test = 0; 
E_Int	kodnst;

// gestion des impression de contrôle : insérer test = 0/1 permet de contrôler les impressions de debug

//=============================================================================
// Recupere une zone (les pointeurs sont partages avec python)
//=============================================================================
// zone doit être de type non-structurée NGON et contenir les grandeurs suivantes :
//    - les coordonnées (CoordinateX, CoordinateY,CoordinateZ)
//    - la densité (Density) 
//    - les variables conservatives (MomentumX, MomentumY, MomentumZ)
//    - l'energie totale , ro.E (EnergyStagnationDensity )
//    - le rapport de la viscosité turbulente sur viscosité moléculaire (ViscosityEddy/ViscosityMolecular)
//    - la distance à la paroi (TurbulentDistance)
//    - les informations de connectivité :
//        - NGonElements { ElementRange, ElementConnectivity, ParentElements}
//        - NFaceElements  { ElementRange, ElementConnectivity }
//    - les informations sur les conditions limites (ZoneBC)
//    - les valeurs des variables au centre des cellules de condition limite dans BCdata
//  RefStat est un numpy. Il doit permettre de passer :
//  xmach,alpha,roinf,vinf,tinf,pinf,reynolds lref,aref,alref,beta dans real_state1
//  et gam,rgas,suth,tsuth dans real_state2
// • xmach: Mach number
// • alpha: angle of attack (degrees)
// • beta: necessarily beta = 0. (yaw angle)
// • roinf: density in the reference state
// • vinf: velocity in the reference state
// • tinf: static temperature in the reference state
// • pinf: static pressure in the reference state
// • reynolds lref: Reynolds number based on reference length
// • alref: reference length
// • aref: reference area
// • gam: ratio of specific heats
// • rgas: gas constant
// • suth, tsuth: constant for the Sutherland law and associated temperature
//
//   FlowEq est un numpy portant les informations nécessaires pour kod_int
//     kod_int =(kodcc, kodsolver,kodnst,kod2d)
//    • kodcc: kodcc = −1 for a cell-vertex solution, kodcc = +1 for a cell-centred solution
//    • kodsolver: kodsolver = 0 for a TAU solution with the Spalart turbulence model, kodsolver = 1 for a
//     TAU solution with the k −! turbulence model or a solver producing spurious turbulence in the far-field,
//    • kodsolver = 2 for an elsA solution
//    • kodnst: kodnst = +1 for a R.A.N.S solution, kodnst = −1 for an Euler solution
//    • kod2d: kod2d = +1 for a 2D solution, kod2d = −1 for a 3D solution
PyObject* K_CONVERTER::convertPyTree2FFD(PyObject* self, PyObject* args)
{
  PyObject* zone;
  PyObject* RefStat;
  PyObject* FlowEq;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  E_Int nd , ii;
  E_Int neq = 5 ;

  if (!PYPARSETUPLE_(args, OOO_ I_ SSS_, 
                     &zone, &RefStat, &FlowEq, &nd,                 
                     &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  if( test == 1 ){
     printf("je rentre dans convertPyTree2FFD\n");
     printf("nd = " SF_D_ "\n", nd);
   }
  vector<PyArrayObject*> hook;
  E_Float real_state1[10];
  E_Float real_state2[4];
  E_Int kod_int [4];
/*---------------------------------------------------------------*/
/* transfert de l'état de référence et des codes dans tableaux pour fortran.  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

     E_Int lrefstat = PyList_Size(RefStat);
     if( test == 1 )printf("longueur de RefStat = " SF_D_ "\n", lrefstat);
     for  ( ii = 0 ; ii < lrefstat ; ii++)
        {
      PyObject* fl   = PyList_GetItem(RefStat, ii);
      if( ii < 10 )real_state1[ii]     = PyFloat_AS_DOUBLE(fl) ;
      if( ii >  9 )real_state2[ii-10]  = PyFloat_AS_DOUBLE(fl) ;
        }
     for  ( ii = 0 ; ii < 4; ii++)
        {
      PyObject* fl   = PyList_GetItem(FlowEq, ii);
      kod_int[ii]     = PyLong_AsLong(fl) ;
        }
 kodnst =  kod_int[2] ;
 if(kodnst > 0)neq=neq+1 ;
  /*--------------------------------------------*/
  /* zone a ecrire                              */
  /*--------------------------------------------*/
  E_Int nnode, nelmt, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int loc = 2;

  E_Int res = K_PYTREE::getFromZone(zone, 1, loc, varString,
                        fields, locs, nnode, nelmt, km,
                        cn, cnSize, cnNfld, eltType, hook, 
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
   if( test == 1 ){
     printf("getFromZone a repondu :\n");
     printf("res = " SF_D_ "\n", res);
     printf("varString : %s\n", varString);
     printf("nnode = " SF_D_ ", nelmt = " SF_D_ ", km = " SF_D_ " \n", nnode, nelmt, km);
     printf("cnSize= " SF_D_ ", cnNfld= " SF_D_ " \n", cnSize, cnNfld);
     //printf("cn[0]= %d, cn2[0]= %d \n",cn[0],cn2[0]);
   }
    /* Plus d'info dans KCore/PyTree/PyTree.h */
    if (res != 2) 
    {
      RELEASESHAREDZ(hook, varString, eltType);
      PyErr_SetString(PyExc_TypeError, 
                      "convertPyTree2FFD: requires a NGON.");
      return NULL;
    }
    else
    {
      if( test == 1 )printf("Zone non structuree de type %s.\n", eltType);
    }
   // extraction des coordonnées
    //printf("on va chercher ...." ) ;
//    cerr << "on va chercher ...." << endl ;
    E_Int posx, posy, posz;
    E_Float* Coord_x ; E_Float* Coord_y ; E_Float* Coord_z ;
//    cerr << "on va chercher CoordinateX" << endl ;
    posx = K_ARRAY::isCoordinateXPresent(varString);
//    cerr << "on va chercher CoordinateY" << endl ;
    posy = K_ARRAY::isCoordinateYPresent(varString);
//    cerr << "on va chercher CoordinateZ" << endl ;
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
         PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD :zones do not have coordinates. Not written.\n");
    return NULL; 
    }
//    cerr <<"Les positions des coordonnées sont x :"<<posx <<" y : "<<posy<<" , z :"<<posz<< endl;
     Coord_x = fields[posx];
     Coord_y = fields[posy];
     Coord_z = fields[posz];
//     cerr <<"c'est bon ....."<<endl ;
   // extraction des variables
    E_Float* Density; E_Float* MomentumX; E_Float* MomentumY; E_Float* MomentumZ; E_Float* Energy;
    E_Float* TurbulentDistance = NULL;
    vector<E_Float> MutsMu(nelmt);
// Recherche des variables  (reutilisation de posx, posy posz) 
    E_Int posd =  K_ARRAY::isDensityPresent(varString);
    if (posd == -1 )
    {
         PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable Density is missing. Not written.\n");
    return NULL; 
    }
//     cerr <<"La position de Density est  :"<<posd << endl;
    Density = fields[posd];
//    cerr <<"c'est bon ....."<<endl ;
//
//    posx = K_ARRAY::isVelocityXPresent(varString);
//    posy = K_ARRAY::isVelocityYPresent(varString);
//    posz = K_ARRAY::isVelocityZPresent(varString);
    posx = K_ARRAY::isMomentumXPresent(varString);
    posy = K_ARRAY::isMomentumYPresent(varString);
    posz = K_ARRAY::isMomentumZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
         PyErr_SetString(PyExc_TypeError,
                      "convertPyTree2FFD :zones do not have Momentum. Not written.\n");
      return NULL; 
    }
    MomentumX = fields[posx];
    MomentumY = fields[posy];
    MomentumZ = fields[posz];
// equation d'energie
    posd =  K_ARRAY::isEnergyStagnationDensityPresent(varString);
//    posd =  K_ARRAY::isTemperaturePresent(varString);
      if (posd == -1 )
      {
         PyErr_SetString(PyExc_TypeError,
                 "convertPyTree2FFD : zone do not have EnergyStagnationDensity. Not written.\n");
      return NULL; 
      }
    Energy = fields[posd];
//    cerr<<"Les coordonnées et les grandeurs du champ moyen ont été lues"<<endl;
// 
// Calcul du champ de Mut/Mu
// existance de champs de viscosité moléculaire ou cinématique
//      cerr<<"kodnst ="<< kodnst <<endl;
      if( kodnst == 1){
//    cerr<<"l y a des grandeurs turbulentes"<<endl;
            E_Int posVM;  E_Int posVK;  E_Int posVE;  E_Int posVEK;
//
           posVM = K_ARRAY::isNamePresent("ViscosityMolecular", varString);  //Mu
           posVK = K_ARRAY::isNamePresent("ViscosityKinetic", varString); // Nu
           if( posVM == -1 && posVK == -1)
             {
              PyErr_SetString(PyExc_TypeError,
               "convertPyTree2FFD : zone do not have neither ViscosityMolecular neither ViscosityKinetic. Not written.\n");
               return  NULL;
              }
           else {
// existance de champs de viscosité turbulente ( Nut ou Mut )
           posVE = K_ARRAY::isNamePresent("ViscosityEddy", varString); // Mut
           posVEK = K_ARRAY::isNamePresent("ViscosityEddyKinetic", varString); //Nut
            if( posVM == -1 && posVK == -1 )
              {
               PyErr_SetString(PyExc_TypeError,
               "convertPyTree2FFD : zone do not have neither ViscosityEddy neither ViscosityEddyKinetic. Not written.\n");
               return  NULL;
              }
              }
             if ( posVM != -1 && posVE != -1) {
                  E_Float* Mut = fields[posVE] ;
                  E_Float* Mu  = fields[posVM] ;
                  K_CONVERTER::calculateMutsMu(nelmt,&MutsMu[0],Mut,Mu);
               }
             else if ( posVEK != -1 && posVK != -1) {
                  E_Float* Mut = fields[posVEK] ;
                  E_Float* Mu  = fields[posVK] ;
                  K_CONVERTER::calculateMutsMu(nelmt,&MutsMu[0],Mut,Mu);
               }
             else if ( posVM != -1 && posVEK != -1) {
                  E_Float* Mut = fields[posVEK] ;
                  E_Float* Mu  = fields[posVM] ;
                  K_CONVERTER::calculateMutsMu(nelmt,&MutsMu[0],Mut,Mu,Density);  //MutsMu= Ro*Vek/vm
               }
             else if ( posVE != -1 && posVK != -1) {
                  E_Float* Mut = fields[posVEK] ;
                  E_Float* Mu  = fields[posVM] ;
                  K_CONVERTER::calculateMutsMu(nelmt,&MutsMu[0],Mut,Mu,NULL,Density); //MutsMu= Ve/(Ro*Vk)
               }
             else
               {             
              PyErr_SetString(PyExc_TypeError,
               "convertPyTree2FFD : cannot calculate MutsMu. Not written.\n");
               return  NULL;
               }
//Distance paroi
    posd = K_ARRAY::isNamePresent("TurbulentDistance", varString);
    if (posd == -1 )
    {
         PyErr_SetString(PyExc_TypeError,
               "convertPyTree2FFD : zone do not have Turbulent Distance. Not written.\n");
    return NULL; 
    }
   TurbulentDistance = fields[posd];
//    cerr<<"Les grandeurs turbulentes ont été lues"<<endl;
    }
//
//
// traitement des BC et des faces
//
//   cerr<<"Début de traitement des interfaces"<<endl;
   E_Int nmtch, s1, nBC = 0;
   PyObject* NG = K_PYTREE::getNodeFromName1(zone, "NGonElements");
   PyObject* PE = K_PYTREE::getNodeFromName1(NG,  "ParentElements");
   E_Int* ielmtmtch1 =  K_PYTREE::getValueAI(PE, nmtch, s1, hook);
//   E_Int* p =  K_PYTREE::getValueAI(PE, nmtch, s1, hook);

//   cerr << "appel de copyPE2Elmtmtch s1="<<s1<<endl;
   E_Int* ielmtmtch2 = ielmtmtch1 + nmtch;
   E_Int* p = NULL;
//   E_Int  ielmtmtch2[nmtch];
//   E_Int  ielmtmtch1[nmtch];
   K_CONVERTER::copyPE2Elmtmtch( p, ielmtmtch1, ielmtmtch2, nmtch, &nBC ) ; 
//   cerr << "sortie de copyPE2Elmtmtch"<<endl;
//
// transfert des BC dans ielmtmtch2.
// transfert des variables de BC vers tableaux pour écriture
//  
// calcul taille des tableaux BC
//    cerr << "appel de scanBC"<<endl ;
   E_Int nlimt = 0 ;
   K_CONVERTER::scanBC(zone, &nlimt) ;
//    cerr << " taille des tableaux BC :"<< nlimt <<endl ;
//
//   E_Float Density_l[nlimt]; E_Float MomentumX_l[nlimt]; E_Float MomentumY_l[nlimt]; E_Float MomentumZ_l[nlimt];
//   E_Float Energy_l[nlimt]; E_Float MutsMu_l[nlimt];
   vector<E_Float> Var_l(neq*nlimt);
   K_CONVERTER::getVarBC(zone, &Var_l[0], ielmtmtch2, &nlimt) ;
//
//  Connectivité - désenlacement de ElementConnectivity
    E_Int npoint ;
    E_Int NElCon = 0 ;
    //E_Int s2;

    // cnNfld est > npoint, mais pas toujours !!!
    //PyObject* NF = K_PYTREE::getNodeFromName1(zone, "NFaceElements");
    //PyObject* EC = K_PYTREE::getNodeFromName1(NF,  "ElementConnectivity");
    //E_Int* q =  K_PYTREE::getValueAI(EC, NElCon, s2, hook);
    //cerr << "NElCon="<<NElCon<<" s2="<<s2<<endl;
   // cnNfld est > npoint 
    vector<E_Int> nodmtch(NElCon) ; vector<E_Int> kpoinmtch(NElCon) ;
    cerr << "appel de splitElementConnectivity"<<endl ;
    K_CONVERTER::splitElementConnectivity(zone, &npoint, &nodmtch[0], &kpoinmtch[0]) ;
    vector<E_Int> izoneznzn (nlimt) ; vector<E_Int> ielmtznzn (nlimt);
    for (E_Int i = 0 ; i < nlimt ; i++) {  izoneznzn [i] = 0 ;  ielmtznzn [i] = 0 ;}
//
// appel du sous-programme d'ecriture fortran
// kod_int =(kodcc, kodsolver,kodnst,kod2d)
    writeffdfilef_(nd, kod_int, nnode, nelmt,
                   real_state1, real_state2,
		   Coord_x, Coord_y, Coord_z,
		   Density, MomentumX, MomentumY, MomentumZ, Energy,
                   &MutsMu[0], TurbulentDistance,
                   nlimt, &Var_l[0],
                   nmtch, ielmtmtch1, ielmtmtch2, &izoneznzn[0], &ielmtznzn[0],
                   npoint, &nodmtch[0], &kpoinmtch[0]) ; 

//   fin de convertPyTree2FFD 
     if( test == 1 ){
     printf("je sors de convertPyTree2FFD\n");
     printf("nd = " SF_D_ "\n", nd);
     }
    RELEASESHAREDZ(hook, varString, eltType);
    Py_INCREF(Py_None);
    return Py_None;
}
//
//
void K_CONVERTER::calculateMutsMu(E_Int nelmt, E_Float* MutsMu, E_Float* Mut, E_Float* Mu, E_Float* Density, E_Float* Densitym1)
{
      for ( E_Int i = 0 ;  i < nelmt ; i++){
           MutsMu[i] = Mut[i]/Mu[i] ;
      if ( Density != NULL )
           MutsMu[i] = MutsMu[i]*Density[i] ;
      if ( Densitym1 != NULL )
           MutsMu[i] = MutsMu[i]/Density[i] ;
      }
}
//
//
void K_CONVERTER::scanBC(PyObject* zone,  E_Int* nlimt)
{
// 
    PyObject* zonebc; 
    E_Int nb_bc=0;
    vector<PyArrayObject*> hook;
    //E_Int type_bc;
    *nlimt = 0;
    vector<PyObject*> data_set;
    K_PYTREE::getNodesFromType1(zone, "ZoneBC_t", data_set);
    zonebc = data_set[0];
//   zonebc   = K_PYTREE::getNodeFromName1(zone, "ZoneBC" );
    if (zonebc == NULL)
    {  printf("pas de ZoneBC \n"); 
    }
    else
    {
     PyObject* list_bc = PyList_GetItem(zonebc, 2); 
//     cerr<< " list_bc ="<< list_bc <<endl ;

     nb_bc = PyList_Size(list_bc);

     PyObject* bc; PyObject* node; char* str;

    // boucle sur les noeuds d'arbre contenus dans zone BC
    for (E_Int ibc = 0; ibc < nb_bc; ibc++)
    {
      bc   = PyList_GetItem(list_bc, ibc);
      node = PyList_GetItem(bc, 3);
      str = NULL;
      if (PyString_Check(node)) str = PyString_AsString(node); // type_bc
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
      if (K_STRING::cmp(str, "BC_t") == 0)
      {
        //E_Int s;
        E_Int PLSize = 0;
        E_Int s1;
        //char* st = K_PYTREE::getValueS(bc, s, hook);
        PyObject* PL = K_PYTREE::getNodeFromName1(bc, "PointList");
        K_PYTREE::getValueAI(PL, s1, PLSize, hook);
        *nlimt += PLSize;
       }
      }
      }
}
//
//void K_CONVERTER::getVarBC(PyObject* zone, E_Float* Density_l, E_Float* MomentumX_l, E_Float* MomentumY_l, E_Float* MomentumZ_l, E_Float* Energy_l, E_Float* MutsMu_l, E_Int* ielmtmtch2, E_Int* nlimt)
void K_CONVERTER::getVarBC(PyObject* zone, E_Float* Var_l, E_Int* ielmtmtch2, E_Int* nlimt)
{
    PyObject* zonebc; E_Int nb_bc=0; E_Int ielmbc = 0;
    vector<PyArrayObject*> hook;
// les types de conditions limites ( http://cgns.github.io/CGNS_docs_current/sids/bc.html#BC §9.7 )
/*  BCType_t := Enumeration(
    BCTypeNull, BCTypeUserDefined, BCAxisymmetricWedge, BCDegenerateLine.
    BCDegeneratePoint, BCDirichlet, BCExtrapolate, BCFarfield, BCGeneral,
    BCInflow, BCInflowSubsonic, BCInflowSupersonic, BCNeumann,
    BCOutflow, BCOutflowSubsonic, BCOutflowSupersonic, BCSymmetryPlane,
    BCSymmetryPolar, BCTunnelInflow, BCTunnelOutflow, BCWall,
    BCWallInviscid, BCWallViscous, BCWallViscousHeatFlux,
    BCWallViscousIsothermal, FamilySpecified ) ;*/
//
    E_Int type_bc ;
    E_Int PLSize = 0;
    zonebc   = K_PYTREE::getNodeFromName1(zone, "ZoneBC" );
    if(zonebc == NULL)
    {  printf("pas de ZoneBC \n"); 
    }
    else
    {
//     cerr <<"récupération des informations dans l'arbre"<<endl ;
     PyObject* list_bc = PyList_GetItem(zonebc, 2); 
     nb_bc =             PyList_Size(list_bc);

     PyObject* bc; PyObject* node; char* str;

    // boucle sur les noeuds d'arbre contenus dans zone BC
//     cerr <<"boucle dans zone : nb_bc="<< nb_bc<<endl ;
    for (E_Int ibc = 0; ibc < nb_bc; ibc++)
    {
      bc   = PyList_GetItem(list_bc, ibc);
      node = PyList_GetItem(bc, 3);
      str = NULL;
      if (PyString_Check(node)) str = PyString_AsString(node); // type_bc
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
#endif
      if (K_STRING::cmp(str, "BC_t") == 0)
      {
        E_Int s;
        E_Int s1;
        char* st = K_PYTREE::getValueS(bc, s, hook);
//        cerr <<" getVarBC : st ="<<st <<endl ;
        // put FFD72 tags for BCs ielmtmtch2[]
        //
        type_bc = -9999 ;
        if      (K_STRING::cmp(st, s, "BCFarfield")              == 0) type_bc = -1;
        else if (K_STRING::cmp(st, s, "BCWall")                  == 0) type_bc = 0;
        else if (K_STRING::cmp(st, s, "BCWallInviscid")          == 0) type_bc = 0;
        else if (K_STRING::cmp(st, s, "BCWallViscous")           == 0) type_bc = 0;
        else if (K_STRING::cmp(st, s, "BCOutflow")               == 0) type_bc =-10;
        else if (K_STRING::cmp(st, s, "BCOutflowSubsonic")       == 0) type_bc =-10;
        else if (K_STRING::cmp(st, s, "BCInflowSubsonic")        == 0) type_bc =-20;
        else if (K_STRING::cmp(st, s, "BCautoperiod")            == 0) type_bc =-3;
        else if (K_STRING::cmp(st, s, "BCSymmetryPlane")         == 0) type_bc =-2;
        else if (K_STRING::cmp(st, s, "BCSymmetryPolar")         == 0) type_bc =-2;


         PyObject* PL = K_PYTREE::getNodeFromName1(bc, "PointList");
         E_Int* p =  K_PYTREE::getValueAI(PL, s1, PLSize, hook);
//        cerr <<" getVarBC : PLSize ="<<PLSize <<endl ;
        for (E_Int i = 0; i < PLSize; i++)
            {
            E_Int ind = *(p++) -1 ;
            //if( test == 1)printf("avant : ielmtmtch2[%d] = %d  ", ind, ielmtmtch2[ind]);
            ielmtmtch2[ind] = type_bc ;
            //if( test == 1){printf("après : ielmtmtch2[%d] = %d  \n", ind, ielmtmtch2[ind]);fflush(stdout) ;}
            } 

           cerr <<" getVarBC : traitement des données= "<<PLSize <<endl ;
           vector<PyObject*> data_set;
           vector<PyObject*> data;
           //PyObject* zonebc;
           K_PYTREE::getNodesFromType1(bc, "BCDataSet_t", data_set );
            if ( data_set.size() != 0 )
            {
            cerr << data_set.size() << " BCDataset present " << endl;
            K_PYTREE::getNodesFromType1(data_set[0], "BCData_t", data );
            }
            else
            {
            cerr<< "BCDataSet absent"<<endl;
            K_PYTREE::getNodesFromType1(bc, "BCData_t", data );
            }
//           K_PYTREE::getNodesFromType1(bc, "BCData_t", data );
           
           if ( data.size() != 0 )
           {
            cerr << data.size() << " BCData présent " << endl;
           PyObject* bcdata = data[0] ;
//           cerr <<" getVarBC :  type_bc = "<< type_bc  <<endl ;
           PyObject* list_bcdata = PyList_GetItem(bcdata, 2); 
//           cerr <<" getVarBC :  list_bcdata = "<< list_bcdata <<endl ;
           E_Int nb_bcdata       = PyList_Size(list_bcdata);
//           cerr <<" getVarBC :  nb_bcdata = "<< nb_bcdata <<endl ;
           if( nb_bcdata == 1 ){ 
//            cerr <<" getVarBC :  nombre de données nb_bcdata = "<< nb_bcdata<< " insuffisant" <<endl ;
           }
           else
           {
           PyObject*          t  = K_PYTREE::getNodeFromName1(bcdata, "Density");
           E_Float* iptro;  E_Float* iptmx; E_Float* iptmy ; E_Float* iptmz; E_Float* ipte; 
           E_Float* iptve = NULL ; E_Float* iptvm = NULL ;
           if (t == NULL )
           {
                PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable Density is missing. Not written.\n");
                  return ; 
           }
           else
           {  iptro    = K_PYTREE::getValueAF(t, hook);
//            cerr <<" pointeur Density récupéré "<<endl ;
           }
           t  = K_PYTREE::getNodeFromName1(bcdata, "MomentumX");
           if (t == NULL )
           {
                PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable MomentumX is missing. Not written.\n");
                  return ; 
           }
           else 
           { iptmx    = K_PYTREE::getValueAF(t, hook);
//             cerr <<" pointeur MomentumX récupéré "<<endl ;
           }
           t  = K_PYTREE::getNodeFromName1(bcdata, "MomentumY");
           if (t == NULL )
           {
                PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable MomentumY is missing. Not written.\n");
                  return ; 
           }
           else 
           { iptmy    = K_PYTREE::getValueAF(t, hook);
//             cerr <<" pointeur MomentumY récupéré "<<endl ;
           }
           t  = K_PYTREE::getNodeFromName1(bcdata, "MomentumZ");
           if (t == NULL )
           {
                PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable MomentumZ is missing. Not written.\n");
                  return ; 
           }
           else 
           { iptmz    = K_PYTREE::getValueAF(t, hook);
//           cerr <<" pointeur MomentumZ récupéré "<<endl ;
           }
           t  = K_PYTREE::getNodeFromName1(bcdata, "EnergyStagnationDensity");
           if (t == NULL )
           {
                PyErr_SetString(PyExc_TypeError,
                    "convertPyTree2FFD : variable EnergyStagnationDensity is missing. Not written.\n");
                  return ; 
           }
           else
           {  ipte     = K_PYTREE::getValueAF(t, hook);
//           cerr <<" pointeur EnergyStagnationDensity récupéré "<<endl ;
           }
           if( kodnst > 0){
               t  = K_PYTREE::getNodeFromName1(bcdata, "ViscosityEddy");
               if (t == NULL )
                          {
                              //PyErr_SetString(PyExc_TypeError,
                              //    "convertPyTree2FFD : variable ViscosityEddy is missing. Not written.\n");
                              //  return ; 
                              iptve = NULL ;
                          }
                else iptve    = K_PYTREE::getValueAF(t, hook);
                t  = K_PYTREE::getNodeFromName1(bcdata, "ViscosityMolecular");
                if (t == NULL )
                          {
                              // PyErr_SetString(PyExc_TypeError,
                              //     "convertPyTree2FFD : variable ViscosityMolecular is missing. Not written.\n");
                              //   return ; 
                              iptvm = NULL ;
                          }
                else iptvm    = K_PYTREE::getValueAF(t, hook);
             }
//           cerr << " Debut transfert variables BC"<<endl;
           for (E_Int i = 0; i < PLSize; i++)
             {
             *(Var_l++) = *(iptro++) ;
             *(Var_l++) = *(iptmx++) ;
             *(Var_l++) = *(iptmy++) ;
             *(Var_l++) = *(iptmz++) ;
             *(Var_l++) = *(ipte++) ;
                if( kodnst > 0 ){
                if ( iptve == NULL){
                if( i == 0)printf("Pas de Mut dans BCData ! Mut/Mu est mis à 1.\n");
                *(Var_l++) = 1. ;
                }
                else {
                E_Float Mut        = *(iptve++) ;
                E_Float Mu         = *(iptvm++) ;
                *(Var_l++) = Mut/Mu ;
                }
                }
             ielmbc += 1 ;
             }
             }
            } //pas de data
//            else //pas de data
//           cerr <<" getVarBC : pas de data = " <<endl ;
          }
      } // loop ibc (ii)

   } // loop bc (i)
//   cerr << "sortie de getVarBC"<<endl ;
}
void K_CONVERTER::copyPE2Elmtmtch(E_Int* p, E_Int *ielmtmtch1, E_Int *ielmtmtch2, E_Int nmtch, E_Int *nBC )
{
//      cerr<<"ielmtmtch1 ="<<ielmtmtch1<<endl ; 
//      cerr<<"ielmtmtch2 ="<<ielmtmtch2<<"  nmtch = "<<nmtch  <<endl ; 
      for ( E_Int i = 1 ;  i < nmtch ; i++ ){
//      cerr<<" ielmtmtch1[i] ="<<ielmtmtch1[i]<<" ielmtmtch2[i] ="<<ielmtmtch2[i]<<"  i = "<< i <<endl ; 
      if(  ielmtmtch2[i] == 0 ){  
                                  ielmtmtch2[i] = -9999;
                                  *nBC += 1 ;
                               }
      }
}
//void K_CONVERTER::copyPE2Elmtmtch(E_Int* p, E_Int *ielmtmtch1, E_Int *ielmtmtch2, E_Int nmtch, E_Int *nBC )
//{
//      cerr<<"p ="<<p<<"  nmtch = "<<nmtch <<endl ; 
//      ielmtmtch1[0] =  *(p++) ;
//      ielmtmtch2[0] =  *(p++) ;
//      for ( E_Int i = 1 ;  i < nmtch ; i++ ){
//      ielmtmtch1[i] =  *(p++);
//      ielmtmtch2[i] =  *(p++);
//      cerr<<" ielmtmtch1[i] ="<<ielmtmtch1[i]<<" ielmtmtch2[i] ="<<ielmtmtch2[i]<<"  i = "<< i <<endl ; 
//      if(  ielmtmtch2[i] == 0 ){  
//                                  ielmtmtch2[i] = -9999;
//                                  *nBC += 1 ;
//                               }
//      }
//}

void K_CONVERTER::splitElementConnectivity(PyObject* zone, E_Int* npoint,E_Int* nodmtch,E_Int* kpoinmtch)
{
    vector<PyArrayObject*> hook;
    E_Int NElCon, s1;
    PyObject* NF = K_PYTREE::getNodeFromName1(zone, "NFaceElements");
    PyObject* EC = K_PYTREE::getNodeFromName1(NF,  "ElementConnectivity");
    E_Int* p =  K_PYTREE::getValueAI(EC, NElCon, s1, hook);
//    cerr<<"p ="<<p<<" NElCon  = "<<NElCon <<endl ; 
    E_Int j ;
    E_Int i = 0 ; E_Int l = 0 ;
    kpoinmtch[0] = 1  ;
    do {
       i++ ;
       kpoinmtch[i] = kpoinmtch[i-1] + *(p++) ;
       l++;
       for ( j = kpoinmtch[i-1] ; j <   kpoinmtch[i] ; j++ )
       {
       nodmtch[j]= *(p++) ;
       l++ ;
       }
       *npoint = kpoinmtch[i] ;
    } while ( l < NElCon ) ;
}
//==============================================================================
// Recherche par nom d'un seul niveau, retourne une liste de noeuds.
// IN: o: objet representant un noeud de pyTree
// IN: type: le type du noeud
// OUT: out: la liste des noeuds trouves
//==============================================================================
//void K_PYTREE::getNodesFromTypeDB(PyObject* o, const char* type,
//                                 vector<PyObject*>& out)
//{
//  PyObject* childrens = PyList_GetItem(o, 2);
//  E_Int n = PyList_Size(childrens);
//  PyObject *l; PyObject *node; char* str;
//  for (E_Int i = 0; i < n; i++)
//  {
//    l = PyList_GetItem(childrens, i);
//    node = PyList_GetItem(l, 3);
//    if (PyString_Check(node)) str = PyString_AsString(node); // type_bc
//#if PY_VERSION_HEX >= 0x03000000
//       else if (PyUnicode_Check(node)) str = PyBytes_AsString(PyUnicode_AsUTF8String(node));
//#endif
//    if (K_STRING::cmp(str, type) == 0) out.push_back(l);
//  }
//}

