/*    
    Copyright 2013-2025 Onera.

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
 
# include "rigidMotion.h"
using namespace K_FLD;
using namespace std;
using namespace K_CONST;

extern "C"
{
  void evalrotfor_(
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&, const E_Float&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&, const E_Float&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&, const E_Float&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&, const E_Float&,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&,
    const E_Float&, const E_Int&, const E_Float*, const E_Float*,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&,
    const E_Float&, const E_Int&, const E_Float*, const E_Float*,
    const E_Float&, const E_Float&, const E_Float&,
    const E_Int&,
    const E_Float&, const E_Int&, const E_Float*, const E_Float*,
    E_Float&, E_Float&,E_Float&,E_Float&,
    E_Float*, E_Float*, E_Float*, E_Float*, E_Float*);
}

//=============================================================================
/* Compute rotor motion matrix, coordinates of the origin of the relative frame
 in the abs frame and the coords of the center of rotation in the relative frame
*/
//=============================================================================
PyObject* K_RIGIDMOTION::_computeRotorMotionInfo(PyObject* self, PyObject* args)
{
  E_Float time, psi0, psi0_b, alp0, rot_omg, del0, bet0, tet0, pre_lag_ang, pre_con_ang;
  PyObject *transl_speed, *alp_pnt0, *alp_vct, *rot_pnt0, *rot_vct,
    *del_pnt0, *del_vct, *delc0, *dels0,
    *bet_pnt0, *bet_vct, *betc0, *bets0,
    *tet_pnt0, *tet_vct, *tetc0, *tets0,
    *span_vct, *pre_lag_pnt0, *pre_lag_vct, *pre_con_pnt0, *pre_con_vct;
  if (!PYPARSETUPLE_(args, 
                     R_ O_ RR_ OO_ R_ OO_ R_ OO_ R_ OOOO_ R_ OOOO_ R_ OOO_ R_ OO_ R_ OO_,
                     &time, &transl_speed, &psi0, &psi0_b, 
                     &alp_pnt0, &alp_vct, &alp0,
                     &rot_pnt0, &rot_vct, &rot_omg,                     
                     &del_pnt0, &del_vct, &del0, &delc0, &dels0,
                     &bet_pnt0, &bet_vct, &bet0, &betc0, &bets0,
                     &tet_pnt0, &tet_vct, &tet0, &tetc0, &tets0,
                     &span_vct,
                     &pre_lag_ang, &pre_lag_pnt0, &pre_lag_vct,
                     &pre_con_ang, &pre_con_pnt0, &pre_con_vct))
    return NULL;

  E_Int blade_span_axis=0, axis0=0, axis1=0, axis2=0, axis3=0, axis4=0, axis5=0, axis6=0;
  E_Float transl[3]; E_Float alp_pnt[3];E_Float rotor_pnt[3];
  E_Float pre_lag_pnt[3]; E_Float pre_con_pnt[3];
  E_Float del_pnt[3]; E_Float bet_pnt[3]; E_Float tet_pnt[3];
  
  //harmonics for lead-lag, pitching and flapping
  E_Int nhdel = PyList_Size(delc0);
  FldArrayF delc(nhdel); FldArrayF dels(nhdel);
  E_Float* ptdelc = delc.begin(); E_Float* ptdels = dels.begin();
  for (E_Int nov=0; nov < nhdel; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(delc0,nov);
    ptdelc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(dels0,nov);
    ptdels[nov] = PyFloat_AsDouble(tpl0);
  }
  E_Int nhbet = PyList_Size(betc0);
  FldArrayF betc(nhbet); FldArrayF bets(nhbet);
  E_Float* ptbetc = betc.begin(); E_Float* ptbets = bets.begin();
  for (E_Int nov=0; nov < nhbet; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(betc0,nov);
    ptbetc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(bets0,nov);
    ptbets[nov] = PyFloat_AsDouble(tpl0);
  }
  E_Int nhtet = PyList_Size(tetc0);
  FldArrayF tetc(nhtet); FldArrayF tets(nhtet);
  E_Float* pttetc = tetc.begin(); E_Float* pttets = tets.begin();
  for (E_Int nov=0; nov<nhtet; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(tetc0,nov);
    pttetc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(tets0,nov);
    pttets[nov] = PyFloat_AsDouble(tpl0);
  }
  
  for (E_Int nov = 0; nov < 3; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(transl_speed,nov);
    E_Float val = PyFloat_AsDouble(tpl0);
    transl[nov] = val;

    tpl0 = PyList_GetItem(alp_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    alp_pnt[nov] = PyFloat_AsDouble(tpl0);
    
    tpl0 = PyList_GetItem(rot_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    rotor_pnt[nov] = val;

    tpl0 = PyList_GetItem(pre_lag_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    pre_lag_pnt[nov] = val;
    
    tpl0 =  PyList_GetItem(pre_con_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    pre_con_pnt[nov] = val;

    tpl0 = PyList_GetItem(del_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    del_pnt[nov] = val;

    tpl0 = PyList_GetItem(bet_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    bet_pnt[nov] = val;
    
    tpl0 =  PyList_GetItem(tet_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    tet_pnt[nov] = val;
        
    tpl0 =  PyList_GetItem(span_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) blade_span_axis = nov+1;

    tpl0 =  PyList_GetItem(alp_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis0 = nov+1;

    tpl0 =  PyList_GetItem(rot_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis1 = nov+1;
    
    tpl0 =  PyList_GetItem(pre_lag_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis2 = nov+1;

    tpl0 =  PyList_GetItem(pre_con_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis3 = nov+1;

    tpl0 =  PyList_GetItem(del_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis4 = nov+1;
    
    tpl0 =  PyList_GetItem(bet_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis5= nov+1;

    tpl0 =  PyList_GetItem(tet_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis6= nov+1;
  }

  FldArrayF rotMat(3,3);//matrice du mouvement
  FldArrayF r0(3);
  FldArrayF x0(3);
  FldArrayF s0(3);
  FldArrayF omega(3);
  E_Float psideg, deldeg, betdeg, tetdeg;
  evalrotfor_(time, psi0, psi0_b,
              blade_span_axis,
              transl[0], transl[1], transl[2],
              alp_pnt[0], alp_pnt[1], alp_pnt[2], axis0, alp0,
              rotor_pnt[0], rotor_pnt[1], rotor_pnt[2], axis1, rot_omg,
              pre_lag_pnt[0], pre_lag_pnt[1], pre_lag_pnt[2],
              axis2, pre_lag_ang,
              pre_con_pnt[0], pre_con_pnt[1], pre_con_pnt[2],
              axis3, pre_con_ang,
              del_pnt[0], del_pnt[1], del_pnt[2],
              axis4, del0, E_Int(nhdel), delc.begin(), dels.begin(),
              bet_pnt[0], bet_pnt[1], bet_pnt[2],
              axis5, bet0, E_Int(nhbet), betc.begin(), bets.begin(),
              tet_pnt[0], tet_pnt[1], tet_pnt[2],
              axis6, tet0, E_Int(nhtet), tetc.begin(), tets.begin(),
              psideg, deldeg, betdeg, tetdeg,
              r0.begin(), rotMat.begin(), x0.begin(),
              s0.begin(), omega.begin());

  // printf(" MATRICE : \n");
  // for (E_Int no = 0; no < 3; no++)
  //   printf(" %f %f %f \n", rotMat(no,1), rotMat(no,2), rotMat(no,3));
 
  // printf("Deplacement : \n");
  // printf(" %g %g %g \n", r0[0], r0[1], r0[2]);
  // printf(" psi=%f, pitch=%f, flap=%f, lag=%f\n",
  //        psideg, tetdeg, betdeg, deldeg); 
  PyObject* l = PyList_New(0);
  PyObject* tpl = K_NUMPY::buildNumpyArray(r0);
  PyList_Append(l,tpl); Py_DECREF(tpl);
  tpl = K_NUMPY::buildNumpyArray(x0);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  tpl = K_NUMPY::buildNumpyArray(rotMat,1);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  tpl = K_NUMPY::buildNumpyArray(s0);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  tpl = K_NUMPY::buildNumpyArray(omega);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  
  return l;
}

//=============================================================================
/* Move grid and compute grid velocity according to rotor_motion */
//=============================================================================
PyObject* K_RIGIDMOTION::_computeRotorMotionZ(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  E_Float time, psi0, psi0_b, alp0, rot_omg, del0, bet0, tet0, pre_lag_ang, pre_con_ang;
  PyObject *zone, *sxo, *syo, *szo, *transl_speed, *alp_pnt0, *alp_vct, *rot_pnt0, *rot_vct,
    *del_pnt0, *del_vct, *delc0, *dels0,
    *bet_pnt0, *bet_vct, *betc0, *bets0,
    *tet_pnt0, *tet_vct, *tetc0, *tets0,
    *span_vct, *pre_lag_pnt0, *pre_lag_vct, *pre_con_pnt0, *pre_con_vct;
  if (!PYPARSETUPLE_(args,
                      OOOO_ R_ O_ RR_ OO_ R_ OO_ R_ OO_ R_ OOOO_ R_ OOOO_ R_ OOO_ R_ OO_ R_ OO_ SSS_,
                     &zone, &sxo, &syo, &szo,
                     &time, &transl_speed, &psi0, &psi0_b, 
                     &alp_pnt0, &alp_vct, &alp0,
                     &rot_pnt0, &rot_vct, &rot_omg,                     
                     &del_pnt0, &del_vct, &del0, &delc0, &dels0,
                     &bet_pnt0, &bet_vct, &bet0, &betc0, &bets0,
                     &tet_pnt0, &tet_vct, &tet0, &tetc0, &tets0,
                     &span_vct,
                     &pre_lag_ang, &pre_lag_pnt0, &pre_lag_vct,
                     &pre_con_ang, &pre_con_pnt0, &pre_con_vct,
                     &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
    return NULL;
                 
  vector<PyArrayObject*> hookz;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(zone, 1, 0, varString, fields, locs, 
                                    im, jm, km, 
                                    cn, cnSize, cnNfld, eltType, hookz, 
                                    GridCoordinates, 
                                    FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDZ(hookz, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "rotorMotionZ: cannot find coordinates in zone.");
    return NULL;
  }    
  E_Int npts;
  if (res == 1) npts = im*jm*km;
  else npts = im;

  E_Int blade_span_axis=0, axis0=0, axis1=0, axis2=0, axis3=0, axis4=0, axis5=0, axis6=0;
  E_Float transl[3]; E_Float alp_pnt[3];E_Float rotor_pnt[3];
  E_Float pre_lag_pnt[3]; E_Float pre_con_pnt[3];
  E_Float del_pnt[3]; E_Float bet_pnt[3]; E_Float tet_pnt[3];
  
  //harmonics for lead-lag, pitching and flapping
  E_Int nhdel = PyList_Size(delc0);
  FldArrayF delc(nhdel); FldArrayF dels(nhdel);
  E_Float* ptdelc = delc.begin(); E_Float* ptdels = dels.begin();
  for (E_Int nov=0; nov<nhdel; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(delc0,nov);
    ptdelc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(dels0,nov);
    ptdels[nov] = PyFloat_AsDouble(tpl0);
  }
  E_Int nhbet = PyList_Size(betc0);
  FldArrayF betc(nhbet); FldArrayF bets(nhbet);
  E_Float* ptbetc = betc.begin(); E_Float* ptbets = bets.begin();
  for (E_Int nov=0; nov<nhbet; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(betc0,nov);
    ptbetc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(bets0,nov);
    ptbets[nov] = PyFloat_AsDouble(tpl0);
  }
  E_Int nhtet = PyList_Size(tetc0);
  FldArrayF tetc(nhtet); FldArrayF tets(nhtet);
  E_Float* pttetc = tetc.begin(); E_Float* pttets = tets.begin();
  for (E_Int nov=0; nov<nhtet; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(tetc0,nov);
    pttetc[nov] = PyFloat_AsDouble(tpl0);
    tpl0 = PyList_GetItem(tets0,nov);
    pttets[nov] = PyFloat_AsDouble(tpl0);
  }
  
  for (E_Int nov = 0; nov < 3; nov++)
  {
    PyObject* tpl0 = PyList_GetItem(transl_speed,nov);
    E_Float val = PyFloat_AsDouble(tpl0);
    transl[nov] = val;

    tpl0 = PyList_GetItem(alp_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    alp_pnt[nov] = PyFloat_AsDouble(tpl0);
    
    tpl0 = PyList_GetItem(rot_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    rotor_pnt[nov] = val;

    tpl0 = PyList_GetItem(pre_lag_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    pre_lag_pnt[nov] = val;
    
    tpl0 =  PyList_GetItem(pre_con_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    pre_con_pnt[nov] = val;

    tpl0 = PyList_GetItem(del_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    del_pnt[nov] = val;

    tpl0 = PyList_GetItem(bet_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    bet_pnt[nov] = val;
    
    tpl0 =  PyList_GetItem(tet_pnt0,nov);
    val = PyFloat_AsDouble(tpl0);
    tet_pnt[nov] = val;
        
    tpl0 =  PyList_GetItem(span_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) blade_span_axis = nov+1;

    tpl0 =  PyList_GetItem(alp_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis0 = nov+1;

    tpl0 =  PyList_GetItem(rot_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis1 = nov+1;
    
    tpl0 =  PyList_GetItem(pre_lag_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis2 = nov+1;

    tpl0 =  PyList_GetItem(pre_con_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis3 = nov+1;

    tpl0 =  PyList_GetItem(del_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis4 = nov+1;
    
    tpl0 =  PyList_GetItem(bet_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis5 = nov+1;

    tpl0 =  PyList_GetItem(tet_vct,nov);
    val = PyFloat_AsDouble(tpl0);
    if (K_FUNC::fEqualZero(val) == false) axis6 = nov+1;
  }

  FldArrayF rotMat(3,3); //matrice du mouvement
  FldArrayF r0(3);
  FldArrayF x0(3);
  FldArrayF s0(3);
  FldArrayF omega(3);
  E_Float psideg, deldeg, betdeg, tetdeg;
  evalrotfor_(time, psi0, psi0_b,
              blade_span_axis,
              transl[0], transl[1], transl[2],
              alp_pnt[0], alp_pnt[1], alp_pnt[2], axis0, alp0,
              rotor_pnt[0], rotor_pnt[1], rotor_pnt[2], axis1, rot_omg,
              pre_lag_pnt[0], pre_lag_pnt[1], pre_lag_pnt[2],
              axis2, pre_lag_ang,
              pre_con_pnt[0], pre_con_pnt[1], pre_con_pnt[2],
              axis3, pre_con_ang,
              del_pnt[0], del_pnt[1], del_pnt[2],
              axis4, del0, E_Int(nhdel), delc.begin(), dels.begin(),
              bet_pnt[0], bet_pnt[1], bet_pnt[2],
              axis5, bet0, E_Int(nhbet), betc.begin(), bets.begin(),
              tet_pnt[0], tet_pnt[1], tet_pnt[2],
              axis6, tet0, E_Int(nhtet), tetc.begin(), tets.begin(),
              psideg, deldeg, betdeg, tetdeg,
              r0.begin(), rotMat.begin(), x0.begin(),
              s0.begin(), omega.begin());

  vector<PyArrayObject*> hook;
  //char* zoneName = K_PYTREE::getNodeName(zone, hook);
  // printf(" MATRICE RMZ : \n");
  // for (E_Int no = 0; no < 3; no++)
  //   printf(" %f %f %f \n", rotMat(no,1), rotMat(no,2), rotMat(no,3));
 
  // printf("Deplacement : \n");
  // printf(" %g %g %g \n", r0[0], r0[1], r0[2]);
  // printf(" Zone %s : psi=%f, pitch=%f, flap=%f, lag=%f\n", zoneName,
  //        psideg, tetdeg, betdeg, deldeg); 
  RELEASEHOOK(hook);
  
  // move zone
  E_Float* xt = fields[posx];
  E_Float* yt = fields[posy];
  E_Float* zt = fields[posz];

  // coordinates of the origin of the relative frame in the abs frame
  E_Float xa = r0[0]; E_Float ya = r0[1]; E_Float za = r0[2];
  // coordinates of the center of rotation in the relative frame
  E_Float xr = x0[0]; E_Float yr = x0[1]; E_Float zr = x0[2];
  //rotation matrix 
  E_Float r11 = rotMat(0,1); E_Float r12 = rotMat(0,2); E_Float r13 = rotMat(0,3);
  E_Float r21 = rotMat(1,1); E_Float r22 = rotMat(1,2); E_Float r23 = rotMat(1,3);
  E_Float r31 = rotMat(2,1); E_Float r32 = rotMat(2,2); E_Float r33 = rotMat(2,3);

#pragma omp parallel default(shared) 
  {
    E_Float x, y, z;
  #pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      x = xt[ind]; y = yt[ind]; z = zt[ind];
      xt[ind] = r11*(x-xr) + r12*(y-yr) + r13*(z-zr) + xa;
      yt[ind] = r21*(x-xr) + r22*(y-yr) + r23*(z-zr) + ya;
      zt[ind] = r31*(x-xr) + r32*(y-yr) + r33*(z-zr) + za;
    }
  }
  //check grid velocity numpys
  E_Int size; E_Int nfld;
  E_Float* sx;
  K_NUMPY::getFromNumpyArray(sxo, sx, size, nfld, true);
  E_Float* sy;
  K_NUMPY::getFromNumpyArray(syo, sy, size, nfld, true);
  E_Float* sz;
  K_NUMPY::getFromNumpyArray(szo, sz, size, nfld, true);
  size = size*nfld;

  E_Float s01 = s0[0]; E_Float s02 = s0[1]; E_Float s03 = s0[2];
  E_Float omg1 = omega[0];
  E_Float omg2 = omega[1];
  E_Float omg3 = omega[2];


  E_Float tx = s01 - (omg2 * zr - omg3 * yr);
  E_Float ty = s02 - (omg3 * xr - omg1 * zr);
  E_Float tz = s03 - (omg1 * yr - omg2 * xr);

#pragma omp parallel
  {
#pragma omp for
    for (E_Int i = 0; i < size; i++)
    {
        sx[i] = tx + (omg2 * zt[i] - omg3 * yt[i]);
        sy[i] = ty + (omg3 * xt[i] - omg1 * zt[i]); 
        sz[i] = tz + (omg1 * yt[i] - omg2 * xt[i]);
    }
         // faux
  }
  
  if (res == 2) delete [] eltType;
  RELEASESHAREDZ(hookz, varString, (char*)NULL);
  Py_INCREF(Py_None); 
  return Py_None;
}
