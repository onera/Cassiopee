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

// Initialization by vortices

# include "initiator.h"

using namespace K_FLD;
using namespace std;

/*
  Lamb-Oseen vortex initialization
*/
void k6lamb(E_Float x0, E_Float y0, E_Float Gamma, E_Float MInf,
            E_Int npts, 
            const E_Float* xc, const E_Float* yc, const E_Float* zc,
            E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  /* Constants */
  const E_Int nd = 50000; /* discretisation */
  
  //E_Float pc, roc;
  E_Float ro0, a0;
  E_Float p0, h0, S0;
  E_Int ind;
  E_Float dr, r, aa;
  E_Float Sinf, roinf, uinf, pinf;
  E_Int i, m;
  E_Float rmax, vr, vrdr;
  E_Float cos_teta, sin_teta;
  E_Float va, ss, ro, pp;
  E_Float gam;
  E_Float pi;

  E_Float* S = new E_Float [ 2*nd ];

  pi = acos(-1.0);

  //pc = 1.0 / 1.4; /* pression critique */
  //roc = 1.0;

  ro0 = 1.0; /* état d'adimensionnement = état d'arrêt */
  p0 = 1.0 / 1.4;
  a0 = 1.0;
  h0 = 1.0 / 0.4;
  S0 = p0 / pow(ro0, 1.4);

  gam = Gamma * a0; /* intensité du tourbillon */

  rmax = 0.0;
  for (ind = 0; ind < npts; ind++)
  {
    rmax = fmax(rmax, sqrt((xc[ind] - x0) * (xc[ind] - x0)
                         + (yc[ind] - y0) * (yc[ind] - y0)));
  }

  dr = rmax / (1.0 * nd);
  S[0] = S0;

  for (i = 2; i <= 2 * nd; i++)
  {
    r = (i - 1) * dr;

    vr = gam * (1.0 - exp(-r * r)) / (2.0 * pi * r);
    vrdr = gam * (1.0 - exp(-(r + dr) * (r + dr))) / (2.0 * pi * (r + dr));
    aa = (r + dr) * vrdr - r * vr;

    aa = -1.4 * vr * aa / r;
    aa = aa / (h0 - 0.5 * vr * vr);
    S[i - 1] = log(S[i - 2]) + aa;
    S[i - 1] = exp(S[i - 1]);
  }

  Sinf = S[2 * nd - 1];
  roinf = 0.4 * h0 / (1.4 * Sinf);
  roinf = pow(roinf, 1.0 / 0.4);
  pinf = roinf * (0.4 / 1.4) * h0;
  uinf = MInf * sqrt(1.4 * pinf / roinf);

  /* pression et densité en chaque point */
  for (ind = 0; ind < npts; ind++)
  {
    r = sqrt((xc[ind] - x0) * (xc[ind] - x0)
           + (yc[ind] - y0) * (yc[ind] - y0));
    if (r <= 1.e-8)
    {
      cos_teta = 0.0;
      sin_teta = 0.0;
      va = 0.0;
    }
    else
    {
      cos_teta = (yc[ind] - y0) / r;
      sin_teta = -(xc[ind] - x0) / r;
      va = gam * (1.0 - exp(-r * r)) / (2.0 * pi * r);
    }

    m = (E_Int)(r * nd / rmax) + 1;
    ss = S[m - 1];
    ro = 0.4 * (h0 - 0.5 * va * va) / (1.4 * ss);
    ro = pow(ro, 1.0 / 0.4);
    pp = ro * 0.4 * (h0 - 0.5 * va * va) / 1.4;

    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = pp / 0.4 + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
  delete [] S;
}

/*
  Initialization of Visbal vortex in field with constant density
*/
void k6visbal(E_Float x0, E_Float y0, E_Float Gamma, E_Float MInf,
              E_Int npts,
              const E_Float* xc, const E_Float* yc, const E_Float* zc,
              E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  /* Constants */
  const E_Int smax = 12000;
  E_Float* S = new E_Float [smax];

  /* Local variables */
  //E_Float ro0, u0, roinf;
  E_Float ainf, hinf, ar, ardr;
  E_Float p0;
  E_Int ind;
  E_Float dr, r;
  E_Float Sinf, uinf, pinf;
  E_Int i, m;
  E_Float rmax, vr, vrdr;
  E_Float cos_teta, sin_teta;
  E_Float va, ss, ro, pp;
  E_Float gam;

  //roinf = 1.0; /* etat d'adimensionnement : etat infini */
  pinf = 1.0 / 1.4;
  ainf = 1.0;
  hinf = 1.0 / 0.4; /* constante */
  uinf = MInf;
  gam = Gamma * ainf; /* intensité du tourbillon */

  rmax = 0.0;
  for (ind = 0; ind < npts; ind++)
  {
    rmax = fmax(rmax, sqrt((xc[ind] - x0) * (xc[ind] - x0)
                         + (yc[ind] - y0) * (yc[ind] - y0)));
  }

  dr = rmax / (1.0 * 6000);
  r = dr;
  vr = gam / (2.0 * 3.14) * exp(0.5) * r * exp(-(r * r) / 2.0);
  S[0] = exp((1.4 / 0.4) * (vr * vr / (r * (hinf - 0.5 * vr * vr))) * dr / 2.0);

  for (i = 2; i <= smax; i++)
  {
    r = (i - 1) * dr;
    vr = gam / (2.0 * 3.14) * exp(0.5) * r * exp(-(r * r) / 2.0);
    vrdr = gam / (2.0 * 3.14) * exp(0.5) * (r + dr) * exp(-((r + dr) * (r + dr)) / 2.0);
    ar = vr * vr / (r * (hinf - 0.5 * vr * vr));
    ardr = vrdr * vrdr / ((r + dr) * (hinf - 0.5 * vrdr * vrdr));
    S[i - 1] = log(S[i - 2]) + ((1.4 / 0.4) * (ardr + ar) / 2.0) * dr;
    S[i - 1] = exp(S[i - 1]);
  }

  Sinf = S[smax - 1];
  p0 = pinf / Sinf;
  //ro0 = (1.4 * p0) / (0.4 * hinf);
  //u0 = MInf;

  /* Pression et densité en chaque point */
  for (ind = 0; ind < npts; ind++)
  {
    r = sqrt((xc[ind] - x0) * (xc[ind] - x0)
           + (yc[ind] - y0) * (yc[ind] - y0));
    if (r <= 1.e-8)
    {
      cos_teta = 0.0;
      sin_teta = 0.0;
      va = 0.0;
    }
    else
    {
      cos_teta = (yc[ind] - y0) / r;
      sin_teta = -(xc[ind] - x0) / r;
      va = gam / (2.0 * 3.14) * exp(0.5) * r * exp(-r * r / 2.0);
    }

    m = (E_Int)(r * 1.0 * 6000 / rmax) + 1;
    ss = S[m - 1];
    pp = ss * p0;
    ro = (1.4 * pp) / (0.4 * (hinf - 0.5 * va * va));

    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = pp / 0.4 + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
  delete [] S;
}

/*
  Initialization of Scully vortex in field
  Adimensionalization with ro_inf = 1 and u_inf = 1
*/
void k6scully(E_Float x0, E_Float y0, E_Float Gamma, E_Float a, E_Float MInf,
              E_Int npts,
              const E_Float* xc, const E_Float* yc, const E_Float* zc,
              E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  /* Local variables */
  E_Int ind;
  E_Float r, p, pa, ptmp, roa, T;
  E_Float roinf, uinf, pinf;
  E_Float cos_teta, sin_teta;
  E_Float va, ro, pi;

  pi = 4.0 * atan(1.0);

  /* Adimensionalization state */
  roinf = 1.0;
  uinf = 1.0;
  pinf = 1.0 / (1.4 * MInf * MInf);

  /* Initialization */
  ptmp = (0.4 / 1.4) * (Gamma * Gamma) / (8.0 * pi * pi * pinf);

  pa = pinf * pow((1.0 - ptmp / (2.0 * a * a)), (1.4 / 0.4));
  roa = pow((pa / pinf), (1.0 / 1.4)) * roinf;
  T = pa / (roa * 287.0);

  for (ind = 0; ind < npts; ind++)
  {
    r = sqrt((xc[ind] - x0) * (xc[ind] - x0)
           + (yc[ind] - y0) * (yc[ind] - y0));

    if (r <= 1.e-12)
    {
      cos_teta = 0.0;
      sin_teta = 0.0;
      va = 0.0;
    }
    else
    {
      cos_teta = (yc[ind] - y0) / r;
      sin_teta = -(xc[ind] - x0) / r;
      va = Gamma / (2.0 * pi) * r / (r * r + a * a);
    }

    p = pinf * pow((1.0 - ptmp / (r * r + a * a)), (1.4 / 0.4));

    if (r >= a)
    {
      ro = pow((1.0 - ptmp / (r * r + a * a)), (1.0 / 0.4)) * roinf;
    }
    else
    {
      ro = p / (287.0 * T);
    }

    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = p / 0.4 + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
}

/*
  Initialization of Scully vortex (variant 2) in field
  Adimensionalization with ro_inf = 1 and u_inf = 1
*/
void k6scully2(E_Float x0, E_Float y0, E_Float Gamma, E_Float a, E_Float MInf,
               E_Int npts,
               const E_Float* xc, const E_Float* yc, const E_Float* zc,
               E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  /* Local variables */
  E_Int ind;
  E_Float r, p, ptmp;
  E_Float roinf, uinf, pinf;
  E_Float cos_teta, sin_teta;
  E_Float va, ro, pi;

  pi = 4.0 * atan(1.0);

  /* Adimensionalization state */
  roinf = 1.0;
  uinf = 1.0;
  pinf = 1.0 / (1.4 * MInf * MInf);

  /* Precompute ptmp */
  ptmp = (0.4 / 1.4) * (Gamma * Gamma) / (8.0 * pi * pi * pinf);

  for (ind = 0; ind < npts; ind++)
  {
    r = sqrt((xc[ind] - x0) * (xc[ind] - x0) +
             (yc[ind] - y0) * (yc[ind] - y0));

    if (r <= 1.e-12)
    {
      cos_teta = 0.0;
      sin_teta = 0.0;
      va = 0.0;
    }
    else
    {
      cos_teta = (yc[ind] - y0) / r;
      sin_teta = -(xc[ind] - x0) / r;
      va = Gamma / (2.0 * pi) * r / (r * r + a * a);
    }

    p = pinf * pow((1.0 - ptmp / (r * r + a * a)), (1.4 / 0.4));
    ro = pow((1.0 - ptmp / (r * r + a * a)), (1.0 / 0.4)) * roinf;

    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = p / 0.4 + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
}

/*
  Initialization of Yee vortex in field
  Adimensionalization: ro_inf = 1, u_inf = Minf
*/
void k6yee(E_Float x0, E_Float y0, E_Float Gamma, E_Float Minf,
           E_Int npts, 
           const E_Float* xc, const E_Float* yc, const E_Float* zc,
           E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  /* Local variables */
  E_Int ind;
  E_Float r2;
  //E_Float u0;
  E_Float roinf, uinf, pinf;
  E_Float ro0, ainf, p0, t0;
  E_Float gam, gma, pi, rgp;
  E_Float cos_teta, sin_teta;
  E_Float va, ro, pp, ta;

  pi = 3.14159265358979323846;
  gma = 1.4;
  rgp = 287.53;

  /* Adimensionalization state */
  roinf = 1.0;          /* far-field density */
  pinf  = 1.0 / gma;    /* far-field pressure */
  ainf  = 1.0;          /* far-field sound speed */
  uinf  = Minf;         /* far-field velocity */
  gam   = Gamma * ainf; /* vortex intensity */

  p0 = pinf;
  ro0 = roinf;
  //u0 = Minf;
  t0 = p0 / (ro0 * rgp);

  /* Loop over points */
  for (ind = 0; ind < npts; ind++)
  {
    r2 = (xc[ind] - x0) * (xc[ind] - x0) +
         (yc[ind] - y0) * (yc[ind] - y0);

    cos_teta = -(yc[ind] - y0);
    sin_teta = +(xc[ind] - x0);

    va = (gam / (2.0 * pi)) * exp((1.0 - r2) / 2.0);
    ta = 1.0 - ((gma - 1.0) / (2.0 * gma * rgp * t0)) * va * va;

    pp = pow(ta, gma / (gma - 1.0)) * p0;
    ro = pow(ta, 1.0 / (gma - 1.0)) * ro0;

    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = pp / (gma - 1.0) + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
}

/*
  Initialization of Wissocq vortex
  Dimensional: ro_inf = 1.1765, p_inf = 101320 Pa
*/
void k6wissocq(E_Float x0, E_Float y0, E_Float Gamma, E_Float MInf,
               E_Int npts,
               const E_Float* xc, const E_Float* yc, const E_Float* zc,
               E_Float* u1, E_Float* u2, E_Float* u3, E_Float* u4, E_Float* u5)
{
  E_Int ind;
  E_Float r, rc, rc2;
  E_Float roinf, uinf, pinf, tinf;
  E_Float cos_teta, sin_teta;
  E_Float va, ro, pp;
  E_Float gam, c0, coef;

  /* Reference state */
  roinf = 1.1765;               /* kg/m^3 */
  pinf = 101320.0;             /* Pa */
  tinf = pinf / (roinf * 287.053); /* K */

  c0 = sqrt(1.4 * 287.053 * tinf); /* speed of sound */
  gam = Gamma * c0;                /* vortex intensity */
  coef = 0.5 * (gam * gam) / (c0 * c0);
  uinf = MInf * c0;
  rc = 0.1; /* vortex core radius */
  rc2 = rc * rc;

  for (ind = 0; ind < npts; ind++)
  {
    r = sqrt((xc[ind] - x0) * (xc[ind] - x0) +
             (yc[ind] - y0) * (yc[ind] - y0));

    cos_teta = -(yc[ind] - y0) / rc;
    sin_teta =  (xc[ind] - x0) / rc;
    va = gam * exp(-0.5 * (r * r) / rc2);

    /* Density and pressure */
    ro = roinf * exp(-coef * exp(-(r * r) / rc2));
    pp = roinf * 287.053 * tinf + (ro - roinf) * (c0 * c0);

    /* Conservative variables */
    u1[ind] = ro;
    u2[ind] = ro * uinf + ro * cos_teta * va;
    u3[ind] = ro * sin_teta * va;
    u4[ind] = 0.0;
    u5[ind] = pp / 0.4 + 0.5 * (u2[ind] * u2[ind] + u3[ind] * u3[ind]) / ro;
  }
}

// ============================================================================
/* Init by a lamb vortex */
// ============================================================================
PyObject* K_INITIATOR::initLamb(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initLamb: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initLamb: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0; //says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }

  PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, api);
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, npts, *cn, eltType, 0, api, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  // init with lamb
  k6lamb(x0, y0, Gam, MInf, npts,
         f->begin(posx), f->begin(posy), f->begin(posz),
         f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
         f2->begin(posrow), f2->begin(posroe));

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

// ===========================================================================
/* Init by a Visbal vortex */
// ============================================================================
PyObject* K_INITIATOR::initVisbal(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initVisbal: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initVisbal: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0;//says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }
  
    PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, api);
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, npts, *cn, eltType, 0, api, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  // init with visbal
  k6visbal(x0, y0, Gam, MInf, npts,
           f->begin(posx), f->begin(posy), f->begin(posz),
           f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
           f2->begin(posrow), f2->begin(posroe));

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

// ============================================================================
/* Init by a scully vortex */
// ============================================================================
PyObject* K_INITIATOR::initScully(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, coreRadius, MInf;
  E_Int model=0;
  if (!PYPARSETUPLE_(args, O_ TRR_ RRR_ I_,
                     &array, &x0, &y0, &Gam, &coreRadius, &MInf, &model))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initScully: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initScully: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0;//says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }

    PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, api);
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, npts, *cn, eltType, 0, api, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  switch (model)
  {
    case 0:
      k6scully(x0, y0, Gam, coreRadius, MInf, npts,
               f->begin(posx), f->begin(posy), f->begin(posz),
               f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
               f2->begin(posrow), f2->begin(posroe));
    default:
      k6scully2(x0, y0, Gam, coreRadius, MInf, npts,
                 f->begin(posx), f->begin(posy), f->begin(posz),
                 f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
                 f2->begin(posrow), f2->begin(posroe));
  }

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

// ===========================================================================
/* Init by a Yee vortex */
// ============================================================================
PyObject* K_INITIATOR::initYee(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, Minf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &Minf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initYee: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initYee: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0;//says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;  
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }
  
    PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, api);
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, npts, *cn, eltType, 0, api, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  k6yee(x0, y0, Gam, Minf, npts,
        f->begin(posx), f->begin(posy), f->begin(posz),
        f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
        f2->begin(posrow), f2->begin(posroe));

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

// ============================================================================
/* Init by Wissocq's vortex */
// ============================================================================
PyObject* K_INITIATOR::initWissocq(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float x0, y0, Gam, MInf;
  
  if (!PYPARSETUPLE_(args, O_ TRR_ RR_, &array, &x0, &y0, &Gam, &MInf))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res !=2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initWissocq: invalid type of array.");
    return NULL;
  }
  char varStringOut[K_ARRAY::VARSTRINGLENGTH];

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "initWissocq: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int defVars = 0; //says if output array is built using default variables or not
  E_Int neqTot = f->getNfld();
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int posro = K_ARRAY::isDensityPresent(varString);
  E_Int posrou = K_ARRAY::isMomentumXPresent(varString);
  E_Int posrov = K_ARRAY::isMomentumYPresent(varString);
  E_Int posrow = K_ARRAY::isMomentumZPresent(varString);
  E_Int posroe = K_ARRAY::isEnergyStagnationDensityPresent(varString);
  if (posro == -1 || posrou == -1 || posrov == -1 || posrow == -1 || posroe == -1)
  {
    defVars = 1;
    neqTot = 8;
    strcpy(varStringOut, "x,y,z,ro,rou,rov,row,roE");
    posro = 4; posrou = 5; posrov = 6; posrow = 7; posroe = 8;
  }
  else
  {
    posro++; posrou++; posrov++; posrow++; posroe++;
    strcpy(varStringOut, varString);
  }

  PyObject* tpl;
  if (res == 1)
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, im, jm, km, api);
  else
    tpl = K_ARRAY::buildArray3(neqTot, varStringOut, npts, *cn, eltType, 0, api, true);

  FldArrayF* f2;  
  K_ARRAY::getFromArray3(tpl, f2);
  
  // copy previous data
  E_Float* fp; E_Float* f2p;
  if (defVars == 1)
  {
    fp = f->begin(posx); f2p = f2->begin(1);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i]; 
    fp = f->begin(posy); f2p = f2->begin(2);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    fp = f->begin(posz); f2p = f2->begin(3);
    for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
  }
  else
  {
    for (E_Int n = 1; n <= neqTot; n++)
    {
      fp = f->begin(n); f2p = f2->begin(n);
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }

  k6wissocq(x0, y0, Gam, MInf, npts,
            f->begin(posx), f->begin(posy), f->begin(posz),
            f2->begin(posro), f2->begin(posrou), f2->begin(posrov), 
            f2->begin(posrow), f2->begin(posroe));

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}
