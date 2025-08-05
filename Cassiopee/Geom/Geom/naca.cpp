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
# include "geom.h"

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6naca2_(const E_Float& e, E_Int& N,
                E_Float* x, E_Float* y, E_Float* z);
  void k6nacas4g_(E_Int& im, E_Int& ip, E_Int& it, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
  void k6nacas5g_(E_Int& il, E_Int& ip, E_Int& iq, E_Int& it, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
  void k6nacas4m_(E_Int& im, E_Int& ip, E_Int& ith, E_Int& it, E_Int& ii, E_Int& sharpte, 
                  E_Int& npt, E_Float* x, E_Float* y, E_Float* z, E_Float* xl);
}

/* 
  naca avec fermeture de Van Rouzaud
*/
void k6naca1(E_Float e, E_Int npt,
             E_Float* x, E_Float* y, E_Float* z)
{
  E_Int np, nps2, n;
  E_Float pi, usnps2, xpn, ypn, en, epais, si, is;

  epais = e / 20.0;
  np = npt;
  nps2 = (np - 1) / 2;
  usnps2 = 1.0 / nps2;
  pi = 3.1415926;

  /* Naca donne analytiquement (d'apres Van Rouzaud) */
  for (n = 1; n <= npt; n++)
  {
    is = (E_Float)(n - 2) / nps2;
    is = 2 * is - 1;
    si = is;
    en = n - nps2 - 1;
    en = sqrt(1.008930411365) * sin(0.5 * pi * en * usnps2);
    xpn = en * en;
    ypn = si * epais * (0.2969 * std::abs(en) - 0.126 * xpn - 0.3516 * xpn * xpn
          + 0.2843 * xpn * xpn * xpn - 0.1015 * xpn * xpn * xpn * xpn);
    x[n - 1] = xpn / 1.008930411365;
    y[n - 1] = ypn / 1.008930411365;
    z[n - 1] = 0.0;
  }

  ypn = 0.5 * (y[0] + y[npt - 1]);
  y[0] = ypn;
  y[npt - 1] = ypn;
}

/*
  Naca mais avec fermeture a la leballeur
*/
void k6naca2(E_Float e, E_Int& npt,
             E_Float* x, E_Float* y, E_Float* z)
{
  E_Int np, nps2, n, no, nr;
  E_Float pi, usnps2, xpn, ypn, en, epais, si, is;
  E_Float hbfs2, bdabf;

  epais = e / 20.0;
  np = npt;
  nps2 = (np - 1) / 2;
  usnps2 = 1.0 / nps2;
  pi = 3.1415926;
  
  /* Naca donne analytiquement */
  nr = 1;
  for (n = 1; n <= npt; n++)
  {
    is = E_Float(n - 2) / E_Float(nps2);
    is = 2 * is - 1;
    si = E_Int(is);
    en = n - nps2 - 1;
    en = sqrt(1.008930411365) * sin(0.5 * pi * en * usnps2);
    xpn = en * en;
    ypn = 0.2969 * std::abs(en) - 0.126 * xpn - 0.3516 * (xpn * xpn)
          + 0.2843 * (xpn * xpn * xpn) - 0.1015 * (xpn * xpn * xpn * xpn);
    ypn = si * epais * ypn;
    if (xpn <= 1.0)
    {
      x[nr - 1] = xpn;
      y[nr - 1] = ypn;
      z[nr - 1] = 0.0;
      nr++;
    }
  }
  npt = nr - 1;
  
  /* Fermeture */
  ypn = 0.5 * (y[0] + y[npt - 1]);
  xpn = 0.5 * (x[0] + x[npt - 1]);
  hbfs2 = ypn - y[0];
  bdabf = 70.0;
  if (bdabf * ypn > epais)
  {
    bdabf = epais / hbfs2;
  }

  for (n = 1; n <= npt / 4; n++)
  {
    no = npt - n + 1;
    y[n - 1] = y[n - 1] + hbfs2 * exp(std::max(bdabf * (x[n - 1] - xpn), -30.0));
    y[no - 1] = y[no - 1] - hbfs2 * exp(std::max(bdabf * (x[no - 1] - xpn), -30.0));
  }

  x[0] = 1.0;
  x[npt - 1] = 1.0;
}

/*
  Naca Serie 4 airfoils by Dr Moens
*/
void k6nacas4g(E_Int im, E_Int ip, E_Int it,
               E_Int sharpte, E_Int npt,
               E_Float* x, E_Float* y, E_Float* z,
               E_Float* xl)
{
  E_Int i, shift;
  E_Float a5, coef, rap, xx, y_thick, yy;
  E_Float teta, y_camb;
  E_Float xm;  /* valeur max de la ligne de cambrure (en % corde) */
  E_Float xp;  /* position du maximum de cambrure (en 10 eme de corde) */
  E_Float xth;

  /* Passage digits -> grandeurs */
  xm = im * 0.010;
  xp = ip * 0.100;
  xth = it * 0.010;

  if (sharpte == 1)
  {
    shift = 1;
  }
  else
  {
    shift = 0;
  }

  // Evolution en X : loi imposee en sin**3
  coef = 0.92 * acos(-1.0) / 2.0;
  rap = 1.0 / pow(sin(coef), 3);
  xl[0] = 0.0;
  for (i = 2; i <= npt - 1; i++)
  {
    xl[i - 1] = rap * pow(sin(coef * (E_Float)(i - 1) / (E_Float)(npt - 1)), 3);
  }
  xl[npt - 1] = 1.0;

  /* *************************** */
  a5 = -0.1015;
  if (sharpte == 1) a5 = -0.1036;

  for (i = 1; i <= npt; i++)
  {
    xx = xl[i - 1];
    y_thick = 0.2969 * sqrt(xx) - 0.1260 * xx - 0.3516 * (xx * xx)
            + 0.2843 * (xx * xx * xx) + a5 * (xx * xx * xx * xx);
    y_thick = y_thick * xth / 0.20;

    if (xx < xp)
    {
      if (xp != 0.0)
      {
        y_camb = (2.0 * xp - xx) * xm * xx / (xp * xp);
        teta = (xp - xx) * 2.0 * xm / (xp * xp);
      }
      else
      {
        y_camb = 0.0;
        teta = 0.0;
      }
      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
    else
    {
      if (xp != 1.0)
      {
        y_camb = (1.0 + xx - 2.0 * xp) * xm * (1.0 - xx) / ((1.0 - xp) * (1.0 - xp));
        teta = (xp - xx) * 2.0 * xm / ((1.0 - xp) * (1.0 - xp));
      }
      else
      {
        y_camb = 0.0;
        teta = 0.0;
      }
      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
  }

  /* Pour le bord de fuite, on impose x/c=1 et l'epaisseur */
  yy = (0.1036 + a5) * xth / 0.20;

  x[npt - 1 + 1 - shift] = 1.0;
  y[npt - 1 + 1 - shift] = yy;
  z[npt - 1 + 1 - shift] = 0.0;
  x[npt - 1] = 1.0;
  y[npt - 1] = -yy;
  z[npt - 1] = 0.0;
}

/*
  Naca Serie 5 airfoils by Dr Moens
*/
void k6nacas5g(E_Int il, E_Int ip, E_Int iq, E_Int it,
               E_Int sharpte, E_Int npt,
               E_Float* x, E_Float* y, E_Float* z,
               E_Float* xl)
{
  E_Float xm[5], xk1[5], xf[5];
  E_Float teta, y_camb, y_thick;
  E_Float a5, coef, yy, rap, xt, xk2sk1, xx, xxk1;
  E_Int i, shift;

  if (sharpte == 1)
  {
    shift = 1;
  }
  else
  {
    shift = 0;
  }

  /* Evolution en X : loi imposee en sin**3 */
  xt = it * 0.010;

  /* Comme les digits sont des entiers, autant rentrer les valeurs tabulees ... */
  xf[0] = 0.05;
  xf[1] = 0.10;
  xf[2] = 0.15;
  xf[3] = 0.20;
  xf[4] = 0.25;

  if (iq == 0)
  {
    xm[0] = 0.0580; xm[1] = 0.1260; xm[2] = 0.2025; xm[3] = 0.2900; xm[4] = 0.3910;
    xk1[0] = 361.4; xk1[1] = 51.65; xk1[2] = 15.65;  xk1[3] = 6.643; xk1[4] = 3.230;
  }
  else
  {
    xm[0] = 0.0580; xm[1] = 0.1300; xm[2] = 0.2170; xm[3] = 0.3180; xm[4] = 0.4410;
    xk1[0] = 361.4; xk1[1] = 51.99; xk1[2] = 15.793; xk1[3] = 6.520; xk1[4] = 3.191;
  }

  /* Calcul coefficient k1 et k2/k1 (Mason) */
  xk2sk1 = (3.0 * pow(xm[ip - 1] - xf[ip - 1], 2) - pow(xm[ip - 1], 3))
         / pow(1.0 - xm[ip - 1], 3);

  /* adaptation du coef k1 selon le premier digit (Cz d'adaptation) */
  xxk1 = xk1[ip - 1] * (E_Float)il / 2.0;

  /* Evolution en X : loi imposee en sin**3 */
  coef = 0.92 * acos(-1.0) / 2.0;
  rap = 1.0 / pow(sin(coef), 3);
  xl[0] = 0.0;
  for (i = 2; i <= npt - 1; i++)
  {
    xl[i - 1] = rap * pow(sin(coef * (E_Float)(i - 1) / (E_Float)(npt - 1)), 3);
  }
  xl[npt - 1] = 1.0;

  a5 = -0.1015;
  if (sharpte == 1) a5 = -0.1036;

  for (i = 1; i <= npt; i++)
  {
    xx = xl[i - 1];
    y_thick = 0.2969 * sqrt(xx) - 0.1260 * xx - 0.3516 * (xx * xx)
            + 0.2843 * (xx * xx * xx) + a5 * (xx * xx * xx * xx);
    y_thick = y_thick * xt / 0.20;

    if (xx < xm[ip - 1])
    {
      if (iq == 0)
      {
        y_camb = pow(xx, 3) - 3.0 * xm[ip - 1] * (xx * xx)
               + xx * (3.0 - xm[ip - 1]) * pow(xm[ip - 1], 2);
        y_camb = y_camb * xxk1 / 6.0;

        teta = 3.0 * (xx * xx) - 6.0 * xm[ip - 1] * xx
             + (3.0 - xm[ip - 1]) * pow(xm[ip - 1], 2);
        teta = teta * xxk1 / 6.0;
      }
      else
      {
        y_camb = (pow(xx - xm[ip - 1], 3) - xk2sk1 * xx * pow(1.0 - xm[ip - 1], 3)
               - xx * pow(xm[ip - 1], 3) + pow(xm[ip - 1], 3)) * xxk1 / 6.0;

        teta = xxk1 / 6.0 * (3.0 * pow(xx - xm[ip - 1], 2)
             - xk2sk1 * pow(1.0 - xm[ip - 1], 3) - pow(xm[ip - 1], 3));
      }

      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
    else
    {
      if (iq == 0)
      {
        y_camb = xxk1 / 6.0 * pow(xm[ip - 1], 3) * (1.0 - xx);
        teta = -xxk1 / 6.0 * pow(xm[ip - 1], 3);
      }
      else
      {
        y_camb = (xk2sk1 * pow(xx - xm[ip - 1], 3)
               - xk2sk1 * xx * pow(1.0 - xm[ip - 1], 3)
               - xx * pow(xm[ip - 1], 3) + pow(xm[ip - 1], 3)) * xxk1 / 6.0;

        teta = xxk1 / 6.0 * (3.0 * xk2sk1 * pow(xx - xm[ip - 1], 2)
             - xk2sk1 * pow(1.0 - xm[ip - 1], 3) - pow(xm[ip - 1], 3));
      }

      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
  }

  /* Pour le bord de fuite, on impose x/c=1 et l'epaisseur */
  yy = (0.1036 + a5) * xt / 0.20;

  x[npt - 1 + 1 - shift] = 1.0;
  y[npt - 1 + 1 - shift] = yy;
  z[npt - 1 + 1 - shift] = 0.0;
  x[npt - 1] = 1.0;
  y[npt - 1] = -yy;
  z[npt - 1] = 0.0;
}

/*
  Naca Serie 4 modified airfoils by Dr Moens
*/
void k6nacas4m(E_Int im, E_Int ip, E_Int ith,
               E_Int it, E_Int ii,
               E_Int sharpte, E_Int npt,
               E_Float* x, E_Float* y, E_Float* z,
               E_Float* xl)
{
  E_Float teta, y_camb;
  E_Float xm, xp, xth, xt;
  E_Float coef, rap, d0, d1, d2, d3;
  E_Float a0, a1, a2, a3;
  E_Float yy, xx, xkhile, y_thick, rho1;
  E_Int i, shift;

  if (sharpte == 1)
  {
    shift = 1;
  }
  else
  {
    shift = 0;
  }

  xm = im * 0.010;
  xp = ip * 0.100;
  xth = ith * 0.010;
  xt = it * 0.10;

  /* Evolution en X : loi imposee en sin**3 */
  coef = 0.92 * acos(-1.0) / 2.0;
  rap = 1.0 / pow(sin(coef), 3);
  xl[0] = 0.0;
  for (i = 2; i <= npt - 1; i++)
  {
    xl[i - 1] = rap * pow(sin(coef * (E_Float)(i - 1) / (E_Float)(npt - 1)), 3);
  }
  xl[npt - 1] = 1.0;

  /* Serie 4 Modifiee : uniquement la loi d'epaisseur (cambrure identique à la Serie 4) */
  /* Coefficients (from Mason) */
  d0 = 0.002;
  if (sharpte == 1) d0 = 0.0;

  d1 = (2.24 - 5.42 * xt + 12.3 * (xt * xt)) / (10.0 * (1.0 - 0.878 * xt));
  d2 = (0.294 - 2.0 * (1.0 - xt) * d1) / pow(1.0 - xt, 2);
  d3 = (-0.196 + (1.0 - xt) * d1) / pow(1.0 - xt, 3);

  xkhile = (E_Float)ii / 6.0;
  if (ii == 9) xkhile = 10.3933;

  a0 = 0.296904 * xkhile;
  rho1 = 0.2 * pow(1.0 - xt, 2) / (0.588 - 2.0 * d1 * (1.0 - xt));
  a1 = 0.3 / xt - 15.0 * a0 / (8.0 * sqrt(xt)) - xt / (10.0 * rho1);
  a2 = -0.3 / (xt * xt) + 5.0 * a0 / (4.0 * pow(xt, 1.5)) + 1.0 / (5.0 * rho1);
  a3 = 0.1 / (xt * xt * xt) - 0.375 * a0 / (xt * xt * sqrt(xt))
     - 1.0 / (10.0 * rho1 * xt);

  /* Generation profil */
  for (i = 1; i <= npt; i++)
  {
    xx = xl[i - 1];
    if (xx < xt)
    {
      y_thick = a0 * sqrt(xx) + a1 * xx + a2 * (xx * xx) + a3 * (xx * xx * xx);
    }
    else
    {
      y_thick = d0 + d1 * (1.0 - xx) + d2 * pow(1.0 - xx, 2)
              + d3 * pow(1.0 - xx, 3);
    }
    y_thick = y_thick * xth / 0.20;

    if (xx < xp)
    {
      if (xp != 0.0)
      {
        y_camb = (2.0 * xp - xx) * xm * xx / (xp * xp);
        teta = (xp - xx) * 2.0 * xm / (xp * xp);
      }
      else
      {
        y_camb = 0.0;
        teta = 0.0;
      }
      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
    else
    {
      teta = 0.0;
      if (xp != 1.0)
      {
        y_camb = (1.0 + xx - 2.0 * xp) * xm * (1.0 - xx) / pow(1.0 - xp, 2);
        teta = (xp - xx) * 2.0 * xm / pow(1.0 - xp, 2);
      }
      else
      {
        y_camb = 0.0;
        teta = 0.0;
      }
      teta = atan(teta);
      x[2 * npt - i - shift] = xx - y_thick * sin(teta);
      y[2 * npt - i - shift] = y_camb + y_thick * cos(teta);
      z[2 * npt - i - shift] = 0.0;
      x[i - 1] = xx + y_thick * sin(teta);
      y[i - 1] = y_camb - y_thick * cos(teta);
      z[i - 1] = 0.0;
    }
  }

  yy = d0 * xth / 0.20;

  x[npt - 1 + 1 - shift] = 1.0;
  y[npt - 1 + 1 - shift] = yy;
  z[npt - 1 + 1 - shift] = 0.0;
  x[npt - 1] = 1.0;
  y[npt - 1] = -yy;
  z[npt - 1] = 0.0;
}

// ============================================================================
/* Create a naca profile line of N points */
// ============================================================================
PyObject* K_GEOM::nacaMesh(PyObject* self, PyObject* args)
{
  E_Int N; E_Float e;
  E_Int im, ip, it, ith, iq, sharpte;
  if (!PYPARSETUPLE_(args, R_ IIII_ III_,
                    &e, &N,
                    &im, &ip, &it, &ith, &iq, &sharpte))
  {
    return NULL;
  }

  PyObject* tpl = NULL;
  if (im == -1) // input avec epaisseur
  {
    // Data check
    if (e < 0)
    {
      PyErr_SetString(PyExc_ValueError, "naca: thickness must be positive.");
      return NULL;
    }

    if ((N/2)*2-N == 0)
    {
      printf("Warning: naca: number of points must be odd.\n");
      printf("Warning: naca: number of points set to " SF_D_ ".\n", N+1);
      N = N+1;
    }
  
    E_Int n = E_Int(N);
    FldArrayF coord(N, 3);
    coord.setAllValuesAtNull();
    k6naca2_(e, n, coord.begin(1), coord.begin(2), coord.begin(3));
    coord.reAllocMat(n, 3);
    tpl = K_ARRAY::buildArray(coord, "x,y,z", n, 1, 1);
  }
  else // input avec digits
  {
    if (sharpte == 1 && (N/2)*2-N == 0)
    {
      // sharp : 2N+1
      printf("Warning: naca: number of points must be odd.\n");
      printf("Warning: naca: number of points set to " SF_D_ ".\n", N+1);
      N = N+1;
    }
    if (sharpte == 0 && (N/2)*2-N == 1)
    {
      // pas sharp : 2N
      printf("Warning: naca: number of points must be even.\n");
      printf("Warning: naca: number of points set to " SF_D_ ".\n", N+1);
      N = N+1;
    }

    E_Int npt = N/2; 
    if (sharpte == 1) npt++;

    FldArrayF xl(npt);
    FldArrayF coord(2*npt, 3);
    coord.setAllValuesAtNull();

    // determination de la serie
    if (im > -0.5 && ip > -0.5 && ith > -0.5 && it > -0.5 && iq > -0.5)
    {
      // iq used as ii
      k6nacas4m_(im, ip, ith, it, iq, sharpte, 
                npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    else if (im > -0.5 && ip > -0.5 && iq > -0.5 && it > -0.5)
    {
      // im used as il
      k6nacas5g_(im, ip, iq, it, sharpte,
                npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    else if (im > -0.5 && ip > -0.5 && it > -0.5)
    {
      k6nacas4g_(im, ip, it, sharpte,
                npt, coord.begin(1), coord.begin(2), coord.begin(3), xl.begin());
    }
    coord.reAllocMat(N, 3);
    tpl = K_ARRAY::buildArray(coord, "x,y,z", N, 1, 1);
  }

  return tpl;
}
