// 2. Une fois calculee la vitesse tangentielle au pt IBC
// on interpole lineairement la composante normale

ucible += alphasbeta*un;
vcible += alphasbeta*vn;
wcible += alphasbeta*wn;

/* Si modele de Spalart-Allmaras: traitement de la viscosite turbulente 
   Methode de Newton pour corriger la viscosite turbulente par fv1  
if (nvars == 6) 
{
  npass = 0;
  nutilde = K_FUNC::E_abs(nutcible);
# include "IBC/fnutilde.h"
  //f1v = fnutilde(nutilde, nutcible, roext, muext);
  while (K_FUNC::E_abs(f1v) > newtoneps && npass <= 30)
  {    
# include "IBC/fnutildeprime.h"
    //f1p = fnutildeprime(nutilde, nutcible, roext, muext);
    if (K_FUNC::E_abs(f1p) < newtonepsprime) break;
    nutilde = K_FUNC::E_abs(nutilde-f1v/f1p);
#   include "IBC/fnutilde.h"
    //f1v = fnutilde(nutilde, nutcible, roext, muext);    
    npass++;
  }
  if ( npass <= 30 ) nutcible = nutilde*signibc;
  varSAOut[indR] = nutcible; // nutilde au pt IBC
}
*/
