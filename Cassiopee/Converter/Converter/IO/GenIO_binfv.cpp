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

// fieldview format (not usable)

// ecriture_bin_fv : ecrit un champ structuré au format binaire fieldview
// entrées: 
//    - nom_fichier : le nom générique du fichier de sortie (il sera suivi de .xyz.bin_fv pour le maillage, .aero.bin_fv pour le champ aéro,
//		    .bin_fv.nam pour le nom des variables aéro)
//    - titre : un titre
//    - nom_variable : le tableau du nom des variables
//    - endi, endj, endk: dimensions du champ aéro
//    - nombre_variable: nombre de variables dans le champ aéro initial (sans les variables composées)
//    - var: champ à écrire, rangé dans un tableau 1D par i puis j puis k puis numero_variable croissants
// 
// Remarque: on écrit toujours en big endian
//
# include <stdio.h>
# include <stdlib.h>

//=============================================================================
void ecriture_bin_fv(string nom_fichier, string title, 
                     string* nom_variable, long endi, long endj, long endk, long nvar, double* var)
{
   bool petit_indien=false; // booleen: la machine est little-endian
   if (MachineLE()==1)
   {
      petit_indien=true;      
   }
   
   string nom_fichier_base=nom_fichier;
   long nvar_xyz=0;   
   long nvar_aero=0;  
   
   // On doit découper le champ en 2 parties: maillage + champ aero
   // On doit créer 2 tableaux qui correspondent à cette découpe
   for (long nv=0;nv<nvar;nv++)
   {
      if (nom_variable[nv]=="x" || nom_variable[nv]=="y" || nom_variable[nv]=="z")
      {
         nvar_xyz=nvar_xyz+1;
      }
      else
      {
         nvar_aero=nvar_aero+1;
      }
   }
   
   double* maillage;
   double* champ_aero;
   
   maillage=new double[endi*endj*endk*3];
   // on initialise le maillage à zero (ça peut servir si le maillage est 1D ou 2D, la/les coordonnée(s) inutile(s) seront alors à 0)
   for (long i=0;i<endi*endj*endk*3;i++)
   {
      maillage[i]=0.0;	      
   }
   
   champ_aero=new double[endi*endj*endk*nvar_aero];
   
   nvar_xyz=0;
   nvar_aero=0;
   
   for (long nv=0;nv<nvar;nv++)
   {
      if (nom_variable[nv]=="x" || nom_variable[nv]=="y" || nom_variable[nv]=="z")
      {
         for (long i=0;i<endi*endj*endk;i++)
	 {
	    maillage[i+endi*endj*endk*nvar_xyz]=var[i+endi*endj*endk*nv];	    
	 }
	 nvar_xyz=nvar_xyz+1;
      }
      else
      {
         for (long i=0;i<endi*endj*endk;i++)
	 {
	    champ_aero[i+endi*endj*endk*nvar_aero]=var[i+endi*endj*endk*nv];	    
	 }
	 nvar_aero=nvar_aero+1;
      }     
   }
   
   if (nvar_aero>0)
   {
      FILE* pFile;
      long si = sizeof(long);
      long sf = sizeof(double);
      int delim;
      long ib;
      nom_fichier=nom_fichier_base+"aero.bin_fv";
      cout << "   ecriture -- " << nom_fichier << endl;      
      pFile = fopen(nom_fichier.c_str(), "w");
      if (pFile == NULL)
      {
         printf("ecriture_bin_fv : echec lors de l'ouverture de %s \n",nom_fichier.c_str());
         exit(0);
      } 
   
      // taille 1 entier
      delim=si;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);   
      // nombre de domaine
      ib=1;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // taille 1 entier
      delim=si;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);   
      
      // taille 4 entiers
      delim=si*4;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      // endi
      ib=endi;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // endj
      ib=endj;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // endk
      ib=endk;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);      
      // nvar_aero
      ib=nvar_aero;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);      
      // taille 4 entiers
      delim=si*4;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      
      // taille champ aero
      delim=sf*endi*endj*endk*nvar_aero;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      float element;
      for (long nv=0;nv<nvar_aero;nv++)
      {
         for (long k=0;k<endk;k++)
         {
	    for (long j=0;j<endj;j++)
            {
	       for (long i=0;i<endi;i++)
               {
	          element=champ_aero[i+j*endi+k*endi*endj+nv*(endi*endj*endk)];
		  if (petit_indien)
        	  {element=DBE(element);}
		  fwrite(&element, sf, 1, pFile);
	       }
	    }
	 }
      }
           
      // taille champ aero
      delim=sf*endi*endj*endk*nvar_aero;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      fclose(pFile);
      delete champ_aero;
      
      nom_fichier=nom_fichier_base+"bin_fv.nam";
      cout << "   ecriture -- " << nom_fichier << endl;      
      pFile = fopen(nom_fichier.c_str(), "w");
      for (long nv=0;nv<nvar;nv++)
      {
         if (nom_variable[nv]!="x" && nom_variable[nv]!="y" && nom_variable[nv]!="z")
         { 
	    fprintf (pFile , "%s\n", nom_variable[nv].c_str());
	 }
      }
      fclose(pFile);
      
   }// fin ecriture champ aero
   
   
   if (nvar_xyz>0)
   {
      FILE* pFile;
      long si = sizeof(int);
      long sf = sizeof(float);
      int delim;
      long ib;
      nom_fichier=nom_fichier_base+"xyz.bin_fv";
      cout << "   ecriture -- " << nom_fichier << endl;
      pFile = fopen(nom_fichier.c_str(), "w");
      if (pFile == NULL)
      {
         printf("ecriture_bin_fv : echec lors de l'ouverture de %s \n",nom_fichier.c_str());
         exit(0);
      } 
   
      // taille 1 entier
      delim=si;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);   
      // nombre de domaine = 1
      ib=1;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // taille 1 entier
      delim=si;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);   
      
      // taille 3 entiers
      delim=si*3;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      // endi
      ib=endi;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // endj
      ib=endj;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);
      // endk
      ib=endk;
      if (petit_indien)
      {ib=LBE(ib);}
      fwrite(&ib, si, 1, pFile);      
      // taille 3 entiers
      delim=si*3;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      
      // taille maillage
      delim=sf*endi*endj*endk*3;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      float element;
      for (long nv=0;nv<3;nv++)
      {
         for (long k=0;k<endk;k++)
         {
	    for (long j=0;j<endj;j++)
            {
	       for (long i=0;i<endi;i++)
               {
	          element=maillage[i+j*endi+k*endi*endj+nv*(endi*endj*endk)];
		  if (petit_indien)
        	  {element=DBE(element);}
		  fwrite(&element, sf, 1, pFile);
	       }
	    }
	 }
      }
           
      // taille maillage
      delim=sf*endi*endj*endk*3;
      if (petit_indien)
      {delim=IBE(delim);}
      fwrite(&delim, si, 1, pFile);
      
      fclose(pFile);
      delete maillage;
   }
   
}

//=============================================================================

