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

// ============================================================================
/*
   Some numerical methods may require that a field defined on a Block 
   (BlkBaseBlock) need to be interpolated on other points (for example, in 
   the Chimera method). If such an 'interpolable' block is required, one has 
   to add a BlkInterpData to this Block and use its methods to get the 
   interpolation cells and interpolation coefficients.
   
   Instances of this class must be explicitely created by the user.
   
   This class is an interface to use various methods of interpolation :
   ICG: for general curvilinear grids, use cartesian preconditioning 
   (long to build, fast to use)
   CARTESIAN: for regular cartesian grids only (fastest technique 
   for these grids)
   ADT: for general curvilinear grids, use a tree preconditioning 
   (fast to build, long to use)
   CYL: for cylindrical grids (fastest technique for these grids),
   IRRCART: for irregular cartesian grids.

   ICG, ADT, CYL, IRRCART requires a in memory extended centers mesh 
   and so derive from BlkInterpExtCenters class.

   One problem when dealing with geometry are the epsilon values (used 
   for example to know if a point is inside a cell or nor, if a cell 
   intersects another...). These values (_EPS_DET, _EPS_TETRA, _EPS_GEOM) 
   are set in the constructor and defined further in this header.
*/
// ============================================================================
class BlkInterpData
{
  public:

    ///+ 0- enum ans struct
    /** Various interpolation types*/
    enum InterpolationType
    {
      O2CF, O3CF, O2ABC, O3ABC, O5ABC
    };
    /** Interpolation mesh type : defined in nodes or extended centers */
    enum InterpMeshType
    {
      NODE, EXT_CENTERS
    };
    /** Various interpolation methods that can be used.
     */
    enum InterpolationMethod
    {
      ICG, CARTESIAN, ADT, ADTP, CYL, IRRCART
    };
    
    struct InterpolationParameters
    {
      /* ICG paramaters */
      E_Int _icgDepth;
      E_Int _icgCut;
      E_Float _icgDimFactor;
      /* Periodicity by rotation (if any) */
      FldArrayF _axisPnt;
      FldArrayF _axisVect;
      E_Int _periodicDir;
      E_Float _theta;
      /* Periodicity by translation */
      FldArrayF _transVect;
    };

    ///-
    
    ///+ 1- Constructor / Destructor
    /** Constructor */
    BlkInterpData();
    
    /** Destructor. */
    virtual
    ~BlkInterpData();
    ///-

    ///+ 2- Get methods
    short getInterpolationCellStruct(E_Float x, E_Float y, E_Float z,
                                     E_Int& ic, E_Int& jc, E_Int& kc, 
                                     FldArrayF& cf, InterpMeshType interpMeshType,
                                     InterpolationType interpType=O2CF);
       
    /** Retourne les indices (ic,ic+1,..,jc, jc+1,kc, kc+1)  de la molecule
        d'interpolation et les coefficients ai, bj, ck correspondants.
        pour un point (x,y,z).  
        Return -1 if not found
                1 if found on the main domain (no periodic interpolation)
                2 if found on the +theta domain
                3 if found on the -theta domain    */
    short getInterpolationCellStruct(E_Float x, E_Float y, E_Float z,
                                     FldArrayI& indi, FldArrayF& cf,
                                     InterpMeshType interpMeshType,
                                     InterpolationType interpType=O2CF);

    /** CAS NON STRUCTURE
        IN : x,y,z coord du pt a interpoler
        OUT : noet : numero du tetraedre d interpolation ds _cn
        OUT : indi : elt donneur + marqueur en indi[0]
        OUT : cf : coefficients d interp. sur le tetraedre
        Retourne 0 si pas trouve
    */
    short getInterpolationCellUnstruct( E_Float x, E_Float y, E_Float z,
                                        E_Int& noet, 
                                        FldArrayI& indi, FldArrayF& cf);
    
    /** Get the interpolation cell (ic,jc,kc) and interpolation coefficients cf
        for a vector of points (x,y,z) given in coord. Only the part
        [istart : iend] of vectors is treated. Usefull for vectorization. */
    void getInterpolationCellStructv(
      FldArrayF& coord, E_Int istart, E_Int iend,
      FldArrayI& indi, FldArrayF& cf, FldArrayIS& found,
      FldArrayI& extrap,
      InterpMeshType interpMeshType,
      InterpolationType interpType=O2CF );
    /** same as before. Return 3 additionnal arrays containing (i,j,k)-indices from interpolation cells in 
     extended mesh*/
    void getInterpolationCellStructv(
      FldArrayF& coord, E_Int istart, E_Int iend,
      FldArrayI& indi, FldArrayI& icv, FldArrayI& jcv, FldArrayI& kcv, 
      FldArrayF& cf, FldArrayIS& found,
      FldArrayI& extrap,
      InterpMeshType interpMeshType,
      InterpolationType interpType=O2CF);
    
    /** Get extrapolation coefficients cf for point (x,y,z) inside cell 
	(ic, jc, kc) (in extended centers).
        
        If order is equal to zero, if at least an interpolable point 
	(cellNatureField=1) exists in this cell, return true and cf are 
	degenerated (order 0) extrapolation coefficients.
        
        If order is equal to one, try to find an interpolable tetrahedra. 
	Return true if OK, and cf are extrapolation cofficients in this 
	tetrahedra (order 1).

        If order is equal to two, we try first to find an interpolable 
	tetrahedra, if it is not possible, we switch to an order 0 formula.
        
        CellNatureField is information relative to the grid we are 
	interpolating from.

        @param x,y,z               IN  coordinates of point to extrapolate
        @param ic,jc,kc            IN  index of cell that furnishes the values
        @param cf                  OUT extrapolation coefficients
        @param cellNatureField     IN  cell nature field on the grid we are 
	                               extrapolating from (optional)
        @param order               IN  order of extrapolation formula (see up) 
	                               (optional)
        @param is,js,ks            OUT if extrapolation or degeneration cannot 
	                               be performed
                                       is,js,ks are the index of a neigbouring 
				       cell of ic,jc,kc nearer of (x,y,z) 
				       (optional).
     */
    virtual 
    short getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                                       E_Int ic, E_Int jc, E_Int kc,
                                       FldArrayF& cf, E_Float cfMax,
                                       BlkInterpData::InterpolationType interpType,
                                       BlkInterpData::InterpMeshType interpMeshType=EXT_CENTERS) = 0;

    virtual 
    short getExtrapolationCoeffForCell(E_Float x, E_Float y, E_Float z,
                                       E_Int ic, E_Int jc, E_Int kc,
                                       FldArrayF& cf,
                                       const FldArrayI& cellNatureField,
                                       E_Int testNature, E_Int order,
                                       E_Int& is, E_Int& js, E_Int& ks, E_Float cfMax,
                                       BlkInterpData::InterpolationType interpType,
                                       BlkInterpData::InterpMeshType interpMeshType=EXT_CENTERS) = 0;

    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) for structured meshes 
        (nearest boundary cell) */
    virtual 
    short getExtrapolationCell(E_Float x, E_Float y, E_Float z,
                               E_Int& ic, E_Int& jc, E_Int& kc,
                               FldArrayF& cf,
                               const FldArrayI& cellNatureField, 
                               E_Int testNature,
                               E_Float& test,
                               BlkInterpData::InterpolationType interpType,
                               BlkInterpData::InterpMeshType interpMeshType=EXT_CENTERS,
                               E_Float cfMax=30.) = 0;
    /** Find the extrapolation elt for point (x,y,z) for unstructured meshes */
    virtual 
    short getExtrapolationCellUnstr(E_Float x, E_Float y, E_Float z,
                                    E_Int& noelt,
                                    FldArrayF& cf,
                                    FldArrayI& indi,
                                    E_Int order,
                                    const FldArrayF& cellNatureField) = 0;

    /** Find the extrapolation cell (ic,jc,kc) for point (x,y,z) on structured meshes */
    virtual short getExtrapolationCellStruct(E_Float x, E_Float y, E_Float z,
                                             FldArrayI& indi,
                                             FldArrayF& cf,
                                             E_Int order, E_Float cfMax,
                                             const FldArrayF& cellNatureField,
                                             InterpolationType interpType) = 0;
   
    /** A partir du pt (ic,jc,kc) reconstruit le tableau indi des 
        indices en i,j,k des points de la molecule d interpolation
        Attention : les indices demarrent a 0 dans chq direction !
        Retourne 0 si extrapolation du bord, 1 sinon */
    virtual E_Int fromExtendedToStandardCenters(
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayI& indTab, 
      InterpolationType interpType) = 0 ;

    /** Meme chose mais sous forme de tableau : points a transformer entre
        istart et iend1-1 inclus. extrap = 0 pour un pt extrapole du bord */
    virtual void fromExtendedToStandardCentersv(
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayI& indTab, FldArrayI& extrap, 
      InterpolationType interpType) = 0 ;
    /* Calcul des indices i,i+1... de la molecule d' interpolation 
     le premier elt indi[0] correspond a l ordre d interpolation local 
     et n est pas calcule ici*/
    virtual void compStandardIndices(
      E_Int ic, E_Int jc, E_Int kc, 
      FldArrayI& indi, 
      InterpolationType interpType)= 0;
    /* Calcul des indices i,i+1... des molecules d' interpolation 
       comprises entre istart et iend-1 inclus
       le premier elt indi(.,1) correspond a l ordre d interpolation local 
       et n est pas calcule ici
    */
    virtual void compStandardIndicesv(
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc, 
      FldArrayI& indi, 
      InterpolationType interpType) = 0;

    /** 2nd order interpolation search using O2ABC interpolation type.  
        IN: x,y,z: point a interpole
        OUT: ic,jc,kc: cellule d'interpolation sur _mesh. */
    virtual 
    short searchInterpolationCellO2ABC(E_Float x, E_Float y, E_Float z, 
                                       E_Int& ic, E_Int& jc, E_Int& kc) = 0;
    /** 3rd order interpolation search using O3ABC interpolation type. 
        IN: x,y,z : point a interpole
        OUT: ic,jc,kc: cellule d'interpolation sur _mesh. */
    virtual short searchInterpolationCellO3ABC(
      E_Float x, E_Float y, E_Float z, E_Int& ic, E_Int& jc, E_Int& kc) = 0;
    /** 5th order interpolation search using O5ABC interpolation type. 
        IN: x,y,z: point a interpole
        OUT: ic,jc,kc: cellule d'interpolation sur _mesh. */
    virtual short searchInterpolationCellO5ABC(
      E_Float x, E_Float y, E_Float z, E_Int& ic, E_Int& jc, E_Int& kc) = 0;

    ///-

    //=========================================================================
  protected:
    virtual void 
    searchInterpolationCellByJumpv( FldArrayF& coord,
                                    E_Int istart, E_Int iend,
                                    FldArrayI& ic, FldArrayI& jc,FldArrayI& kc,
                                    FldArrayF& cf, FldArrayIS& found ) = 0;
    virtual short 
    searchInterpolationCellByJump(E_Float x, E_Float y, E_Float z,
                                  E_Int& ic, E_Int& jc, E_Int& kc,
                                  FldArrayF& cf) = 0;
    //=========================================================================
   /** Interpolation d ordre 2 non structuree */
    virtual short searchInterpolationCellUnstruct(
      E_Float x, E_Float y, E_Float z,
      E_Int& noelt, FldArrayI& indi, FldArrayF& cf) = 0;

    /** 2nd order interpolation search using O2CF interpolation type : 
        IN: x,y,z : point a interpole
        OUT: ic,jc,kc : cellule d'interpolation sur _mesh
        OUT: cf : coefficients d'interpolation (8)
    */
    virtual
    short searchInterpolationCellO2CF( E_Float x, E_Float y, E_Float z,
                                       E_Int& ic, E_Int& jc, E_Int& kc,
                                       FldArrayF& cf ) = 0;
    
    /** 3rd order interpolation search using O3CF interpolation type :
        IN: x,y,z : point a interpole
        OUT: ic,jc,kc : cellule d'interpolation sur _mesh
        OUT: cf : coefficients d'interpolation (27)
    */
    virtual
    short searchInterpolationCellO3CF( E_Float x, E_Float y, E_Float z,
                                       E_Int& ic, E_Int& jc, E_Int& kc,
                                       FldArrayF& cf ) = 0;

    /** O2CF interpolation search. Vectorized. */
    virtual void searchInterpolationCellO2CFv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayF& cf, FldArrayIS& found) = 0;
   
    /** O2ABC interpolation search. Vectorized. */
    virtual void searchInterpolationCellO2ABCv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found) = 0;
   
    /** O3CF interpolation search. Vectorized.
     O3CF is only done for InterpCart type. Return -1 for all the others 
    interpolation data types (with extended centers).*/
    virtual void searchInterpolationCellO3CFv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayF& cf, FldArrayIS& found) = 0;

    /** O3ABC interpolation search. Vectorized. */
    virtual void searchInterpolationCellO3ABCv( 
      FldArrayF& coord,
      E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found) = 0;

    /** O5ABC interpolation search. Vector. */
    virtual void searchInterpolationCellO5ABCv( 
      FldArrayF& coord, E_Int istart, E_Int iend,
      FldArrayI& ic, FldArrayI& jc, FldArrayI& kc,
      FldArrayIS& found) = 0;
  
    /** Calcule les polynomes de Lagrange pour un point (x,y,z) donné,
        avec les points de controle definis par le tableau indi
        si interpType = O2CF indi=[i,i+1,j,j+1,k,k+1], 
        si interpType = O3CF indi=[i,i+1,i+2,j,j+1,j+2,k,k+1,k+2] ...
        IN : x,y,z : point interpole
        IN : ic,jc,kc : indice de la cellule d'interpolation sur _mesh
        IN : interpType : type d'interpolation
        IN : interpMeshType
        OUT : cf : coefficients d'interpolation (alloue par l'appelant) 
        Retourne 0 si l interpolation est OK
        Retourne 1 si l ordre est degrade jusqu a 2 */
    virtual short compLagrangeCoefs(
      E_Float x, E_Float y, E_Float z,
      E_Int ic, E_Int jc, E_Int kc,
      FldArrayF& cf,
      InterpolationType interpType,
      InterpMeshType interpMeshType) = 0;
    /* Calcule les coefficients de Lagrange - version tableau 
       IN : coord : coordonnees des points contenant les points interpoles
       IN : istart : premier point interpole de coord
       IN : iend : dernier point interpole de coord
       IN : interpType : type d'interpolation
       IN : interpMeshType
       IN : icv, jcv, kcv : indices des cellules d'interpolation sur _mesh
       IN : found : vaut 1 si point interpolable de _mesh, -1 sinon
       OUT : cf : coefficients d interpolation 
       OUT : corr2 : tableau : vaut 1 si correction d ordre 2, 0 sinon
    */
    void compLagrangeCoefsv(FldArrayF& coord, E_Int istart, E_Int iend,
                            FldArrayI& icv, FldArrayI& jcv, FldArrayI& kcv,
                            FldArrayIS& found, FldArrayF& cf,
                            FldArrayIS& corr2,
                            InterpolationType interpType, 
                            InterpMeshType interpMeshType);

    /* Corrige a l'ordre 2 O2CF les coefs d interpolation et les indi
       Attention: les tableaux ne sont pas redimensionnes
       ex: ordre 5: cf de taille 15 toujours, mais les coefs 0 a 7 seult
       sont valables, idem pour indi de taille 16 a l ordre 5.
       Retourne le found du searchInterpolationCellO2CF 
       cad > 0 si point x,y,z interpolable a l ordre 2*/
    short correctInterpToO2CF( E_Float x, E_Float y, E_Float z, 
                               FldArrayI& indi, FldArrayF& cf,
                               InterpMeshType interpMeshType);
    
    /* Correction des coefficients d interpolation a l ordre 2 de type O2CF
       pour les points dont l'approximation n a pas marche
       IN: coord: coordonnees des points contenant les points interpoles
       IN: istart: premier point interpole de coord
       IN: iend: dernier point interpole de coord
       IN/OUT: indiv: indices des cellules d'interpolation sur _mesh
       IN/OUT: foundv: vaut 1 si point interpolable de _mesh, -1 sinon
       IN/OUT: cfv: coefficients d interpolation 
       IN: corrv: tableau : vaut 1 si correction d ordre 2 si OiABC, 0 sinon
       Attention : les tableaux ne sont pas redimensionnes  */
    void correctInterpToO2CFv(FldArrayF& coord, E_Int istart, E_Int iend,
                              FldArrayI& indiv,
                              FldArrayIS& foundv, FldArrayF& cfv,
                              FldArrayIS& corrv,
                              InterpMeshType interpMeshType);

    //=========================================================================
    /* 
       Find the interpolation coefficients for point (x,y,z) in an hexahedral
       cell. Hexahedra is described by its eight points coordinates:          
       (x0,y0,z0),(x1,y1,z1),... Numerotation of points in an hexahedral   
       cell is given by:                                                      
                        4   
                        o*************o 6
                     *  .          *  *
                  *     .       *     *
               *        .    *        *
          5 o*************o 7         *            K
            *           . *           *            .
            *           . *           *            .
            *           . *           *            .
            *         0 o.*...........o 2          ........ J
            *        .    *        *            .
            *     .       *     *            .
            *  .          *  *            I
          1 o*************o 3
    
    
          Interpolation coefficients are stored in the cf[8] vector and should 
          be used to interpolate a function F defined on the cell edges by:
          F(x,y,z)=cf[0]*F(xt[0],yt[0],zt[0])+cf[1]*F(xt[1],yt[1],zt[1])+...
          This function returns true if the point (x,y,z) is inside the cell,
          E_False otherwise.
    */
    E_Boolean coeffInterpHexa(E_Float x, E_Float y, E_Float z,
                              E_Float* xt, E_Float* yt, E_Float* zt,
                              FldArrayF& cf);
    E_Boolean coeffInterpHexav(E_Int index, E_Int index2,
                               E_Float x, E_Float y, E_Float z,
                               FldArrayF& xtv, FldArrayF& ytv, FldArrayF& ztv,
                               FldArrayF& cfv);
    void
    coeffInterpHexav(const E_Int istart, const E_Int nbI,
                     const FldArrayF& coord, 
                     FldArrayF& xtv, FldArrayF& ytv, FldArrayF& ztv,
                     FldArrayF& cfv, FldArrayIS& found, 
                     FldArrayI& indv, E_Int& c);
    /* 
       Find the interpolation coefficients for point (x,y,z) in a tetrahedra.
       Tetrahedra is described by its four points coordinates:                
       (xp,yp,zp),(xq,yq,zq),(xr,yr,zr),(xs,ys,zs),                            
       linear interpolation coefficient are: (xi,yi,zi)                       
       Point (x,y,z) is inside the tetrahedra if (xi>=0,yi>=0,zi>=0 and       
       xi+yi+zi<=1.                                                           
    */
    void coeffInterpTetra(E_Float x, E_Float y, E_Float z,
                          E_Float xp, E_Float yp, E_Float zp,
                          E_Float xq, E_Float yq, E_Float zq,
                          E_Float xr, E_Float yr, E_Float zr,
                          E_Float xs, E_Float ys, E_Float zs,
                          E_Float &xi, E_Float &yi, E_Float &zi);
    
    /* Find if the cell contains the point interpolated */
    E_Boolean getCellJump(E_Float x, E_Float y, E_Float z,
                          E_Float* xt, E_Float* yt, E_Float* zt,
                          E_Int& isomm,
                          E_Float& xi, E_Float& yi, E_Float& zi);

    /* Get the interp coeff in an hexa cell */
    E_Boolean getCoeffInterpHexa(E_Float x, E_Float y, E_Float z,
                                 E_Int isomm,
                                 E_Float xi, E_Float yi, E_Float zi, 
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& cf);

    const E_Float _EPS_DET;    // EPS for determinant in tetrahedra inversion
    const E_Float _EPS_TETRA;  // EPS for tetrahedra belonging test
    const E_Float _EPS_GEOM;   // EPS for Geometric intersection
};

// ======= Interp/BlkInterpData.h === Last line ==========================
