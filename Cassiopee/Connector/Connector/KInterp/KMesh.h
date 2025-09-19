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
#ifndef _CONNECTOR_KMESH_H_
#define _CONNECTOR_KMESH_H_

# include "Def/DefTypes.h"
# include "Fld/FldArray.h"

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI
namespace K_KINTERP
{
// ============================================================================
/*
>   KMesh definit : 
>   -soit un maillage structure 3D 
>   -soit un maillage non structure de type TETRA               
>   
>   Dans le cas structure:
>    - (im,jm,km) : number of mesh planes in each direction
>    - Coordinates of points ( for construction, order of points is
>       (i-1)+ (j-1)*im + (k-1)*im*jm )
>    Pay attention :
>    For SetX, SetY, SetZ, SetPoint, getPos methods:
>    the needed integers are included in (1,im)x(1,jm)x(1,km).
>  
>   Dans le cas non structure :
>   - cn designe la connectivite. cn(.,.) demarre a 1 
>   C est la taille du tableau de connectivite qui permet de dire si le 
>   maillage est structure ou non
> 
*/
//=============================================================================
class KMesh
{
  public: 
    ///+ 1- Constructors / Destructor
    
    /** Empty constructor */
    KMesh();
  
    /** Maillage structure
        Uses sizes. Constructor is available for tests. */
    KMesh(E_Int im, E_Int jm, E_Int km);

    /** Maillage structure
        Uses sizes and FldArrayF. (all x then y then z. order
        of points ni+nj*npti+nk*npti*nptj) */
    KMesh(E_Int im, E_Int jm, E_Int km, const FldArrayF& coord);
    
    /** Maillage tetraedrique*/ 
    KMesh(FldArrayF& coord, FldArrayI& cn);
    
    /** Destructor */
    virtual ~KMesh();
    ///-
        
    ///+ 3- GET methods communes a un KMesh structure et non structure :
    
    /** Maillage structure ou non ? base sur la taille de _cn */
    E_Bool isStructured();

    /** Get x-coordinate of (i,j,k)  point (l = getPos(i,j,k)) */
    E_Float getX(E_Int l) const;
    
    /** Get y-coordinate of (i,j,k)  point (l = getPos(i,j,k)) */
    E_Float getY(E_Int l) const;
    
    /** Get z-coordinate of (i,j,k)  point (l = getPos(i,j,k)) */
    E_Float getZ(E_Int l) const;

    /** Get X vector */
    E_Float* getXVector();
    
    /** Get Y vector */
    E_Float* getYVector();
    
    /** Get Z vector */
    E_Float* getZVector();
    
    /** Get X vector */
    const E_Float* getXVector() const;
    
    /** Get Y vector */
    const E_Float* getYVector() const;
    
    /** Get Z vector */
    const E_Float* getZVector() const;
    
    /** Get read/write coordinate array */
    FldArrayF& getCoord();
    
    /** Get read coordinate array */
    const FldArrayF& getCoord() const;

    /** Return the cell volume */
    FldArrayF& getCellVol();

    /** Return the cell center coordinates */
    FldArrayF& getCellCenter();

    /** Find the bounding box of a mesh.
        @param xmax;ymax;zmax   coordinates of the upper-right point of 
                                bounding box
        @param xmin;ymin;zmin   coordinates of the lower-left point of 
                                bounding box
    */
    void boundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                     E_Float& xmin, E_Float& ymin, E_Float& zmin);
    
    /** Find the bounding box of a mesh in the absolute frame.
        @param m   rotation matrix of mesh relatively to absolute frame
        @param r0  deplacement of mesh origin in the absolute frame
    */
    void boundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                     E_Float& xmin, E_Float& ymin, E_Float& zmin,
                     const FldArrayF& m, const FldArrayF& r0,
                     const FldArrayF& xc0 ) const;

    ///+ 4- methodes specifiques a un KMesh structure
    void computeCellGreatLength(E_Float& dxmax, E_Float& dymax,
                                E_Float& dzmax);

    /** Build a BlkMesh from a BlkMesh.
        The original BlkMesh is supposed to be the nodes of a grid. The built
        BlkMesh contains the center of cells and points on the centers of 
        boundary interfaces.
        @param originalMesh Mesh (nodes) to be transformed to centers. */
    void createExtendedCenterMesh(const KMesh& origMesh);

    /** Build a BlkMesh from a BlkMesh.
        The original BlkMesh origMesh is supposed to be the nodes of the grid.
        The built BlkMesh is a duplication of origMesh on the left or right 
        sides following teta. The resulting mesh contains centers and
        points on the centers of boundary interface.
        @param origMesh   Mesh (nodes) to be transformed.
        @param axisVct    Vector of rotation
        @param axisPnt    Center of rotation
        @param teta       Angle of rotation (in radian)
    */
    void createDuplicatedExtendedPeriodMesh(const KMesh& origMesh,
                                            FldArrayF& axisVct,
                                            FldArrayF& axisPnt,
                                            E_Float theta);

    /** Find the cartesian elements bounding box (CEBB) for varying
        elements in Z.
        @param IN  ni, nj -> discretization number of points in x and y dir
        @param OUT xmin, ymin, zmin -> bounding box of mesh 
        @param OUT xmax, ymax, zmax -> bounding box of mesh
        @param OUT deltax -> discretization step in x
        @param OUT deltay -> discretization step in y
        @param OUT cartZmin -> min of cartesian element
        @param OUT cartZmax -> max of cartesian element.

        A cartesian element (i,j) is :
        (xmin + i*deltax, ymin + j*deltay, cartZmin(i,j) ) x
        (xmin + (i+1)*deltax, ymin + (j+1)*deltay, cartZmax(i,j) ).
    */
    void computeZCEBB(E_Int ni, E_Int nj,
                      E_Float& xmin, E_Float& ymin, E_Float& zmin,
                      E_Float& xmax, E_Float& ymax, E_Float& zmax,
                      E_Float& deltax, E_Float& deltay,
                      FldArrayF& cartZmin, FldArrayF& cartZmax);

    /** Find the cartesian elements bounding box (CEBB) for varying
        elements in Y.
        @param IN  ni, nk -> discretization number of points in x and z dir
        @param OUT xmin, ymin, zmin -> bounding box of mesh 
        @param OUT xmax, ymax, zmax -> bounding box of mesh
        @param OUT deltax -> discretization step in x
        @param OUT deltaz -> discretization step in y
        @param OUT cartYmin -> min of cartesian element
        @param OUT cartYmax -> max of cartesian element.

        A cartesian element (i,j) is :
        (xmin + i*deltax, cartYmin(i,j), zmin + k*deltaz ) x
        (xmin + (i+1)*deltax, cartYmax(i,j), zmax + (k+1)*deltaz ).
    */
    void computeYCEBB(E_Int ni, E_Int nk,
                      E_Float& xmin, E_Float& ymin, E_Float& zmin,
                      E_Float& xmax, E_Float& ymax, E_Float& zmax,
                      E_Float& deltax, E_Float& deltaz,
                      FldArrayF& cartYmin, FldArrayF& cartYmax);

    /** Find the cartesian elements bounding box (CEBB) for varying
        elements in X.
        @param IN  nj, nk -> discretization number of points in x and z dir
        @param OUT xmin, ymin, zmin -> bounding box of mesh 
        @param OUT xmax, ymax, zmax -> bounding box of mesh
        @param OUT deltay -> discretization step in x
        @param OUT deltaz -> discretization step in y
        @param OUT cartXmin -> min of cartesian element
        @param OUT cartXmax -> max of cartesian element.

        A cartesian element (i,j) is :
        (cartXmin(i,j), ymin + j*deltay, zmin + k*deltaz ) x
        (cartXmax(i,j), ymin + (j+1)*deltay, zmax + (k+1)*deltaz ).
    */
    void computeXCEBB(E_Int nj, E_Int nk,
                      E_Float& xmin, E_Float& ymin, E_Float& zmin,
                      E_Float& xmax, E_Float& ymax, E_Float& zmax,
                      E_Float& deltay, E_Float& deltaz,
                      FldArrayF& cartXmin, FldArrayF& cartXmax);

    /* Retourne la connectivite dans le cas non structure */
    FldArrayI& getConnectivity();

    /* Retourne le nombre de points ds le KMesh*/
    E_Int getNumberOfPts();

    /** Get im, number of mesh-planes in i direction */
    E_Int getIm() const;
    
    /** Get jm, number of mesh-planes in j direction */
    E_Int getJm() const;
    
    /** Get km, number of mesh-planes in k direction */
    E_Int getKm() const;

    /** Get (i-1)+ (j-1)*_im + (k-1)*_im*_jm position of (i,j,k) 
        point in arrays */
    E_Int getPos(E_Int i, E_Int j, E_Int k) const;


    ///+ 5- methodes specifiques a un KMesh non structure  

    ///-
  private:                       // methods

    KMesh(const KMesh&);
    KMesh& operator=(const KMesh&);

  private:                       // data
    
    E_Int    _im;                // number of mesh-planes in i direction
    E_Int    _jm;                // number of mesh-planes in j direction
    E_Int    _km;                // number of mesh-planes in k direction
    E_Int    _imjm;              // _im*_jm for optimization
    E_Bool _isStruct;         // KMesh struct ou non 
    FldArrayF _coord;            // coordinates array(nnodes,3)
    E_Int    _npts;              // taille de _coord 
    FldArrayF _cellVol;          // volume of cells 
    FldArrayF _bary;             // coordinates of  cell centers  
    FldArrayI _cn;               // non structure : connectivite tetra->sommets
};
}
//-----------------------------------------------------------------------------
inline E_Int K_KINTERP::KMesh::getNumberOfPts()
{
  return _npts;
}
inline E_Bool K_KINTERP::KMesh::isStructured()
{
  return _isStruct;
}
inline E_Int K_KINTERP::KMesh::getIm() const
{
  return _im;
}

inline E_Int K_KINTERP::KMesh::getJm() const 
{
  return _jm;
}

inline E_Int K_KINTERP::KMesh::getKm() const
{
  return _km;
}

inline E_Float* K_KINTERP::KMesh::getXVector()
{
  return _coord.begin();
}

inline E_Float* K_KINTERP::KMesh::getYVector()
{
  return _coord.begin() + _npts;
}

inline E_Float* K_KINTERP::KMesh::getZVector()
{
  return _coord.begin() + 2*_npts;
}

inline const E_Float* K_KINTERP::KMesh::getXVector() const
{
  return _coord.begin();
}

inline const E_Float* K_KINTERP::KMesh::getYVector() const
{
  return _coord.begin() + _npts;
}

inline const E_Float* K_KINTERP::KMesh::getZVector() const
{
  return _coord.begin() + 2*_npts;
}

inline FldArrayF& K_KINTERP::KMesh::getCoord()
{
  return _coord;
}

inline const FldArrayF& K_KINTERP::KMesh::getCoord() const
{
  return _coord;
}

inline E_Float K_KINTERP::KMesh::getX(E_Int l) const
{
  return _coord(l,1);
}

inline E_Float K_KINTERP::KMesh::getY(E_Int l) const 
{
  return _coord(l,2);
}

inline E_Float K_KINTERP::KMesh::getZ(E_Int l) const
{
  return _coord(l,3);
}


#undef FldArrayF
#undef FldArrayI

#endif
// ====================== Interp/KMesh.h === Last line ===================
