// ============================================================================
// THIS FILE IS COPYRIGHTED - SEE Kernel/COPYRIGHT.txt
// ============================================================================
// File   : Distributor/distributor.h
// SVN    : $Rev:$ $Date$
// Cksum  :
// ============================================================================
// ===========================================================================
// Project : Mesh distributor for Chimera.
// ===========================================================================
#ifndef _DISTRIBUTOR_DISTRIBUTOR_H_
# define _DISTRIBUTOR_DISTRIBUTOR_H_

# include "Python.h"

# define E_DOUBLEREAL
//# define E_DOUBLEINT
# define _DEF_USE_ISO_
# define _DEF_USE_USING_
# define __STD std
# define _ISO_LIST_
# define _ISO_MAP_
# define _ISO_VECTOR_

# include "Spl/Api/SplImpl.h"

// ===========================================================================
// @Name SplDistributor
// @Memo
/* @Text

Synopsis
> Purpose :
>   Distribute original blocks amongst computing nodes.
*/
// ===========================================================================
class SplDistributor : public SplImpl
{
  public:
    ///+ 1- Constructors and destructor.
    SplDistributor();
    ///-
    ///+ 2- Accessors
    virtual const list<DesBlock*> getSplittedBlocks()  const;
    virtual const list<DesMesh*>  getSplittedMeshes()  const;
    virtual const list<DesMask*>  getSplittedMasks()   const;
    virtual const list<DesBoundary*> getSplittedBnds() const;
    ///-
    ///+ 3- Split and other methods.
    virtual void split();
    virtual void printSplitInfo();
    virtual E_Int getIdOfOriginalBlock( E_Int idSpltBlock );
    virtual SplImpl* clone() const;
    ///-
  private:
    E_Float _meanPointsPerProc;
};

#endif
