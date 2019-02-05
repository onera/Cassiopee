/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <float.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_morton.h"
#include "pdm_mpi.h"
#include "pdm_box.h"
#include "pdm_sort.h"
#include "pdm_box_priv.h"
#include "pdm_box_tree.h"
#include "pdm_dbbtree.h"
#include "pdm_dbbtree_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

/**
 * \brief  Initialize box_tree statistics
 *
 * \param [inout]  bts  pointer to box tree statistics structure
 *
 */

static void
_init_bt_statistics
(
_box_tree_stats_t  *bts
)
{
  size_t i;

  assert(bts != NULL);

  bts->dim = 0;

  for (i = 0; i < 3; i++) {
    bts->depth[i] = 0;
    bts->n_leaves[i] = 0;
    bts->n_boxes[i] = 0;
    bts->n_threshold_leaves[i] = 0;
    bts->n_leaf_boxes[i] = 0;
    bts->mem_used[i] = 0;
    bts->mem_required[i] = 0;
  }
}


/**
 * \brief Update box-tree statistics.
 *
 * For most fields, we replace previous values with the current ones.
 *
 * For memory required, we are interested in the maximum values over time
 * (i.e. algorthm steps); this is the case even for the minimal memory
 * required, we is thus the time maximum of the rank minimum.
 *
 * \param [inout]   bts   Pointer to box tree statistics structure
 * \param [inout]   bt    Pointer to box tree structure
 *
 */

static void
_update_bt_statistics
(
_box_tree_stats_t     *bts,
const PDM_box_tree_t  *bt
)
{
  int dim;
  size_t i;
  size_t mem_required[3];

  assert(bts != NULL);

  dim = PDM_box_tree_get_stats (bt,
                                bts->depth,
                                bts->n_leaves,
                                bts->n_boxes,
                                bts->n_threshold_leaves,
                                bts->n_leaf_boxes,
                                bts->mem_used,
                                mem_required);

  bts->dim = dim;

  for (i = 0; i < 3; i++)
    bts->mem_required[i] = _MAX(bts->mem_required[i], mem_required[i]);
}


/**
 * \brief Distribute bounding boxes over the ranks according to a Morton encoding
 * index. Try to get a well-balanced distribution and spatially coherent.
 *
 * \param [inout]  n     <-> pointer to neighborhood management structure
 * \param [inout]  boxes <-> box set to redistribute
 *
 */

static void
_redistribute_boxes
(
_PDM_dbbtree_t      *dbbt
)
{

  /* Sanity checks */

  assert (dbbt != NULL);

  PDM_box_tree_t  *coarse_tree = PDM_box_tree_create (dbbt->maxTreeDepthCoarse,
                                                      dbbt->maxBoxesLeafCoarse,
                                                      dbbt->maxBoxRatioCoarse);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (coarse_tree,
                          dbbt->boxes,
                          PDM_BOX_TREE_SYNC_LEVEL);

  _update_bt_statistics(&(dbbt->btsCoarse), coarse_tree);

  if (1 == 0) {
    PDM_printf ("-- dump stats\n");
    
    PDM_box_tree_dump_statistics(coarse_tree);
    
    PDM_printf ("-- fin dump stats\n");

    PDM_printf ("-- dump \n");
    
    PDM_box_tree_dump(coarse_tree);

    PDM_printf ("-- fin dump\n");
  }

  /*
   * Compute an index based on Morton encoding to ensure a good distribution
   * of bounding boxes among the ranks.
   */

  PDM_box_distrib_t  *distrib = PDM_box_tree_get_distrib (coarse_tree, dbbt->boxes);

  PDM_box_tree_destroy (&coarse_tree);

  if (1 == 0) {
    PDM_box_distrib_dump_statistics (distrib, dbbt->comm);
  }

  /* Define a new distribution of boxes according to the Morton
     encoding index */

  if (1 == 0) {
    PDM_printf("affichage 1\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 1\n");
  }

  PDM_box_set_redistribute (distrib, dbbt->boxes);

  if (1 == 0) {
    PDM_printf("affichage 2\n");
    PDM_box_set_dump( dbbt->boxes,1);
    PDM_printf("fin affichage 2\n");
  }

  /* Delete intermediate structures */

  PDM_box_distrib_destroy (&distrib);
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm          comm,
const int         dim
)
{
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) malloc(sizeof(_PDM_dbbtree_t));

  _dbbt->comm                 = comm;
  _dbbt->dim                  = dim;
  _dbbt->maxTreeDepth         = 30;
  _dbbt->maxBoxRatio          = 10.;
  _dbbt->maxBoxesLeaf         = 30;

  _dbbt->maxTreeDepthShared   = 30;
  _dbbt->maxBoxRatioShared    = 6;
  _dbbt->maxBoxesLeafShared   = 10;

  _dbbt->maxTreeDepthCoarse   = 30;
  _dbbt->maxBoxRatioCoarse    =  6.;
  _dbbt->maxBoxesLeafCoarse   = 30;

  _dbbt->rankBoxes            = NULL;
  _dbbt->btShared             = NULL;

  _dbbt->boxes                = NULL;
  _dbbt->btLoc                = NULL;

  _init_bt_statistics (&(_dbbt->btsShared));
  _init_bt_statistics (&(_dbbt->btsLoc));
  _init_bt_statistics (&(_dbbt->btsCoarse));

  return (PDM_dbbtree_t *) _dbbt;
}


/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
PDM_dbbtree_t     *dbbt
)
{
  if (dbbt != NULL) {
    _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;
    
    PDM_box_set_destroy (&_dbbt->rankBoxes);
    PDM_box_set_destroy (&_dbbt->boxes);
    
    PDM_box_tree_destroy (&_dbbt->btShared);
    PDM_box_tree_destroy (&_dbbt->btLoc);
    
    free (_dbbt);
  }

  return NULL;
}



/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  nPart    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to 
 * the tree location
 *
 */


PDM_box_set_t  *
PDM_dbbtree_boxes_set
(
PDM_dbbtree_t     *dbbt,
const int          nPart,
const int         *nElts,
const double     **extents,
const PDM_g_num_t **gNum
)
{
  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < nPart; i++) {
    nEltsProc += nElts[i];
  }

  PDM_printf ("nEltsProc : %d \n", nEltsProc);

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents      = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  for (int i = 0; i < nPart; i++) {

    for (int j = 0; j < nElts[i]; j++) {
      _boxGnum[idx++] = gNum[i][j];

      for (int k = 0; k < sExtents; k++) {
        _extents[idx1++] = extents[i][sExtents*j+k];
      }

      _initLocation[idx2++] = myRank;
      _initLocation[idx2++] = i;
      _initLocation[idx2++] = j;
    }
  }

  /*
   * Redistribute boxes of mesh A
   */ 

  _dbbt->boxes = PDM_box_set_create(3,
                                    0,  // No normalization to preserve initial extents
                                    0,  // No projection to preserve initial extents
                                    nEltsProc,
                                    _boxGnum,
                                    _extents,
                                    nPart,
                                    nElts,
                                    _initLocation,
                                    _dbbt->comm);

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  if (lComm > 1) {
    _redistribute_boxes(_dbbt);

    /*
     * Compute processus extents 
     */ 

    int nBoxes = PDM_box_set_get_size (_dbbt->boxes);
    const double *extents2 = PDM_box_set_get_extents (_dbbt->boxes);
    
    double gExtents[sExtents];
    for (int i = 0; i < _dbbt->dim; i++) {
      gExtents[i]   =  DBL_MAX;
      gExtents[_dbbt->dim+i] = -DBL_MAX;
    }
    
    
    for (int i = 0; i < nBoxes; i++) {
      for (int k1 = 0; k1 < _dbbt->dim; k1++) {
        gExtents[k1]   = _MIN (gExtents[k1], extents2[sExtents * i + k1]);
        gExtents[_dbbt->dim+k1] = _MAX (gExtents[_dbbt->dim+k1], extents2[sExtents * i
                                                                          + _dbbt->dim + k1]);          
      }
    }
    
    /*
     * Exchange extents to build shared boxes
     */ 
    
    int *allNBoxes = (int *) malloc (sizeof(int) * lComm); 
    PDM_MPI_Allgather (&nBoxes, 1, PDM_MPI_INT, 
                   allNBoxes, 1, PDM_MPI_INT, 
                   _dbbt->comm);
    
    int nUsedRank = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        nUsedRank += 1;
      }
    }
    
    
    double *allGExtents = (double *) malloc (sizeof(double) * sExtents * lComm); 
    PDM_MPI_Allgather (gExtents, sExtents, PDM__PDM_MPI_REAL, 
                   allGExtents, sExtents, PDM__PDM_MPI_REAL, 
                   _dbbt->comm);
    
    PDM_printf ("_extents shared :");
    idx1 = 0;
    for (int i = 0; i < lComm; i++) {
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf (" %12.5e", allGExtents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");


    PDM_g_num_t *numProc = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nUsedRank);
    idx = 0;
    for (int i = 0; i < lComm; i++) {
      if (allNBoxes[i] > 0) {
        numProc[idx] = i+1;
        for (int j = 0; j < sExtents; j++) {
          allGExtents[idx*sExtents + j] = allGExtents[i*sExtents + j];
        }
        idx += 1;
      }
    }

    allGExtents = (double *) realloc (allGExtents, sizeof(double) * sExtents * lComm);
    
    int *initLocationProc = (int *) malloc (sizeof(int) * 3 * nUsedRank);
    for (int i = 0; i < 3 * nUsedRank; i++) {
      initLocationProc[i] = 0;
    }

    _dbbt->rankBoxes = PDM_box_set_create(3,
                                          0,  // No normalization to preserve initial extents
                                          0,  // No projection to preserve initial extents
                                          nUsedRank,
                                          numProc,
                                          allGExtents,
                                          1,
                                          &nUsedRank,
                                          initLocationProc,
                                          _dbbt->comm);
  
    _dbbt->btShared = PDM_box_tree_create (_dbbt->maxTreeDepthShared,
                                           _dbbt->maxBoxesLeafShared,
                                           _dbbt->maxBoxRatioShared);

    /* Build a tree and associate boxes */

    PDM_box_tree_set_boxes (_dbbt->btShared,
                            _dbbt->rankBoxes,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    _update_bt_statistics(&(_dbbt->btsShared), _dbbt->btShared);

    free (allGExtents);
    free (numProc);
    free (initLocationProc);

  }

  /*
   * Build local bt
   */ 

  _dbbt->btLoc = PDM_box_tree_create (_dbbt->maxTreeDepth,
                                      _dbbt->maxBoxesLeaf,
                                      _dbbt->maxBoxRatio);

  /* Build a tree and associate boxes */

  PDM_box_tree_set_boxes (_dbbt->btLoc,
                          _dbbt->boxes,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  _update_bt_statistics(&(_dbbt->btsLoc), _dbbt->btLoc);

  return _dbbt->boxes;

}


/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  nPart    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according 
 * to the tree intersection
 *
 */

PDM_box_set_t  *
PDM_dbbtree_intersect_boxes_set
(
PDM_dbbtree_t    *dbbt,
const int         nPart,
const int        *nElts,
const double     **extents,
const PDM_g_num_t **gNum,
int              *box_index[],
int              *box_l_num[]
)
{

  assert (dbbt != NULL);
  _PDM_dbbtree_t *_dbbt = (_PDM_dbbtree_t *) dbbt;

  int myRank;
  PDM_MPI_Comm_rank (_dbbt->comm, &myRank);
  int lComm;
  PDM_MPI_Comm_size (_dbbt->comm, &lComm);

  const int nInfoLocation = 3;
  const int sExtents = _dbbt->dim * 2;

  int nEltsProc = 0;
  for (int i = 0; i < nPart; i++) {
    nEltsProc += nElts[i];
  }

  PDM_g_num_t *_boxGnum = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * nEltsProc);
  double *_extents     = (double *) malloc (sizeof(double) * nEltsProc * sExtents);
  int *_initLocation   = (int *) malloc (sizeof(int) * nEltsProc * nInfoLocation);

  int idx = 0;
  int idx1 = 0;
  int idx2 = 0;

  for (int i = 0; i < nPart; i++) {

    for (int j = 0; j < nElts[i]; j++) {
      _boxGnum[idx++] = gNum[i][j];

      for (int k = 0; k < sExtents; k++) {
        _extents[idx1++] = extents[i][sExtents*j+k];
      }

      _initLocation[idx2++] = myRank;
      _initLocation[idx2++] = i;
      _initLocation[idx2++] = j;
    }

  }

  if (1 == 0) {

    PDM_printf ("nEltsProc : %d", nEltsProc);
    
    PDM_printf ("_boxGnum :");
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" "PDM_FMT_G_NUM, _boxGnum[i]);
    }
    PDM_printf ("\n");
    
    PDM_printf ("_extents :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf (" %12.5e", _extents[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
    
    PDM_printf ("_initLocation :");
    idx1 = 0;
    for (int i = 0; i < nEltsProc; i++) {
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf (" %d", _initLocation[idx1++]);
      PDM_printf ("\n");
    }
    PDM_printf ("\n");
  }

  PDM_box_set_t  *boxes = PDM_box_set_create (3,
                                              0,  // No normalization to preserve initial extents
                                              0,  // No projection to preserve initial extents
                                              nEltsProc,
                                              _boxGnum,
                                              _extents,
                                              nPart,
                                              nElts,
                                              _initLocation,
                                              _dbbt->comm);

  free (_boxGnum);
  free (_extents);
  free (_initLocation);

  /*
   * Intersection boxes whith shared tree
   */ 

  if (_dbbt->btShared != NULL) {

    PDM_box_tree_get_boxes_intersects (_dbbt->btShared,
                                       boxes,
                                       box_index,
                                       box_l_num);
    
    int nUsedRank = PDM_box_set_get_size (_dbbt->rankBoxes);
    const PDM_g_num_t *usedRanks = PDM_box_set_get_g_num (_dbbt->rankBoxes);
    
    /*
     * Distribute boxes on intersection ranks
     */
    
    PDM_printf ("box_l_num_shared : ");
    for (int i = 0; i < nUsedRank; i++) {
      for (int j = (*box_index)[i]; j < (*box_index)[i+1]; j++) {
        PDM_printf (" %d", (*box_l_num)[j]);
      }
      PDM_printf ("\n");
    }
    
    PDM_box_distrib_t  *distrib = NULL;
    distrib = PDM_box_distrib_create (boxes->n_boxes,
                                      boxes->n_g_boxes,
                                      1, // Don't use in this case
                                      boxes->comm);
    
    for (int i = 0; i < lComm + 1; i++) {
      distrib->index[i] = 0;
    }

    for (int i = 0; i < nUsedRank; i++) {
      distrib->index[usedRanks[i]] = (*box_index)[i+1] - (*box_index)[i];
    }
    
    for (int i = 0; i < lComm; i++) {
      distrib->index[i+1] = distrib->index[i+1] + distrib->index[i];
    }
    
    distrib->list = (int *) malloc (sizeof(int) * distrib->index[lComm]);
    
    for (int i = 0; i < distrib->index[lComm]; i++) {
      distrib->list[i] = (*box_l_num)[i];
    }
    
    /*
     * Redistribute boxes on intersecting ranks
     */ 
    
    PDM_box_set_redistribute (distrib,
                              boxes);
  
    /*
     * Free
     */ 

    PDM_box_distrib_destroy (&distrib);

    free (*box_l_num);
    free (*box_index);
  
  }
  
  /*
   * Intersection boxes whith local tree
   */ 

  PDM_box_tree_get_boxes_intersects (_dbbt->btLoc,
                                     boxes,
                                     box_index,
                                     box_l_num);

  /*
   * Sort boxes and remove double boxes
   */ 

  int nBoxesA = _dbbt->boxes->n_boxes;

  int *newIndex = (int *) malloc (sizeof(int) * (nBoxesA + 1));
  for (int i = 0; i < nBoxesA + 1; i++) {
    newIndex[i] = 0;
  }

  int *_box_index = *box_index;
  int *_box_l_num = *box_l_num;

  idx = 0;
  for (int i = 0; i < nBoxesA; i++) {

    int *ideb = _box_l_num + _box_index[i];
    int length = _box_index[i+1] - _box_index[i];

    PDM_sort_int (ideb, NULL, length);

    int pre = -1;
    for (int k = 0; k < length; k++) {
      if (pre != ideb[k]) {
        (*box_l_num)[idx++] = ideb[k];
        newIndex[i+1] += 1;
        pre = ideb[k];
      }
    }
  }

  for (int i = 0; i < nBoxesA; i++) {
    newIndex[i+1] += newIndex[i];
  }

  free (*box_index);
  *box_index = newIndex;

  *box_l_num = (int *) realloc (*box_l_num, sizeof (int) * newIndex[nBoxesA]);

  return boxes;

}

#undef _MIN
#undef _MAX
 
#ifdef __cplusplus
}
#endif /* __cplusplus */
