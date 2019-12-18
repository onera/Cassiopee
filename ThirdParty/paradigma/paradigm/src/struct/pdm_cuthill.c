/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_part.h"
#include "pdm_cuthill.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/

/**
 *
 * \brief Quick sort
 *
 * \param [inout]   a     Array to sort
 * \param [in]      l     First element
 * \param [in]      r     Last  element
 *
 */

static void
_quickSort_int
(
 int a[],
 int l,
 int r
)
{
  if (l < r) {
    int j = r+1;
    int t;
    int pivot = a[l];
    int i = l;

    while(1) {
      do ++i; while (a[i] <= pivot && i < r);
      do --j; while (a[j] > pivot);
      if (i >= j) break;

      t    = a[i];
      a[i] = a[j];
      a[j] = t;

    }
    t    = a[l];
    a[l] = a[j];
    a[j] = t;

    _quickSort_int(a, l  , j-1);
    _quickSort_int(a, j+1,   r);
  }
}


/**
 *
 * \brief Builds dual graph face cell connectivity
 *
 * \param [inout] part_ini                 Part object - fine mesh partition
 *
 * \param [inout] cellCellIdxCompressed    Array of indexes of the dual graph
 * \param [inout] cellCellCompressed       Dual graph
 *
 */

static void
_dual_graph_firstrank
(
  _part_t        *part_ini,
  int           **cellCellIdxCompressed,
  int           **cellCellCompressed
)
{
  //cellCellN: array of counters of the numbers of connectivities
  //cellCell: dual graph to be built
  //cellCellIdx: array of indexes of the dual graph (same as cellFaceIdx)

  int *cellCellN = (int *) malloc(part_ini->nCell * sizeof(int));
  for (int i = 0; i < part_ini->nCell; i++) {
    cellCellN[i] = 0;
  }

  int *cellCell = (int *) malloc(part_ini->cellFaceIdx[part_ini->nCell] * sizeof(int));
  for (int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
    cellCell[i] = -1;
  }

  int *cellCellIdx = (int *) malloc((part_ini->nCell + 1) * sizeof(int));
  for(int i = 0; i < part_ini->nCell + 1; i++) {
    cellCellIdx[i] = part_ini->cellFaceIdx[i];
  }

  for (int i = 0; i < part_ini->nFace; i++) {
    int iCell1 = PDM_ABS (part_ini->faceCell[2*i    ]) - 1;
    int iCell2 = PDM_ABS (part_ini->faceCell[2*i + 1]) - 1;
    //Only the non-boundary faces are stored
    if (iCell2 > 0) {
      int idx1 = cellCellIdx[iCell1] + cellCellN[iCell1];
      cellCell[idx1] = iCell2 + 1;
      cellCellN[iCell1] += 1;

      int idx2 = cellCellIdx[iCell2] + cellCellN[iCell2];
      cellCell[idx2] = iCell1 + 1;
      cellCellN[iCell2] += 1;
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cellCellN after looping over cellFace: ");
    for(int i = 0; i < part_ini->nCell; i++) {
      PDM_printf(" %d ", cellCellN[i]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCell after looping over cellFace: ");
    for(int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
      PDM_printf(" %d ", cellCell[i]);
    }
    PDM_printf("\n");
  }

  //cellCellIdx is rebuilt
  *cellCellIdxCompressed = malloc((part_ini->nCell + 1) * sizeof(int));

  (*cellCellIdxCompressed)[0] = 0;
  for(int i = 0; i < part_ini->nCell; i++) {
    (*cellCellIdxCompressed)[i + 1] = (*cellCellIdxCompressed)[i] + cellCellN[i];
  }

  //We compress the dual graph since cellCellIdx was built from cellFaceIdx
  //We have then nFace elements in cellCell whereas it needs to be composed of nCell elements

  //    PDM_printf("(*cellCellIdxCompressed)[part_ini->nCell] : %d \n", (*cellCellIdxCompressed)[part_ini->nCell]);
  *cellCellCompressed = malloc((*cellCellIdxCompressed)[part_ini->nCell] * sizeof(int));

  int cpt_cellCellCompressed = 0;
  for(int i = 0; i < part_ini->cellFaceIdx[part_ini->nCell]; i++) {
    //        PDM_printf("I am testing a value for the %d time! \n", i);

    //We have an information to store when a neighboring cell exists
    if(cellCell[i] > -1){
      //            PDM_printf("I am storing a value for the %d time! \n", i);
      //We add a -1 to have the graph vertices numbered from 0 to n (C numbering)
      (*cellCellCompressed)[cpt_cellCellCompressed++] = cellCell[i] - 1;
      //            PDM_printf("Valeur stockee : %d \n ", (*cellCellCompressed)[cpt_cellCellCompressed - 1]);
    }
  }

  if( 0 == 1) {
    PDM_printf("Content of cellCellCompressed after compression and renumbering: ");
    for(int i = 0; i < (*cellCellIdxCompressed)[part_ini->nCell]; i++) {
      PDM_printf(" %d ", (*cellCellCompressed)[i]);
    }
    PDM_printf("\n");
  }

  /* Free temporary arrays*/

  free(cellCellN);
  free(cellCell);
  free(cellCellIdx);

  //Remove duplicate cells of the dual graph
  //We use the following scheme:
  //We loop over the indexes for the whole array to subdivide it into subarrays
  //We sort locally each subarray (determined thanks to cellCellIdxCompressed)
  //We loop over each subarray
  //We store the first value of each subarray anyway
  //We store each non-duplicated value and increment the writing index
  //We update the index array at each iteration

  int idx_write = 0;
  int tabIdxTemp = 0;

  for (int i = 0; i < part_ini->nCell; i++) {
    _quickSort_int((*cellCellCompressed), tabIdxTemp, (*cellCellIdxCompressed)[i + 1] - 1);

    int last_value = -1;

    for (int j = tabIdxTemp; j < (*cellCellIdxCompressed)[i + 1]; j++) {
      //We need to have a local index (between 0 and nFace)
      //If the value is different from the previous one (higher than is the same as different since the array is sorted)

      if(last_value != (*cellCellCompressed)[j]) {
        (*cellCellCompressed)[idx_write++] = (*cellCellCompressed)[j];
        last_value = (*cellCellCompressed)[j];
      }
    }

    if (0 == 1) {
      PDM_printf("\n Contenu de cellCellCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < (*cellCellIdxCompressed)[part_ini->nCell]; i1++) {
        PDM_printf(" %d ", (*cellCellCompressed)[i1]);
      }
      PDM_printf("\n");
    }

    tabIdxTemp = (*cellCellIdxCompressed)[i + 1];
    (*cellCellIdxCompressed)[i + 1] = idx_write;

    if (0 == 1) {
      PDM_printf("\n Contenu de cellCellIdxCompressed apres reecriture: \n");
      for(int i1 = 0; i1 < part_ini->nCell + 1; i1++) {
        PDM_printf(" %d ", (*cellCellIdxCompressed)[i1]);
      }
      PDM_printf("\n");
    }
  }

  if (0 == 1) {
    PDM_printf("Content of cellCellIdxCompressed after compression: ");
    for(int i1 = 0; i1 < part_ini->nCell + 1; i1++) {
      PDM_printf(" %d ", (*cellCellIdxCompressed)[i1]);
    }
    PDM_printf("\n");

    PDM_printf("Content of cellCellCompressed after compression: ");
    for(int i1 = 0; i1 < (*cellCellIdxCompressed)[part_ini->nCell]; i1++) {
      PDM_printf(" %d ", (*cellCellCompressed)[i1]);
    }
    PDM_printf("\n");
  }

  //We reallocate the memory in case of duplicated values removed
  //The new array size is idx_write (stored in (*cellCellIdxCompressed)[part_ini->nCell])
  *cellCellCompressed = realloc(*cellCellCompressed,
                                (*cellCellIdxCompressed)[part_ini->nCell] * sizeof(int));

}

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Compute Bandwidth of a graph
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     ADJ_BANDWIDTH,      The bandwidth of the adjacency matrix.
 */

static int
_adj_bandwidth
(
int node_num,
int adj_row[],
int adj[]
)
{
  int band_hi;
  int band_lo;
  int i,j,col;
  int value;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j-1] - 1;
      band_lo = _MAX( band_lo, i - col );
      band_hi = _MAX( band_hi, col - i );
    }
  }

  value = band_lo + 1 + band_hi;

  return value;

}

/**
 *
 * \brief TODOUX
 *
 */

static void
_level_set
(
int  root,
int  adj_row[],
int  adj[],
int  mask[],
int *level_num,
int  level_row[],
int  level[]
)
{
  int i;
  int iccsze;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  mask[root-1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
  /*
   *  LBEGIN is the pointer to the beginning of the current level, and
   *  LVLEND points to the end of this level.
   */
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num-1] = lbegin;
   /*
    *   Generate the next level by finding all the masked neighbors of nodes
    *   in the current level.
    */
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = level[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          iccsze = iccsze + 1;
          level[iccsze-1] = nbr;
          mask[nbr-1] = 0;
        }
      }
    }
    /*
     *   Compute the current level width (the number of nodes encountered.)
     *   If it is positive, generate the next level.
     */
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
    {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
  /**  Reset MASK to 1 for the nodes in the level structure. **/
  for ( i = 0; i < iccsze; i++ )
  {
    mask[level[i]-1] = 1;
  }

  return;
}

/**
 *
 * \brief ROOT_FIND finds a pseudo-peripheral node.
 *
 */

static void
_root_find
(
int *root,
int adj_row[],
int adj[],
int mask[],
int *level_num,
int level_row[],
int level[]
)
{
  int iccsze;
  int j;
  int jstrt;
  int k;
  int kstop;
  int kstrt;
  int level_num2;
  int mindeg;
  int nabor;
  int ndeg;
  int node;

  /** Determine the level structure rooted at ROOT. **/

  _level_set ( *root, adj_row, adj, mask, level_num,
              level_row, level);

    /** Count the number of nodes in this level structure. **/
  iccsze = level_row[*level_num] - 1;

  /*  Extreme case:
   *    A complete graph has a level set of only a single level.
   *    Every node is equally good (or bad).
   */

  if ( *level_num == 1 )
  {
    return;
  }
  /*   Extreme case:
   *     A "line graph" 0--0--0--0--0 has every node in its only level.
   *     By chance, we've stumbled on the ideal root.
   */
  if ( *level_num == iccsze )
  {
    return;
  }
  /*
   *   Pick any node from the last level that has minimum degree
   *   as the starting point to generate a new level set.
   */
  for ( ; ; )
  {
    mindeg = iccsze;

    jstrt = level_row[*level_num-1];
    *root = level[jstrt-1];

    if ( jstrt < iccsze )
    {
      for ( j = jstrt; j <= iccsze; j++ )
      {
        node = level[j-1];
        ndeg = 0;
        kstrt = adj_row[node-1];
        kstop = adj_row[node] - 1;

        for ( k = kstrt; k <= kstop; k++ )
        {
          nabor = adj[k-1];
          if ( 0 < mask[nabor-1] )
          {
            ndeg = ndeg + 1;
          }
        }

        if ( ndeg < mindeg )
        {
          *root = node;
          mindeg = ndeg;
        }
      }
    }

   /** Generate the rooted level structure associated with this node. **/
   _level_set( *root, adj_row, adj, mask, &level_num2,
                level_row, level);

   /** If the number of levels did not increase, accept the new ROOT. **/

    if ( level_num2 <= *level_num )
    {
      break;
    }

    *level_num = level_num2;
    /*
     *   In the unlikely case that ROOT is one endpoint of a line graph,
     *   we can exit now.
     */
    if ( iccsze <= *level_num )
    {
      break;
    }
  }

  return;

}

/**
 *
 * \brief Reverse an integer array
   \param [in]     The number of entries in the array.
   \param [in,out] int A(N), the array to be reversed.
 *
 */

static void
_i4vec_reverse
(
int n,
int a[]
)
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}

/**
 *
 * \brief TODOUX
 *
 */

static void
_degree
(
int root,
int adj_row[],
int adj[],
int mask[],
int deg[],
int *iccsze,
int ls[]
)
{
  int i;
  int ideg;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  /** The sign of ADJ_ROW(I) is used to indicate if node I has been considered. **/

  ls[0] = root;
  adj_row[root-1] = -adj_row[root-1];
  lvlend = 0;
  *iccsze = 1;
  /*
   *   LBEGIN is the pointer to the beginning of the current level, and
   *   LVLEND points to the end of this level.
   */
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
    /*
     *   Find the degrees of nodes in the current level,
     *   and at the same time, generate the next level.
     */
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = ls[i-1];
      jstrt = -adj_row[node-1];
      jstop = abs ( adj_row[node] ) - 1;
      ideg = 0;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          ideg = ideg + 1;

          if ( 0 <= adj_row[nbr-1] )
          {
            adj_row[nbr-1] = -adj_row[nbr-1];
            *iccsze = *iccsze + 1;
            ls[*iccsze-1] = nbr;
          }
        }
      }
      deg[node-1] = ideg;
    }

    /** Compute the current level width. **/

    lvsize = *iccsze - lvlend;

    /**  If the current level width is nonzero, generate another level. **/

    if ( lvsize == 0 )
    {
      break;
    }
  }

  /** Reset ADJ_ROW to its correct sign and return. **/

  for ( i = 0; i < *iccsze; i++ )
  {
    node = ls[i] - 1;
    adj_row[node] = -adj_row[node];
  }

  return;
}

/**
 *
 * \brief Compute reverse Cut-Hill Mac-Kee ordering
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_num[NODE_NUM+1] The number of adjacency entries
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     perm                The RCM ordering
 */

static void
_rcm
(
int root,
int adj_row[],
int adj[],
int mask[],
int perm[],
int *iccsze,
int node_num
)
{

  int fnbr;
  int i;
  int j;
  int jstop;
  int jstrt;
  int k;
  int l;
  int lbegin;
  int lnbr;
  int lperm;
  int lvlend;
  int nbr;
  int node;

  /** If node_num out of bounds, something is wrong. **/

  if ( node_num < 1 )
  {
    PDM_printf( "\n");
    PDM_printf( "RCM - Fatal error!\n");
    PDM_printf( "  Unacceptable input value of NODE_NUM = %d \n ", node_num);
    exit ( 1 );
  }

  /** If the root is out of bounds, something is wrong. **/

  if ( root < 1 || node_num < root )
  {
    PDM_printf("\n");
    PDM_printf("RCM - Fatal error!\n");
    PDM_printf("  Unacceptable input value of ROOT = %d \n ",root);
    PDM_printf("  Acceptable values are between 1 and %d inclusive \n",node_num);
    exit ( 1 );
  }

  /** Allocate memory for the degree array. **/
  int *deg = (int *) malloc( sizeof(int) * node_num); // deg = new int[node_num];

  /** Find the degrees of the nodes in the component specified by MASK and ROOT. **/

  _degree ( root, adj_row, adj, mask, deg, iccsze, perm);

  /** If the connected component size is less than 1, something is wrong. **/

  if ( *iccsze < 1 )
  {
    PDM_printf("\n");
    PDM_printf("RCM - Fatal error!\n");
    PDM_printf("  Connected component size ICCSZE returned from DEGREE as %d \n ", *iccsze);
    exit ( 1 );
  }

  /** Set the mask value for the root. **/

  mask[root-1] = 0;

  /** If the connected component is a singleton, there is no ordering necessary. **/

  if ( *iccsze == 1 )
  {
    free(deg);// delete [] deg;
    return;
  }
  /*   Carry out the reordering.
   *
   *   LBEGIN and LVLEND point to the beginning and
   *  the end of the current level respectively.
   */
  lvlend = 0;
  lnbr = 1;

  while ( lvlend < lnbr )
  {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for ( i = lbegin; i <= lvlend; i++ )
    {

      /** For each node in the current level...**/

      node = perm[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      /* Find the unnumbered neighbors of NODE.
       *   FNBR and LNBR point to the first and last neighbors
       *   of the current node in PERM.
       */
      fnbr = lnbr + 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          lnbr = lnbr + 1;
          mask[nbr-1] = 0;
          perm[lnbr-1] = nbr;
        }
      }

      /** If no neighbors, skip to next node in this level. **/

      if ( lnbr <= fnbr )
      {
        continue;
      }
     /*
      *   Sort the neighbors of NODE in increasing order by degree.
      *   Linear insertion is used.
      */
      k = fnbr;

      while ( k < lnbr )
      {
        l = k;
        k = k + 1;
        nbr = perm[k-1];

        while ( fnbr < l )
        {
          lperm = perm[l-1];

          if ( deg[lperm-1] <= deg[nbr-1] )
          {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }

  /*
  *   We now have the Cuthill-McKee ordering.
  *   Reverse it to get the Reverse Cuthill-McKee ordering.
  */
  _i4vec_reverse ( *iccsze, perm );

  /**  Free memory. **/
  free(deg);// delete [] deg;


  return;
}


/**
 *
 * \brief Compute reverse Cut-Hill Mac-Kee ordering
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_num[NODE_NUM+1] The number of adjacency entries
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     perm                The RCM ordering
 */

void
PDM_genrcm
(
int node_num,
int adj_row[],
int adj[],
int perm[]
)
{
  int i;
  int iccsze;
  int level_num;
  int num;
  int root;

  int *level_row  = (int *) malloc(sizeof(int) * node_num + 1); //level_row = new int[node_num+1];
  int *mask       = (int *) malloc(sizeof(int) * node_num    ); //mask = new int[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    mask[i] = 1;
  }

  num = 1;

  for ( i = 0; i < node_num; i++ )
  {
    /*
     * For each masked connected component...
     */
    // printf(" mask[%i] = %i \n", i, mask[i]);
    if ( mask[i] != 0 )
    {
      root = i + 1;
      /*
       *  Find a pseudo-peripheral node ROOT.  The level structure found by
       *   ROOT_FIND is stored starting at PERM(NUM).
       */
      _root_find(&root, adj_row, adj, mask, &level_num,
                level_row, perm+num-1);
      /*
       *   RCM orders the component using ROOT as the starting node.
       */
      _rcm(root, adj_row, adj, mask, perm+num-1, &iccsze, node_num );

      num = num + iccsze;

      /*
       *   We can stop once every node is in one of the connected components.
       */
      if ( node_num < num )
      {
        free(level_row); //delete [] level_row;
        free(mask);      //delete [] mask;
        return;
      }
    }
  }

  free(level_row); //delete [] level_row;
  free(mask);      //delete [] mask;

  return;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compute bandwidth of a mesh partition
 *
 * \param [in]  part ppart structure
 */

int
PDM_cuthill_checkbandwidth
(
 _part_t           *ppart
)
{
  int *dualGraphIdx = NULL;
  int *dualGraph    = NULL;

  _dual_graph_firstrank(ppart,
                        (int **) &dualGraphIdx,
                        (int **) &dualGraph);

  /** Offset Graph and Arr **/
  for (int i = 0; i < ppart->nCell; i++){
    dualGraphIdx[i] = dualGraphIdx[i]+1;
  }
  for (int i = 0; i < dualGraphIdx[ppart->nCell]; i++){
    dualGraph[i] = dualGraph[i]+1;
  }

  int dualBandWidth = _adj_bandwidth(dualGraphIdx[ppart->nCell], dualGraphIdx, dualGraph);

  /** Free memory **/
  free(dualGraphIdx);
  free(dualGraph);

  return dualBandWidth;
}

/**
 *
 * \brief Compute reverse CutHill-MecKee reordering
 *
 * \param [in]  part ppart structure
 */

void
PDM_cuthill_generate
(
 _part_t           *ppart,
 int               *perm
)
{

  /** Graph computation (in the new partition ) **/
  int *dualGraphIdx = NULL;
  int *dualGraph    = NULL;

  _dual_graph_firstrank(ppart,
                        (int **) &dualGraphIdx,
                        (int **) &dualGraph);

  /** Offset Graph and Arr **/
  for (int i = 0; i < ppart->nCell; i++){
    dualGraphIdx[i] = dualGraphIdx[i]+1;
  }
  for (int i = 0; i < dualGraphIdx[ppart->nCell]; i++){
    dualGraph[i] = dualGraph[i]+1;
  }

  /** Apply rcm to current Graph **/
  PDM_genrcm(ppart->nCell, dualGraphIdx, dualGraph, perm);

  /** Offset Permutation array **/
  for (int i = 0; i < ppart->nCell; i++)
  {
    perm[i] = perm[i]-1;
  }

  /** Verbose **/
  if (0 == 1) {
      PDM_printf("\n Contenu de perm : \n");
      for(int i = 0; i < ppart->nCell; i++) {
        PDM_printf(" %d ", perm[i]);
    }
    PDM_printf("\n");
  }

  /** Free memory **/
  free(dualGraphIdx);
  free(dualGraph);
}

#ifdef  __cplusplus
}
#endif
